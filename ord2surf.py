#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2023, F. Hoffmann La-Roche Ltd.

import logging
import re
from collections import defaultdict

import pandas as pd
import typer
from ord_schema import message_helpers, validations
from ord_schema.proto import dataset_pb2
from rdkit.Chem import Descriptors, MolFromSmiles
from rich.progress import track
from typing_extensions import Annotated

from surf_utils.mappings import (
    const_cols,
    last_cols,
    mapping_analyses_ord,
    mapping_role_ord,
    mapping_stirring_ord,
    mapping_time,
    mapping_vol_mol_mass,
)

logging.basicConfig(format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

app = typer.Typer()


@app.command()
def surf2ord(
    input_file: Annotated[str, typer.Argument(help="name of the input file in ORD format; suffixes: .pbtxt or .pb")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in SURF format")],
    validate: Annotated[
        bool, typer.Option(help="whether the input ORD reactions should be validated for correctness first")
    ] = True,
):
    """Translate reaction data from the tabular SURF format into the protocol buffer format which is used in the open reaction database."""
    assert input_file.split(".")[-1] in [
        "pbtxt",
        "pb",
    ], "Unknown file extension for input_file! (Known: .pbtxt or .pb)"

    assert output_file.split(".")[-1] in [
        "txt",
        "tsv",
    ], "Unknown file extension for output_file! (accepted: .txt or .tsv)"

    logger.info(f"Reading dataset from file {input_file} ...")
    dataset = message_helpers.load_message(input_file, dataset_pb2.Dataset)

    if validate:
        logger.info("Validating dataset before processing...")
        validations.validate_datasets({input_file: dataset})

    surf = pd.DataFrame()

    for i, reaction in track(
        enumerate(dataset.reactions), total=len(dataset.reactions), description="Transforming ORD to SURF..."
    ):
        row = pd.Series()
        err = False

        if len(reaction.reaction_id):
            row["rxn_id"] = reaction.reaction_id
        else:
            logger.warning(f"No reaction identifier present in reaction {i+1}! Adding number {i+1} as rxn_id.")
            row["rxn_id"] = f"rxn_{i+1}"
            
        if len(reaction.identifiers):
            for id in reaction.identifiers:
                if id.type == 5:
                    row["rxn_type"] = id.value
                else:
                    row["rxn_name"] = id.value

        try:  # add origin (person) if present
            provenance = " ".join(
                [
                    reaction.provenance.record_created.person.name,
                    reaction.provenance.record_created.person.username,
                    reaction.provenance.record_created.person.email,
                    reaction.provenance.record_created.person.orcid,
                    reaction.provenance.record_created.person.organization,
                ]
            )
            row["provenance"] = re.sub(" {2,}", " ", provenance).strip()
        except AttributeError:
            row["provenance"] = "unknown"

        # add DOI / URL as source ID if present
        if len(reaction.provenance.doi):
            row["source_id"] = reaction.provenance.doi
            row["source_type"] = "publication"
        elif len(reaction.provenance.publication_url):
            row["source_id"] = reaction.provenance.publication_url
            row["source_type"] = "publication"

        # add temperature if present
        if reaction.conditions.temperature.setpoint.value:
            temp = reaction.conditions.temperature.setpoint.value
            if reaction.conditions.temperature.setpoint.units == 2:
                temp = (temp - 32) * 5 / 9  # F to °C
            row["temperature_deg_c"] = temp
        else:
            logger.warning(f"No temperature found in reaction {row['rxn_id']}!")

        try:  # add atmosphere if present
            row["atmosphere"] = reaction.conditions.pressure.atmosphere.__str__().strip().replace("type: ", "")
        except AttributeError:
            logger.warning(f"No atmosphere found in reaction {row['rxn_id']}!")

        try:  # add stirring info if present
            if reaction.conditions.stirring:
                row["stirring_shaking"] = mapping_stirring_ord[reaction.conditions.stirring.type]
        except AttributeError:
            pass

        try:  # add procedure if present
            row["procedure"] = reaction.notes.procedure_details
        except AttributeError:
            pass

        role_cnt = defaultdict(int)

        for inp in reaction.inputs:
            for cpd in reaction.inputs[inp].components:
                # process roles with smiles and names / cas
                if cpd.reaction_role == 5:  # workup
                    continue
                elif cpd.reaction_role == 0:  # unspecified
                    logger.warning(f"Not all reaction roles are specified in reaction {row.rxn_id}; skipping reaction")
                    err = True
                    continue
                role = mapping_role_ord[cpd.reaction_role]
                role_cnt[role] += 1
                for id in cpd.identifiers:
                    if id.value.lower() in ["unknown", "unspecified", "na", "n/a", "#n/a"]:
                        continue
                    if id.type == 2:  # SMILES
                        row[f"{role}_{role_cnt[role]}_smiles"] = id.value
                    elif id.type == 3:  # INCHI
                        row[f"{role}_{role_cnt[role]}_inchi"] = id.value
                    elif id.type == 6 and re.findall("^[0-9]+-[0-9]+-[0-9]+$", id.value):  # NAME
                        row[f"{role}_{role_cnt[role]}_cas"] = id.value
                    elif id.type == 7 and re.findall("^[0-9]+-[0-9]+-[0-9]+$", id.value):  # CAS
                        row[f"{role}_{role_cnt[role]}_cas"] = id.value
                    elif id.type != 4:
                        row[f"{role}_{role_cnt[role]}_name"] = id.value
                if f"{role}_{role_cnt[role]}_cas" not in row and f"{role}_{role_cnt[role]}_name" not in row:
                    row[f"{role}_{role_cnt[role]}_name"] = inp

                # process amounts and volumes
                if cpd.amount.volume.value:
                    row[f"{role}_{role_cnt[role]}_fraction"] = float(
                        cpd.amount.volume.value * mapping_vol_mol_mass[cpd.amount.volume.units]
                    )  # transform to ml
                elif cpd.amount.moles.value:
                    row[f"{role}_{role_cnt[role]}_eq"] = float(
                        cpd.amount.moles.value * mapping_vol_mol_mass[cpd.amount.moles.units]
                    )  # transform to mmol
                elif cpd.amount.mass.value and f"{role}_{role_cnt[role]}_smiles" in row:
                    try:  # mg to mmol
                        mol = MolFromSmiles(row[f"{role}_{role_cnt[role]}_smiles"])
                        row[f"{role}_{role_cnt[role]}_eq"] = (
                            cpd.amount.mass.value * mapping_vol_mol_mass[cpd.amount.mass.units]
                        ) / Descriptors.MolWt(mol)
                    except Exception:
                        logger.warning(
                            f"Could not calculate amount from given mass for {role}_{role_cnt[role]}_smiles in reaction {row.rxn_id}"
                        )
                else:
                    logger.warning(f"No amount found for {role}_{role_cnt[role]} in reaction {row.rxn_id}")

        # transform amounts to equvalents
        eq, vol = 0, 0  # mmol, ml
        smats = re.findall("startingmat_[0-9]_eq", " ".join(row.keys()))
        if smats:
            eq = min(row[smats])  # mmol
            for c in re.findall("[a-zA-Z0-9]+_[0-9]_eq", " ".join(row.keys())):
                row[c] = row[c] / eq

        # get solvent fractions
        solvs = re.findall("solvent_[0-9]_fraction", " ".join(row.keys()))
        if solvs:
            vol = sum(row[solvs])  # ml
            for solv in solvs:
                row[solv] = row[solv] / vol

        # get scale and concentration
        if eq:
            row["scale_mol"] = eq / 1000.0
            if vol:
                row["concentration_mol_l"] = eq / vol  # mmol/ml = mol/l

        # handle reaction outcomes
        for out in reaction.outcomes:
            try:  # add reaction time from outcome
                time = out.reaction_time.value * mapping_time[out.reaction_time.units]
                row["time_h"] = time
            except AttributeError:
                logger.warning(f"No reaction time found in reaction {row['rxn_id']}!")

            for prod in out.products:
                try:
                    role = mapping_role_ord[prod.reaction_role]
                except KeyError:
                    logger.warning(
                        f"Unknown role for product {prod.identifiers[0].value} in reaction {row['rxn_id']}; assumeing product!"
                    )
                    role = "product"
                role_cnt[role] += 1
                for id in prod.identifiers:
                    if id.type == 2:
                        row[f"{role}_{role_cnt[role]}_smiles"] = id.value
                    elif id.type == 3:
                        row[f"{role}_{role_cnt[role]}_inchi"] = id.value
                    elif id.type == 6 and re.findall("^[0-9]+-[0-9]+-[0-9]+$", id.value):
                        row[f"{role}_{role_cnt[role]}_cas"] = id.value
                    elif id.type == 7 and re.findall("^[0-9]+-[0-9]+-[0-9]+$", id.value):
                        row[f"{role}_{role_cnt[role]}_cas"] = id.value
                    elif id.type != 4:
                        row[f"{role}_{role_cnt[role]}_name"] = id.value

                # find yield
                for m in prod.measurements:
                    if m.type == 3:  # yield
                        row[f"{role}_{role_cnt[role]}_yield"] = m.percentage.value / 100
                        row[f"{role}_{role_cnt[role]}_yieldtype"] = mapping_analyses_ord[
                            out.analyses[m.analysis_key].type
                        ]
                    else:  # process nmr and ms analyses for product
                        ana = out.analyses[m.analysis_key]
                        if "NMR" in ana.__str__().upper():
                            try:
                                row[f"{role}_{role_cnt[role]}_nmr"] = ana.data["peaks"].string_value
                            except AttributeError:
                                pass
                        elif "MS" in ana.__str__().upper():
                            mass = ana.__str__().split("\n")[0].replace("type: ", "")
                            if "expected" in ana.data:
                                mass += f" Expected: {ana.data['expected'].float_value:.4f}"
                            if "found" in ana.data:
                                mass += f" Found: {ana.data['found'].float_value:.4f}"
                            row[f"{role}_{role_cnt[role]}_ms"] = mass.strip()

        if not err:
            surf = pd.concat((surf, row.to_frame().T))

    # save to tabular SURF format
    surf = surf[list(sorted(surf.columns.tolist()))]
    surf = surf.reindex(
        columns=[c for c in const_cols if c in surf.columns]
        + [c for c in surf.columns if (c not in const_cols + last_cols and "product" not in c)]
        + [c for c in surf.columns if "product" in c]
        + [c for c in last_cols if c in surf.columns]
    )
    surf = (
        surf.dropna(how="all", axis=1).convert_dtypes().round(8)
    )  # drop empty columns and round to max 8 decimal points
    surf.to_csv(output_file, sep="\t", index=False, float_format="%f")


if __name__ == "__main__":
    app()
