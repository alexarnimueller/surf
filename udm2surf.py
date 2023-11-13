#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2023, F. Hoffmann La-Roche Ltd.

import logging

import pandas as pd
import typer
import xmltodict
from bs4 import BeautifulSoup
from lxml import etree
from rdkit.Chem import MolFromMolBlock, MolToSmiles
from rich.progress import track
from typing_extensions import Annotated

from surf_utils.mappings import const_cols, last_cols

app = typer.Typer()

# TODO: get complete UDM example and finalize

logging.basicConfig(
    format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def validate_udm_file(xml_file):
    """Function to validate an xml file (.xml) against the UDM schema (.xsd)"""
    errors = []
    xmlschema = etree.XMLSchema(etree.parse("udm_schema/udm_6_0_0.xsd"))
    validation_rslt = xmlschema.validate(etree.parse(xml_file))
    if not validation_rslt:
        errors = [(e.path, e.message) for e in xmlschema.error_log]
    return errors


@app.command()
def udm2smiles(
    input_file: Annotated[
        str, typer.Argument(help="name of the input file in UDM format")
    ],
    output_file: Annotated[
        str,
        typer.Argument(
            help="name of the output file in SURF format; .tsv or .txt format"
        ),
    ],
    validate: Annotated[
        bool,
        typer.Option(
            help="Whether to validate the UDM file structure versus the UDM schema."
        ),
    ] = True,
):
    if validate:  # validate the XML file versus the UDM schema
        errors = validate_udm_file(input_file)
        if errors:
            logger.error(
                f"UDM schema validation found the following {len(errors)} errors:"
            )
            for e in errors:
                logger.error(e)
        else:
            logger.info("UDM validation successful")

    # read UDM data
    with open(input_file, "r") as f:
        raw = f.read()
    xml = BeautifulSoup(raw, "xml")

    # read UDM schema
    with open("udm_schema/udm_6_0_0.xsd", "r") as f:
        raw_xsd = f.read()
    xsd = BeautifulSoup(raw_xsd, "xml")

    udm = xml.find("UDM")
    ver = xml.find("UDM_VERSION")
    legal = xml.find_all("LEGAL")

    # read all citations
    citations = {}
    for cit in xml.find("CITATIONS").find_all("CITATION"):
        citations[cit.attrs["ID"]] = {}
        if cit.find("DOI"):
            citations[cit.attrs["ID"]].update({"DOI": cit.DOI.text})
        if cit.find("TYPE"):
            citations[cit.attrs["ID"]].update({"TYPE": cit.TYPE.text})

    # read all molecules
    molecules = xml.find("MOLECULES")
    mol_dict = {}
    for mol in track(
        molecules.find_all("MOLECULE"), description="Reading molecule structures..."
    ):
        m = MolFromMolBlock(mol.MOLSTRUCTURE.text)
        d = xmltodict.parse(mol.prettify())["MOLECULE"]
        if m:
            d.update({"MOLSTRUCTURE": MolToSmiles(m)})
            id = d.pop("@ID")
            mol_dict[id] = d
        else:
            logger.warning(f"Could not read Molblock of molecule {d['@ID']}")

    surf = pd.DataFrame()

    # read all reactions
    reactions = xml.find("REACTIONS")
    for reaction in track(
        reactions.find_all("REACTION"), description="Reading reactions..."
    ):
        row = pd.Series()
        rxn_id = reaction.attrs["ID"]
        row["rxn_id"] = rxn_id
        # Products
        try:
            prod_ids = [p.text for p in reaction.find_all("PRODUCT_ID")]
            products = {i: mol_dict[i] for i in prod_ids}
        except KeyError:
            logger.error(
                f"One of these product IDs {prod_ids} was not found in molecules! Skipping reaction {rxn_id}"
            )
            continue

        # Reactants -> Startingmaterials
        try:
            react_ids = [r.text for r in reaction.find_all("REACTANT_ID")]
            startmats = {i: mol_dict[i] for i in react_ids}
        except KeyError:
            logger.error(
                f"One of these reactant IDs {react_ids} was not found in molecules! Skipping reaction {rxn_id}"
            )
            continue

        # Variations (first update all reactant amounts to get the scale)
        for variation in reaction.find_all("VARIATION"):
            # update reactants and scale
            for react in variation.find_all("REACTANT"):
                react_id = react.MOLECULE.attrs["MOL_ID"]
                try:
                    startmats[react_id] = mol_dict[react_id]
                    startmats[react_id]["AMOUNT"] = react.find("AMOUNT").text
                except KeyError:
                    logger.error(
                        f"Reactant ID {react_id} not found in molecules! Skipping this reagent in {rxn_id}"
                    )
        row["scale_mol"] = (
            min(
                [
                    float(v["AMOUNT"].replace("mmol", "").strip())
                    for v in startmats.values()
                ]
            )
            * 0.001
        )
        # Variations continued
        for variation in reaction.find_all("VARIATION"):
            if "CIT_ID" in variation.attrs:
                if "DOI" in citations[variation.attrs["CIT_ID"]]:
                    row["source_id"] = citations[variation.attrs["CIT_ID"]]["DOI"]
                if "TYPE" in citations[variation.attrs["CIT_ID"]]:
                    row["source_type"] = citations[variation.attrs["CIT_ID"]]["TYPE"]
            # update product yields
            for p in variation.find_all("PRODUCT"):
                products[p.find("MOLECULE").attrs["MOL_ID"]].update(
                    {"YIELD": float(p.find("YIELD").text)}
                )
            # parse reagents
            for i, reag in enumerate(variation.find_all("REAGENT")):
                reag_id = reag.MOLECULE.attrs["MOL_ID"]
                try:
                    row[f"reagent_{i+1}_name"] = reag_id
                    row[f"reagent_{i+1}_smiles"] = mol_dict[reag_id]["MOLSTRUCTURE"]
                    if "CAS" in mol_dict[reag_id]:
                        row[f"reagent_{i+1}_cas"] = mol_dict[reag_id]["CAS"]
                    row[f"reagent_{i+1}_eq"] = float(
                        reag.AMOUNT.text.replace("mmol", "").strip()
                    ) / (row["scale_mol"] * 1000)
                except KeyError:
                    logger.error(
                        f"Reagent ID {reag_id} not found in molecules! Skipping this reagent in {rxn_id}"
                    )
            # parse catalysts
            for i, cat in enumerate(variation.find_all("CATALYST")):
                cat_id = cat.MOLECULE.attrs["MOL_ID"]
                try:
                    row[f"catalyst_{i+1}_name"] = cat_id
                    row[f"catalyst_{i+1}_smiles"] = mol_dict[cat_id]["MOLSTRUCTURE"]
                    if "CAS" in mol_dict[cat_id]:
                        row[f"catalyst_{i+1}_cas"] = mol_dict[cat_id]["CAS"]
                    row[f"catalyst_{i+1}_eq"] = float(
                        cat.AMOUNT.text.replace("mmol", "").strip()
                    ) / (row["scale_mol"] * 1000)
                except KeyError:
                    logger.error(
                        f"Catalyst ID {cat_id} not found in molecules! Skipping this catalyst in {rxn_id}"
                    )
            # build solvents
            solvents = {
                r.MOLECULE.attrs["MOL_ID"]: mol_dict[r.MOLECULE.attrs["MOL_ID"]]
                for r in variation.find_all("SOLVENT")
            }
            for r in variation.find_all("SOLVENT"):
                solvents[r.MOLECULE.attrs["MOL_ID"]].update({"AMOUNT": r.AMOUNT.text})
            # parse conditions
            conditions = variation.find_all("CONDITIONS")
            if conditions:
                for cond in conditions:
                    if cond.find("PREPARATION"):
                        row["procedure"] = cond.find("PREPARATION").text
                    if cond.find("TEMPERATURE"):
                        row["temperature_deg_c"] = float(
                            cond.find("TEMPERATURE").text.strip()
                        )
                    if cond.find("TIME"):
                        row["time_h"] = float(cond.find("TIME").text.strip())
                    if cond.find("ATMOSPHERE"):
                        row["atmosphere"] = cond.find("ATMOSPHERE").text

        # write products
        for i, (k, v) in enumerate(products.items()):
            row[f"product_{i+1}_name"] = k
            row[f"product_{i+1}_smiles"] = v["MOLSTRUCTURE"]
            row[f"product_{i+1}_yield"] = v["YIELD"]
            # TODO: check how to get "TYPE" working with UDM
            if "CAS" in v:
                row[f"product_{i+1}_cas"] = v["CAS"]

        # write starting materials
        # get reaction scale and transform mmol to mol
        row["scale_mol"] = (
            max(
                [
                    float(v["AMOUNT"].replace("mmol", "").strip())
                    for v in startmats.values()
                ]
            )
            * 0.001
        )
        for i, (k, v) in enumerate(startmats.items()):
            row[f"startingmat_{i+1}_name"] = k
            row[f"startingmat_{i+1}_smiles"] = v["MOLSTRUCTURE"]
            row[f"startingmat_{i+1}_eq"] = float(
                v["AMOUNT"].replace("mmol", "").strip()
            ) / (row["scale_mol"] * 1000)
            if "CAS" in v:
                row[f"startingmat_{i+1}_cas"] = v["CAS"]

        # write solvents
        v_tot = sum(
            [float(v["AMOUNT"].replace("ml", "").strip()) for v in solvents.values()]
        )
        row["concentration_mol_l"] = v_tot * 1000.0 * row["scale_mol"]
        for i, (k, v) in enumerate(solvents.items()):
            row[f"solvent_{i+1}_name"] = k
            row[f"solvent_{i+1}_smiles"] = v["MOLSTRUCTURE"]
            if "CAS" in v:
                row[f"solvent_{i+1}_cas"] = v["CAS"]
            row[f"solvent_{i+1}_fraction"] = (
                float(v["AMOUNT"].replace("ml", "").strip()) / v_tot
            )

        # add reaction to SURF format
        surf = pd.concat((surf, pd.DataFrame(row).T))

    # save to tabular SURF format
    surf = surf[list(sorted(surf.columns.tolist()))]
    surf = surf.reindex(
        columns=[c for c in const_cols if c in surf.columns]
        + [c for c in surf.columns if (c not in const_cols + last_cols and "product" not in c)]
        + [c for c in surf.columns if "product" in c]
        + [c for c in last_cols if c in surf.columns]
    )
    d_round = {'temperature_deg_c': 1, 'time_h': 0, 'concentration_mol_l': 6, 'scale_mol': 8}
    d_round.update({k: 1 for k in surf.columns if "_yield" in k})
    d_round.update({k: 3 for k in surf.columns if "_eq" in k})
    surf = (
        surf.dropna(how="all", axis=1).convert_dtypes().round(d_round)
    )  # drop empty columns and round to max 8 decimal points
    surf.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    app()
