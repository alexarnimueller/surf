#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2023, F. Hoffmann La-Roche Ltd.

import logging
import re
from datetime import datetime

import pandas as pd
import typer
from ord_schema import message_helpers, validations
from ord_schema.proto import dataset_pb2, reaction_pb2
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rich.progress import track
from typing_extensions import Annotated

from surf_utils.helpers import pubchem_property_from_cas
from surf_utils.mappings import (
    doi_pattern,
    email_pattern,
    mapping_analyses_ord,
    mapping_atmo,
    mapping_role,
    mapping_stirring,
    orcid_pattern,
)

logging.basicConfig(format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

app = typer.Typer()


@app.command()
def surf2ord(
    input_file: Annotated[str, typer.Argument(help="name of the input file in SURF format")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in ORD format; suffixes: .pbtxt or .pb")],
    delimiter: Annotated[str, typer.Argument(help="delimiter of the input file")] = "\t",
    username: Annotated[str, typer.Option(help="Name of the person submitting the reaction")] = None,
    email: Annotated[str, typer.Option(help="E-mail of the person submitting the reaction")] = None,
    orcid: Annotated[str, typer.Option(help="ORCID of the person submitting the reaction")] = None,
    overwrite_provenance: Annotated[
        bool, typer.Option(help="Whether to overwrite existing provenance with provided info.")
    ] = False,
    organization: Annotated[
        str, typer.Option(help="Organization the person submitting the reaction is associated with")
    ] = None,
    validate_cat_smls: Annotated[
        bool,
        typer.Option(
            help="whether SMILES of should be validated and if invalid, scraped from PubChem using the CAS number"
        ),
    ] = False,
    validate: Annotated[
        bool, typer.Option(help="whether the created ORD reactions should be validated for correctness")
    ] = True,
):
    """Translate reaction data from the tabular SURF format into the protocol buffer format which is used in the open reaction database."""

    if not any([username, orcid, email]):
        logger.warning(
            "No username, orcid or email provided. If reactions don't contain according provenance, they will be rejected! Suggest to add provenance options."
        )
    elif any([username, orcid]) and not email:
        raise ValueError("Both email and username / orcid are required if provenance is provided by user input!")

    if overwrite_provenance:
        if username and email:
            logger.info("Overwriting existing provenance with provided user information")
        else:
            raise ValueError("To overwrite existing provenance, at least a username and email need to be provided!")

    if validate_cat_smls:  # create cache for catalyst lookup
        logger.info("Using validation")
        smls_dict = dict()

    assert output_file.split(".")[-1] in [
        "pbtxt",
        "pb",
    ], "Unknown file extension for output_file! (Known: .pbtxt or .pb)"

    df = pd.read_csv(input_file, delimiter=delimiter)

    # check that headers are all lowercase
    if not all([c.islower() for c in df.columns]):
        logger.warning("Not all headers are lowercase; transforming them to lowercase.")
        df.columns = [c.lower() for c in df.columns]

    # check if yields are 0-1 or 0-100
    if df["product_1_yield"].max() <= 3.0:
        logger.warning("Yields seem to be provided as fractions instead of percentages. Multiplying by 100")
        for c in df.columns:
            if c.endswith("yield"):
                df[c] = df[c] * 100.0

    reactions, invalid = [], []

    for idx, row in track(df.iterrows(), total=len(df), description="Transforming SURF to ORD..."):
        row = row.dropna()
        err = False

        # row = row[np.where(row != 0.0)[0]]  # drop all null values

        reaction = reaction_pb2.Reaction()

        if "rxn_id" not in row:
            logger.error(f"Reaction ID missing! Can't process reaction {idx}")
            continue
        else:
            reaction.identifiers.add(type="NAME", value=row.rxn_id)
        
        # add further identifiers
        rxn_name = [row[n] for n in row.keys() if n in ("rxn_name", "rxn_tech", "rxn_type")]
        if rxn_name:
            reaction.identifiers.add(
                type="CUSTOM",
                value=re.sub(" {2,}", " ", " ".join(rxn_name)).strip(),
                details=re.sub(" {2,}", " ", " ".join(rxn_name)).strip(),
            )

        # Temperature
        if "temperature_deg_c" in row:
            try:
                temp = float(row.temperature_deg_c)
                if temp == 0.0:  # fixes null
                    temp = 0.0001
                cond_t = reaction.conditions.temperature
                cond_t.setpoint.CopyFrom(reaction_pb2.Temperature(units="CELSIUS", value=temp))
            except ValueError:
                try:
                    # find things like 25-100 etc., take highest temp
                    temp = max([float(t) for t in re.findall("[0-9\.]+", row.temperature_deg_c)])
                    cond_t = reaction.conditions.temperature
                    cond_t.setpoint.CopyFrom(reaction_pb2.Temperature(units="CELSIUS", value=temp))
                except ValueError:
                    logger.warning(
                        f"Invalid temperature number value {row.temperature_deg_c} provided for reaction {row.rxn_id}; ignoring temperature!"
                    )

        # Atmosphere
        if "atmosphere" in row:
            cond_p = reaction.conditions.pressure
            cond_p.control.type = reaction.conditions.pressure.PressureControl.SEALED
            cond_p.atmosphere.type = mapping_atmo.get(
                row.atmosphere.upper().replace(" ", "_"), reaction.conditions.pressure.Atmosphere.UNSPECIFIED
            )
        # Stirring
        if "stirring_shaking" in row:
            cond_s = reaction.conditions.stirring
            if row["stirring_shaking"].upper() in mapping_stirring.keys():
                cond_s.type = mapping_stirring[row["stirring_shaking"].upper()]

        # procedure text
        if "procedure" in row and len(row.procedure):
            reaction.notes.procedure_details = row["procedure"]

        # provenance either from SURF or overwritten from user input
        reaction.provenance.record_created.time.value = datetime.today().strftime("%Y-%m-%d")
        if "source_id" in row and re.match(doi_pattern, row.source_id):  # add DOI if present
            reaction.provenance.doi = re.findall(doi_pattern, row.source_id)[0]

        # add provided provenance
        if username:
            reaction.provenance.record_created.person.username = username
        if email:
            reaction.provenance.record_created.person.email = email
        if orcid and re.match(orcid_pattern, orcid):
            reaction.provenance.record_created.person.orcid = orcid
        if organization:
            reaction.provenance.record_created.person.organization = organization

        # if provenance is present, overwrite
        if "provenance" in row and not overwrite_provenance:
            if re.findall(email_pattern, row.provenance):  # add Email if present
                mail = re.findall(email_pattern, row.provenance)[0]
                reaction.provenance.record_created.person.email = mail
                row["provenance"] = row.provenance.replace(mail, "").strip()
            if re.findall(orcid_pattern, row.provenance) and mail:  # add ORCID if present
                orc = re.findall(orcid_pattern, row.provenance)[0]
                reaction.provenance.record_created.person.orcid = orc
                row["provenance"] = row.provenance.replace(orc, "").strip()
            if re.findall("(?:[A-Za-z,.'-]+ ){1,3}", row.provenance) and mail:  # find something like a name
                user = re.findall("(?:[A-Za-z,.'-]+ ){1,3}", row.provenance)[0].strip()
                reaction.provenance.record_created.person.username = user
                row["provenance"] = row.provenance.replace(user, "").strip()
            if len(row.provenance) and mail and user:  # if still something left, add as organization
                reaction.provenance.record_created.person.organization = row.provenance

        # outcome template
        outcome = reaction.outcomes.add()

        # time
        if "time_h" in row:
            try:
                time_h = float(row.time_h)
                outcome.reaction_time.CopyFrom(reaction_pb2.Time(value=time_h, units="HOUR"))
            except ValueError:
                try:
                    # find things like 4-96 etc., take longest time
                    time_h = max([float(t) for t in re.findall("[0-9\.]+", row.time_h)])
                    outcome.reaction_time.CopyFrom(reaction_pb2.Time(value=time_h, units="HOUR"))
                    logger.warning(
                        f"Guessing time value '{time_h}' from '{row.time_h}' provided for reaction {row.rxn_id}"
                    )
                except ValueError:
                    logger.warning(
                        f"Invalid time number value '{row.time_h}' provided for reaction {row.rxn_id}; ignoring time!"
                    )

        # main loop through SURF headers
        for cpd in sorted(set(re.findall("\w+_[1-9]", " ".join(row.index.tolist())))):
            if f"{cpd}_smiles" not in row and f"{cpd}_cas" not in row and f"{cpd}_name" not in row:
                continue

            # product side
            if "product" in cpd:
                # Identifier
                product = outcome.products.add(identifiers=[{"type": "SMILES", "value": row[f"{cpd}_smiles"]}])
                product.reaction_role = reaction_pb2.ReactionRole.PRODUCT
                if f"{cpd}_cas" in row or f"{cpd}_name" in row:
                    product.identifiers.add(type="NAME", value=row.get(f"{cpd}_cas", row.get(f"{cpd}_name")))

                # MS analytics
                if f"{cpd}_ms" in row:
                    product.measurements.add(analysis_key=f"MS of {cpd}", type="IDENTITY")
                    ms_type = "MS"
                    if isinstance(row[f"{cpd}_ms"], str) and "HRMS" in row[f"{cpd}_ms"].upper():
                        ms_type = "HRMS"
                    outcome.analyses[f"MS of {cpd}"].CopyFrom(
                        reaction_pb2.Analysis(type=ms_type, details=f'{row[f"{cpd}_ms"]}')
                    )
                    outcome.analyses[f"MS of {cpd}"].data["found"].description = "Observed m/z"
                    try:
                        mass = float(row[f"{cpd}_ms"])
                        outcome.analyses[f"MS of {cpd}"].data["found"].float_value = mass
                        outcome.analyses[f"MS of {cpd}"].data["found"].description = "Observed m/z"
                    except ValueError:
                        try:
                            # try to match things like "Found: XX.XX"
                            mass = float(
                                "".join(
                                    re.findall(
                                        "[0-9\.]",
                                        re.findall("[Ff]ound[:\s\w\(\)\+\-]+[0-9\.]+", row[f"{cpd}_ms"])[-1].strip(
                                            "."
                                        ),
                                    )
                                )
                            )
                            outcome.analyses[f"MS of {cpd}"].data["found"].float_value = mass
                            outcome.analyses[f"MS of {cpd}"].data["found"].description = "Observed m/z"
                            logger.warning(
                                f"Guessing MS value '{mass}' from '{row[f'{cpd}_ms']}' provided for reaction {row.rxn_id}"
                            )
                            try:  # try to match things like "Expected: XX.XX"
                                xpctd = float(
                                    "".join(
                                        re.findall(
                                            "[0-9\.]",
                                            re.findall("[Ee]xpect[ed]*[:\s\w\(\)\+\-]+[0-9\.]+", row[f"{cpd}_ms"])[
                                                -1
                                            ].strip("."),
                                        )
                                    )
                                )
                                outcome.analyses[f"MS of {cpd}"].data["expected"].float_value = xpctd
                                outcome.analyses[f"MS of {cpd}"].data["expected"].description = "Expected m/z"
                            except (IndexError, ValueError):
                                pass
                        except (IndexError, ValueError):
                            logger.warning(
                                f"Invalid MS number value {row[f'{cpd}_ms']} provided for reaction {row.rxn_id}; ignoring value!"
                            )

                # NMR analytics
                if f"{cpd}_nmr" in row:
                    product.measurements.add(analysis_key=f"NMR of {cpd}", type="IDENTITY")
                    outcome.analyses[f"NMR of {cpd}"].CopyFrom(reaction_pb2.Analysis(type="NMR_1H"))
                    outcome.analyses[f"NMR of {cpd}"].data["peaks"].string_value = row[f"{cpd}_nmr"]

                # Yield and according measurements
                if f"{cpd}_yield" in row:
                    yield_type = "CUSTOM"
                    if f"{cpd}_yieldtype" in row:
                        if row[f"{cpd}_yieldtype"].upper() in mapping_analyses_ord.keys():
                            yield_type = row[f"{cpd}_yieldtype"].upper()
                        elif "NMR" in row[f"{cpd}_yieldtype"].upper():
                            yield_type = "NMR_1H"
                        elif "MS" in row[f"{cpd}_yieldtype"].upper():
                            yield_type = "MS"
                        elif "GC" in row[f"{cpd}_yieldtype"].upper():
                            yield_type = "GC"
                        outcome.analyses[f"{cpd}_{row[f'{cpd}_yieldtype']}"].CopyFrom(
                            reaction_pb2.Analysis(type=yield_type, details=row[f"{cpd}_yieldtype"])
                        )
                        # outcome.analyses[f"{cpd}_{row[f'{cpd}_yieldtype']}"].is_of_isolated_species = True
                        product.measurements.add(
                            type="YIELD",
                            analysis_key=f"{cpd}_{row[f'{cpd}_yieldtype']}",
                            percentage={"value": 100 * float(row[f"{cpd}_yield"])},
                        )
                    else:
                        # logger.warning(f"Yield type not defined for {cpd}; using 'unknown'")
                        outcome.analyses[f"{cpd}_unkown"].CopyFrom(
                            reaction_pb2.Analysis(type=yield_type, details="unkown")
                        )
                        # outcome.analyses[f"{cpd}_unkown"].is_of_isolated_species = True
                        product.measurements.add(
                            type="YIELD", analysis_key=f"{cpd}_unkown", percentage={"value": 100 * row[f"{cpd}_yield"]}
                        )

            # Educt side
            else:
                # Solvent and volume (ml)
                if "solvent" in cpd and "scale_mol" in row and "concentration_mol_l" in row:
                    try:
                        scale = float(row.scale_mol)
                        conc = float(row.concentration_mol_l)
                        if f"{cpd}_fraction" in row:  # check if fraction available (multiple solvents present)
                            vol = f"{row[f'{cpd}_fraction'] * scale / conc * 1000} ml"
                        else:
                            vol = f"{scale / conc * 1000} ml"
                        name = row.get(f"{cpd}_cas", row.get(f"{cpd}_name", "solvent"))
                        reaction.inputs[re.sub("_[0-9]+", "", cpd)].components.add().CopyFrom(
                            message_helpers.build_compound(
                                smiles=row[f"{cpd}_smiles"],
                                name=name,
                                role="solvent",
                                amount=vol,
                            )
                        )
                    except (ValueError, ZeroDivisionError):
                        logger.warning(
                            f"Invalid scale '{row.scale_mol}' or conc. '{row.concentration_mol_l}' provided in reaction {row.rxn_id}; unable to process reaction entry!"
                        )
                        err = True

                # Educt with amount (mmol)
                elif f"{cpd}_eq" in row and "scale_mol" in row and row[f"{cpd}_eq"] > 0 and row.scale_mol > 0:
                    try:
                        smls, inchi = row.get(f"{cpd}_smiles", None), None
                        if "catalyst" in cpd and validate_cat_smls and smls:
                            if not MolFromSmiles(smls) and f"{cpd}_cas" in row:
                                if smls not in smls_dict:  # try messing with charges
                                    tmp1 = MolFromSmiles(smls.replace("--", "-2").replace("++", "+2"))
                                    tmp2 = MolFromSmiles(smls.replace("+", "").replace("-", ""))
                                    if tmp1 or tmp2:
                                        tmp = MolToSmiles(tmp1) if tmp1 else MolToSmiles(tmp2)
                                        smls_dict[row[f"{cpd}_smiles"]] = (tmp, "smiles")
                                        smls = tmp
                                    elif re.match("^[0-9]+-[0-9]+-[0-9]+$", row[f"{cpd}_cas"]):
                                        smls = pubchem_property_from_cas(row[f"{cpd}_cas"], "SMILES")
                                        if not smls or not MolFromSmiles(smls):  # validate if obtained SMILES is valid
                                            inchi = pubchem_property_from_cas(row[f"{cpd}_cas"], "InChI")
                                            if inchi:
                                                smls_dict[row[f"{cpd}_smiles"]] = (
                                                    inchi,
                                                    "inchi",
                                                )  # store to not call again
                                        else:
                                            smls_dict[row[f"{cpd}_smiles"]] = (smls, "smiles")
                                else:  # compound was already queried before, get from cache
                                    rslt, prop = smls_dict[smls]
                                    if prop == "smiles":
                                        smls = rslt
                                    else:
                                        inchi = rslt

                        eq = float(row[f"{cpd}_eq"])
                        scale = float(row.scale_mol)
                        name = row.get(f"{cpd}_cas", row.get(f"{cpd}_name", "reagent"))
                        if inchi and name:  # when PubChem returned InChI
                            reaction.inputs[re.sub("_[0-9]+", "", cpd)].components.add().CopyFrom(
                                message_helpers.build_compound(
                                    inchi=inchi,
                                    name=name,
                                    role=mapping_role.get(cpd[:-2], "UNSPECIFIED"),
                                    amount=f"{scale * eq * 1000} mmol",
                                )
                            )
                        elif smls or name:  # use SMILES
                            reaction.inputs[re.sub("_[0-9]+", "", cpd)].components.add().CopyFrom(
                                message_helpers.build_compound(
                                    smiles=smls if smls else "",
                                    name=name,
                                    role=mapping_role.get(cpd[:-2], "UNSPECIFIED"),
                                    amount=f"{scale * eq * 1000} mmol",
                                )
                            )
                    except ValueError:
                        logger.warning(
                            f"Invalid {cpd}_eq number value {row[f'{cpd}_eq']} or scale number value {row.scale_mol} in reaction {row.rxn_id}; unable to process reaction entry!"
                        )
                        err = True

        # store reaction if all good
        if not err:
            if validate:
                try:
                    validations.validate_message(reaction)
                except Warning as w:
                    logger.warning(f"Not adding reaction invalid reaction {row.rxn_id} du to: {w} !")
                    invalid.append(reaction.identifiers[0].value)
                    continue
            reactions.append(reaction)

    # Validation and writing
    dataset = dataset_pb2.Dataset(reactions=reactions)

    if validate:
        logger.info("Running final ORD dataset validation...")
        validations.validate_datasets({input_file: dataset})
    #     logger.warning("The following reactions were invalid:")
    #     logger.warning(", ".join(invalid))

    logger.info(f"Writing ORD file {output_file}")
    message_helpers.write_message(dataset, output_file)


if __name__ == "__main__":
    app()
