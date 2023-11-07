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
from typing_extensions import Annotated

app = typer.Typer()

# TODO: get complete UDM example and finalize

logging.basicConfig(format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
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
    input_file: Annotated[str, typer.Argument(help="name of the input file in UDM format")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in ORD format; suffixes: .pbtxt or .pb")],
    delimiter: Annotated[str, typer.Argument(help="delimiter of the input file")] = "\t",
    validate: Annotated[
        bool, typer.Option(help="Whether to overwrite existing provenance with provided info.")
    ] = True,
):
    if validate:  # validate the XML file versus the UDM schema
        errors = validate_udm_file(input_file)
        if errors:
            logger.error(f"UDM schema validation found the following {len(errors)} errors:")
            for e in errors:
                logger.error(e)
        else:
            logger.info("Validation successful")

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
    citations = xml.find("CITATIONS")

    # read all molecules
    molecules = xml.find("MOLECULES")
    mol_dict = {}
    for mol in molecules.find_all("MOLECULE"):
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
    for reaction in reactions.find_all("REACTION"):
        row = pd.Series()
        rxn_id = reaction.attrs["ID"]
        row["rxn_id"] = rxn_id
        # Reactant
        try:
            react_ids = [r.text for r in reaction.find_all("REACTANT_ID")]
            startmats = {i: mol_dict[i] for i in react_ids}
        except KeyError:
            logger.error(
                f"One of these reactant IDs {react_ids} was not found in molecules! Skipping reaction {rxn_id}"
            )
            continue
        # Product
        try:
            prod_ids = [p.text for p in reaction.find_all("PRODUCT_ID")]
            products = {i: mol_dict[i] for i in prod_ids}
        except KeyError:
            logger.error(f"One of these product IDs {prod_ids} was not found in molecules! Skipping reaction {rxn_id}")
            continue
        # Variations
        for variation in reaction.find_all("VARIATION"):
            yields = dict()
            for p in variation.find_all("PRODUCT"):
                yields[p.find("MOLECULE").attrs["MOL_ID"]] = float(p.find("YIELD").text)
            reagents = dict()
            for reagent in variation.find_all("REAGENT"):
                reag_id = reagent.MOLECULE.attrs["MOL_ID"]
                try:
                    reagents[reag_id] = mol_dict[reag_id]
                    reagents[reag_id]["AMOUNT"] = reagent.find("AMOUNT").text
                except KeyError:
                    logger.error(f"Reagent ID {reag_id} not found in molecules! Skipping this reagent in {rxn_id}")
            catalysts = dict()
            for cat in variation.find_all("CATALYST"):
                cat_id = cat.MOLECULE.attrs["MOL_ID"]
                try:
                    catalysts[cat_id] = mol_dict[cat_id]
                    catalysts[cat_id]["AMOUNT"] = cat.find("AMOUNT").text
                except KeyError:
                    logger.error(f"Catalyst ID {cat_id} not found in molecules! Skipping this catalyst in {rxn_id}")
            reactants = dict()
            for react in variation.find_all("REACTANT"):
                react_id = react.MOLECULE.attrs["MOL_ID"]
                try:
                    reactants[react_id] = mol_dict[react_id]
                    reactants[react_id]["AMOUNT"] = react.find("AMOUNT").text
                except KeyError:
                    logger.error(f"Reactant ID {react_id} not found in molecules! Skipping this reagent in {rxn_id}")
            conditions = variation.find("CONDITIONS")
            if conditions:
                if conditions.find("PREPARATION"):
                    row["procedure"] = conditions.find("PREPARATION").text
                if conditions.find("TEMPERATURE"):
                    row["temperature_deg_c"] = float(conditions.find("TEMPERATURE").text.strip())
                if conditions.find("TIME"):
                    row["time_h"] = float(conditions.find("TIME").text.strip())

    # list(xsd.find("xs:element", {"name": "VARIATION"}).find("xs:complexType").find("xs:sequence").children)


if __name__ == "__main__":
    # raise NotImplementedError("Not yet implemented, coming soon!")
    app()
