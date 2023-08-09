#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2023, F. Hoffmann La-Roche Ltd.

import logging
import re
import xml.etree.ElementTree as ET
from datetime import datetime

import pandas as pd
import typer
from rdkit.Chem import MolFromInchi, MolFromSmiles, MolToInchiKey, MolToMolBlock
from rich.progress import track
from typing_extensions import Annotated

from surf_utils.mappings import doi_pattern

logging.basicConfig(format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


app = typer.Typer()


def remove_none_attrib(elem):
    """Remove None attributes from XML tree"""
    elem.attrib = {k: elem.attrib[k] for k in elem.attrib if elem.attrib[k] is not None}
    for subelem in elem:
        remove_none_attrib(subelem)


@app.command()
def surf2ord(
    input_file: Annotated[str, typer.Argument(help="name of the input file in SURF format")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in XML format; suffix: .xml")],
    delimiter: Annotated[str, typer.Argument(help="delimiter of the input file")] = "\t",
):
    # read the SURF file
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

    # define UDM headers and main structure
    udm = ET.Element("UDM", DATABASE="SURF", SEQUENCE="001", TIMESTAMP=datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))
    udm.set("xmlns:xsi", '"http://www.w3.org/2001/XMLSchema-instance"')
    udm.set("xmlns:xsd", "http://www.w3.org/2001/XMLSchema")
    version = ET.SubElement(udm, "UDM_VERSION", MAJOR="6", MINOR="0", REVISION="0", VERSIONTEXT="6.0.0")
    legal = ET.SubElement(udm, "LEGAL")
    citations = ET.SubElement(udm, "CITATIONS")
    molecules = ET.SubElement(udm, "MOLECULES")
    reactions = ET.SubElement(udm, "REACTIONS")

    mol_dict = dict()
    cit_dict = dict()
    cnt_cit = 1

    invalid = []
    for idx, row in track(df.iterrows(), total=len(df), description="Reading molecules..."):
        row = row.dropna()
        err = False

        if "rxn_id" not in row:
            logger.error(f"Reaction ID missing! Can't process reaction {idx}")
            continue

        # reading compounds of this reaction and adding them to the molecule dictionary
        for cpd in sorted(set(re.findall("\w+_[1-9]", " ".join(row.index.tolist())))):
            mol, cas = None, row.get(f"{cpd}_cas", None)
            if f"{cpd}_smiles" in row:  # process SMILES if present
                mol = MolFromSmiles(row[f"{cpd}_smiles"])
            elif f"{cpd}_inchi" in row:  # process InChI if present
                mol = MolFromInchi(row[f"{cpd}_inchi"])
            if mol:
                mol_dict[MolToInchiKey(mol)] = {
                    "MOLSTRUCTURE": MolToMolBlock(mol)
                }  # store the compound using the inchikey as index
                if cas:
                    cas = re.sub(r"-0([0-9])$", r"-\g<1>", cas)  # replace potential leading zero for last CAS pos
                    mol_dict[MolToInchiKey(mol)].update({"CAS": cas})
            else:
                logger.error(f"Can't process SMILES/InChI of {cpd} in reaction {row.rxn_id}")
                err = True
                continue

        # reading the provenance / citations
        if "source_id" in row and re.match(doi_pattern, row.source_id):  # catch DOI if present
            doi = re.findall(doi_pattern, row.source_id)[0]
            if doi not in cit_dict:
                cit_dict[doi] = cnt_cit
                cit = ET.SubElement(citations, "CITATION", ID=str(cnt_cit))
                ET.SubElement(cit, "DOI").text = doi
                cnt_cit += 1

        # exclude reaction if error observed in one of the molecules of this reaction
        if err:
            invalid.append(row.rxn_id)
            logger.error(f"Couln't process all SMILES in reaction {row.rxn_id}, skipping reaction")

    # processing the individual, unique compounds
    for k, v in track(mol_dict.items(), total=len(mol_dict), description="Adding molecules to UDM..."):
        mol = ET.SubElement(molecules, "MOLECULE", ID=k)
        ET.SubElement(mol, "MOLSTRUCTURE").text = "<![CDATA[" + v["MOLSTRUCTURE"] + "]]>"
        if "CAS" in v:
            ET.SubElement(mol, "CAS").text = v["CAS"]

    if invalid:
        logger.warning(f"Ignoring {len(invalid)} reactions where SMILES could not be read for one of the components")
        df = df[~df.rxn_id.isin(invalid)]  # exclude invalid reactions

    # processing the individual reactions
    for idx, row in track(df.iterrows(), total=len(df), description="Processing reactions..."):
        if "rxn_id" not in row:
            logger.error(f"Reaction ID missing! Can't process reaction {idx}")
            continue
        row = row.dropna()
        err = False

        reaction = ET.SubElement(reactions, "REACTION", ID=row.rxn_id.replace("_", "-"))

        # add a variation, sections and conditions
        if "source_id" in row and re.match(doi_pattern, row.source_id):  # catch DOI if present
            doi = re.findall(doi_pattern, row.source_id)[0]
            variation = ET.SubElement(reaction, "VARIATION", CIT_ID=str(cit_dict[doi]))
        else:
            variation = ET.SubElement(reaction, "VARIATION")
        section = ET.SubElement(variation, "SECTION")
        conditions = ET.SubElement(section, "CONDITIONS")

        if "procedure" in row:
            cond_glob = ET.SubElement(section, "CONDITIONS")
            ET.SubElement(cond_glob, "PREPARATION").text = row.procedure
        if "temperature_deg_c" in row:
            temp = ET.SubElement(conditions, "TEMPERATURE")  # default unit is degC
            ET.SubElement(temp, "exact").text = f"{row.temperature_deg_c:.2f}"
        if "time_h" in row:
            tme = ET.SubElement(conditions, "TIME")  # default unit is hour
            ET.SubElement(tme, "exact").text = f"{row.time_h:.1f}"

        for cpd in sorted(set(re.findall("\w+_[1-9]", " ".join(row.index.tolist())))):
            if (
                f"{cpd}_smiles" not in row
                and f"{cpd}_inchi" not in row
                and f"{cpd}_cas" not in row
                and f"{cpd}_name" not in row
            ):
                logger.error(f"no SMILES, InChI, CAS or Name forund for {cpd} in {row.rxn_id}, skipping reaction")
                err = True
                continue

            # get the compound name
            cpd_name = row.get(
                f"{cpd}_cas", row.get(f"{cpd}_name", row.get(f"{cpd}_smiles", row.get(f"{cpd}_inchi", "unknown")))
            )

            # create an InChI Key as the identifier to link molecules to reactions
            if f"{cpd}_smiles" in row:
                inchikey = MolToInchiKey(MolFromSmiles(row[f"{cpd}_smiles"]))
            elif f"{cpd}_inchi" in row:
                inchikey = MolToInchiKey(MolFromInchi(row[f"{cpd}_inchi"]))
            else:
                logger.error(f"no SMILES/InChI forund for {cpd} in {row.rxn_id}, skipping reaction")
                err = True
                continue

            # handle the compounds in the reaction
            if "product" in cpd:
                ET.SubElement(section, "PRODUCT_ID").text = inchikey
                prod = ET.SubElement(section, "PRODUCT")
                ET.SubElement(prod, "MOLECULE", MOL_ID=inchikey).text = ""
                ET.SubElement(prod, "NAME").text = cpd_name
                if f"{cpd}_yield" in row:
                    ET.SubElement(prod, "YIELD").text = f"{row[f'{cpd}_yield']:.1f}"
            elif "startingmat" in cpd:
                ET.SubElement(section, "REACTANT_ID").text = inchikey
                reactant = ET.SubElement(section, "REACTANT")
                ET.SubElement(reactant, "MOLECULE", MOL_ID=inchikey).text = ""
                ET.SubElement(reactant, "NAME").text = cpd_name
                # check if equivalents and a reaction scale are provided
                if f"{cpd}_eq" in row and "scale_mol" in row and row.scale_mol > 0:
                    rct_amnt = ET.SubElement(reactant, "AMOUNT")
                    rct_amnt.text = f"{row.scale_mol * row[f'{cpd}_eq'] * 1000:.4f} mmol"  # mmol
            else:
                if "solvent" in cpd:
                    solvent = ET.SubElement(section, "SOLVENT")
                    ET.SubElement(solvent, "MOLECULE", MOL_ID=inchikey).text = ""
                    ET.SubElement(solvent, "NAME").text = cpd_name
                    if (
                        "scale_mol" in row
                        and "concentration_mol_l" in row
                        and row["concentration_mol_l"] > 0
                        and row["scale_mol"] > 0
                    ):
                        try:
                            scale = float(row.scale_mol)
                            conc = float(row.concentration_mol_l)
                            solv_amnt = ET.SubElement(solvent, "AMOUNT")
                            if f"{cpd}_fraction" in row:  # check if fraction available (multiple solvents present)
                                solv_amnt.text = f'{row[f"{cpd}_fraction"] * scale / conc * 1000:.4f} ml'  # ml
                            else:
                                solv_amnt.text = f"{scale / conc * 1000:.4f} ml"  # ml
                        except ValueError:
                            logger.warning(f"Can't process {cpd}_fraction to amount in reaction {row.rxn_id}")
                elif "catalyst" in cpd:
                    cat = ET.SubElement(section, "CATALYST")
                    ET.SubElement(cat, "MOLECULE", MOL_ID=inchikey).text = ""
                    ET.SubElement(cat, "NAME").text = cpd_name
                    if f"{cpd}_eq" in row and "scale_mol" in row and row.scale_mol > 0:
                        cat_amnt = ET.SubElement(cat, "AMOUNT")
                        cat_amnt.text = f"{row.scale_mol * row[f'{cpd}_eq'] * 1000:.4f} mmol"  # mmol
                else:
                    reagent = ET.SubElement(section, "REAGENT")
                    ET.SubElement(reagent, "MOLECULE", MOL_ID=inchikey).text = ""
                    ET.SubElement(reagent, "NAME").text = cpd_name
                    if f"{cpd}_eq" in row and "scale_mol" in row and row.scale_mol > 0:
                        rgnt_amnt = ET.SubElement(reagent, "AMOUNT")
                        rgnt_amnt.text = f"{row.scale_mol * row[f'{cpd}_eq'] * 1000:.4f} mmol"  # mmol

    remove_none_attrib(udm)
    ET.indent(udm)
    with open(output_file, "w") as f:
        out = ET.tostring(udm, method="xml", encoding="utf-8", xml_declaration=True, short_empty_elements=True)
        # write the UDM to XML, replacing the escaped CDATA characters again
        f.write(
            out.decode("utf-8")
            .replace("<MOLSTRUCTURE>&lt;", "<MOLSTRUCTURE><")
            .replace("&gt;</MOLSTRUCTURE>", "></MOLSTRUCTURE>")
            .replace("&quot;", "")
        )


if __name__ == "__main__":
    app()
