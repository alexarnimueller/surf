#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description: This file contains utility functions to help other functions in this package.
Author: Alex Müller <alex.mueller@roche.com>
Last updated: 2023-02-21
"""

import logging
import os

import pandas as pd
import rdkit.RDLogger as rkl
import requests
from rdkit.Chem import MolFromSmiles, rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(levelname)-8s %(message)s")
rdklogger = rkl.logger()
rdklogger.setLevel(rkl.CRITICAL)


def std_smiles2mol(smls: str) -> rdchem.Mol:
    """Molecule standardization function to standardize all structures uniformly

    :param smls: SMILES string representing a molecule
    :type smls: str
    :return: standardized molecule object
    :rtype: rdkit.Chem.rdchem.Mol
    """
    un = rdMolStandardize.Uncharger()
    te = rdMolStandardize.TautomerEnumerator()
    mol = MolFromSmiles(smls)
    if mol:
        mol = rdMolStandardize.Cleanup(mol)
        mol = un.uncharge(mol)  # will raise Value Error if invalid SMILES
        return te.Canonicalize(mol)
    return None


def pubchem_name_from_inchikey(df: pd.DataFrame) -> pd.DataFrame:
    """Query PubChem to obtain synonyms for a compound by providing its InChIkey

    :param df: dataframe created in the `surf2db.py` file for compounds
    :type df: pd.DataFrame
    :return: input dataframe but with an updated NAME column
    :rtype: pd.DataFrame
    """
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Calling PubChem..."):
        if pd.isna(row["NAME"]):
            r = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{row['INCHIKEY']}/synonyms/json"
            )
            d = r.json()
            try:
                df.at[idx, "NAME"] = d["InformationList"]["Information"][0]["Synonym"][0]
            except KeyError:
                continue
    return df


def pubchem_property_from_cas(cas: str, property: str) -> str:
    """Query PubChem to obtain InChI strings for a compound by providing its CAS number

    :param inchi: cas number as a string
    :type df: str
    :param property: PubChem compound property to obtain, e.g. smiles, inchi, inchikey etc.
    :return: InChI as a string
    :rtype: str
    """
    try:
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/synonyms/json")
        d = r.json()
        cid = d["InformationList"]["Information"][0]["CID"]
        r = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/json")
        c = r.json()
        rslt = [
            p["value"]["sval"] for p in c["PC_Compounds"][0]["props"] if p["urn"]["label"].upper() == property.upper()
        ][0]
    except (KeyError, IndexError):
        rslt = None
    return rslt
