#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2023, F. Hoffmann La-Roche Ltd.

import logging
import re

import pandas as pd
import typer
from typing_extensions import Annotated

from surf_utils.helpers import const_cols, last_cols, rxn_smiles_left, rxn_smiles_right

app = typer.Typer()


@app.command()
def surf2rxn(
    input_file: Annotated[str, typer.Argument(help="name of the input file in SURF format")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in reaction SMILES format")],
    delimiter: Annotated[str, typer.Argument(help="delimiter of the input file")] = "\t",
):
    """Translate a SURF file into reaction SMILES format"""
    # read SURF reaction data file
    df = pd.read_csv(input_file, delimiter=delimiter)
    df = df.replace("n/a", "").replace("N/A", "").replace("#N/A", "").replace("nan", "")  # remove NA placeholders
    df.dropna(how="all", axis=1, inplace=True)  # drop completely empty columns
    logging.info(f"{len(df)} entries read from SURF file {input_file}")

    # check that headers are all lowercase
    if not all([c.islower() for c in df.columns]):
        logging.warning("Not all headers are lowercase; transforming them to lowercase.")
        df.columns = [c.lower() for c in df.columns]

    # check if yields are 0-1 or 0-100
    if df["product_1_yield"].max() <= 3.0:
        logging.warning("Yields seem to be provided as fractions instead of percentages. Multiplying by 100")
        for c in df.columns:
            if c.endswith("yield"):
                df[c] = df[c] * 100.0

    dots = re.compile("[\.]{2,}")

    # find reaction components for the left side of the arrow
    r = re.compile(rxn_smiles_left)
    tmp = df[list(filter(r.match, df.columns.tolist()))].astype(str)
    tmp = tmp.replace("nan", "")
    left = tmp.agg(".".join, axis=1).str.replace(dots, ".", regex=True).str.strip(".")

    # find solvents
    r = re.compile("solvent_[1-9]_smiles")
    tmp = df[list(filter(r.match, df.columns.tolist()))].astype(str)
    tmp = tmp.replace("nan", "")
    solv = tmp.agg(".".join, axis=1).str.replace(dots, ".", regex=True).str.strip(".")

    # find reaction components for the right side of the arrow (products)
    r = re.compile(rxn_smiles_right)
    tmp = df[list(filter(r.match, df.columns.tolist()))].astype(str)
    tmp = tmp.replace("nan", "")
    right = tmp.agg(".".join, axis=1).str.replace(dots, ".", regex=True).str.strip(".")

    # keep only relevant columns and add left and right components of rxn smiles
    cols = [c for c in const_cols if c in df.columns] + [c for c in last_cols if c in df.columns]
    df = df[cols]
    df["rxn_smiles"] = left + ">" + solv + ">" + right
    df["rxn_smiles"] = df["rxn_smiles"].str.replace(dots, ".", regex=True)

    # drop empty columns and round to max 8 decimal points
    df = df.dropna(how="all", axis=1).convert_dtypes().round(8)

    # save to rxn smiles file
    logging.info(f"writing {len(df)} entries to RXN SMILES file {output_file}...")
    df.to_csv(output_file, sep=delimiter, index=False)


if __name__ == "__main__":
    app()
