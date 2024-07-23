#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os

import pandas as pd
import typer
from rich.progress import track
from typing_extensions import Annotated

from surf_utils.mappings import const_cols, last_cols

app = typer.Typer()


@app.command()
def concat_surfs(
    input_folder: Annotated[str, typer.Argument(help="name of the input folder containing SURF files")],
    output_file: Annotated[str, typer.Argument(help="name of the output file in SURF format; suffixes: .txt or .tsv")],
    delimiter: Annotated[str, typer.Argument(help="delimiter of the input file")] = "\t",
):
    """Concatenate multiple SURF files in the same directory into a single output SURF file."""

    surf = pd.DataFrame()

    for f in track(os.listdir(input_folder), description=f"Reading files in {input_folder}..."):
        if f.endswith(".tsv") or f.endswith(".txt"):
            tmp = pd.read_csv(os.path.join(input_folder, f), delimiter=delimiter)
            surf = pd.concat((surf, tmp))

    # save to tabular SURF format
    surf = surf[list(sorted(surf.columns.tolist()))]
    surf = surf.reindex(
        columns=[c for c in const_cols if c in surf.columns]
        + [c for c in surf.columns if c not in const_cols + last_cols]
        + [c for c in last_cols if c in surf.columns]
    )
    surf = (
        surf.dropna(how="all", axis=1).convert_dtypes().round(8)
    )  # drop empty columns and round to max 8 decimal points
    surf.to_csv(output_file, sep="\t", index=False, float_format="%f")


if __name__ == "__main__":
    app()
