#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (©) 2024, F. Hoffmann La-Roche Ltd.

import os
import subprocess
import threading
from queue import Queue

import typer
from typing_extensions import Annotated

app = typer.Typer()


def process_one_file(input_file: str, output_file: str, validate: str):
    subprocess.run(["python", "ord2surf.py", input_file, output_file, validate, "--verbose", "0"])
    Q.get()


@app.command()
def whole_ord_to_surf(
    ord_data_path: Annotated[str, typer.Argument(help="path to the ORD data folder")],
    output_path: Annotated[str, typer.Argument(help="path to the output folder")],
    validate: Annotated[
        bool, typer.Option(help="whether the input ORD reactions should be validated for correctness first")
    ] = False,
    n_jobs: Annotated[int, typer.Option(help="number of parallel jobs to run")] = 4,
    ignore_existing: Annotated[bool, typer.Option(help="ignore inputs with existing result files")] = True,
):
    """Transform all ORD data into SURF format.
    The output is saved in the same folders as the input pb.gz files as SURF .txt files
    """
    msg = "Not v" if not validate else "V"
    print(msg + "alidating input ORD reactions")
    validate = "--validate" if validate else "--no-validate"
    print(f"Using {n_jobs} parallel threads for processing")

    # create output folder if it does not exist
    os.makedirs(output_path, exist_ok=True)

    # read input files
    all_input_files = []
    for subdir, _, files in os.walk(ord_data_path):
        for file in files:
            if file.endswith(".pb.gz"):
                all_input_files.append(os.path.join(subdir, file))
    all_input_files = list(sorted(all_input_files))  # sorting alphabetically

    if ignore_existing:  # ignore input files with existing output files
        print("Ignoring input files with existing output files")
        all_input_files = [
            f
            for f in all_input_files
            if not os.path.exists(os.path.join(output_path, os.path.basename(f.replace(".pb.gz", ".txt"))))
        ]
    print(f"Found {len(all_input_files)} ORD input files to process")

    # multithreading, process n_jobs files at a time
    fidx = 0
    threads = []
    global Q
    Q = Queue(n_jobs)
    for i, input_file in enumerate(all_input_files):
        out_file = os.path.join(output_path, os.path.basename(input_file.replace(".pb.gz", ".txt")))
        t = threading.Thread(
            target=process_one_file,
            args=(input_file, out_file, validate),
        )
        Q.put(i)
        threads.append(t)
        t.start()
        fidx += 1

    [t.join() for t in threads]  # wait for all threads to terminate

    # check how many output files are present
    all_output_files = []
    for subdir, _, files in os.walk(output_path):
        for file in files:
            if file.endswith(".txt"):
                all_output_files.append(os.path.join(subdir, file))

    # check if all input files were processed, report the ones that were not
    print(f"Generated {len(all_output_files)} SURF output files in {output_path}")
    if len(all_output_files) < len(all_input_files):
        print("WARNING: Some input files failed to be processed:")
        for f in all_input_files:
            out_file = os.path.basename(f.replace(".pb.gz", ".txt"))
            if not os.path.exists(os.path.join(output_path, out_file)):
                print(f"\t{f}")


if __name__ == "__main__":
    app()
