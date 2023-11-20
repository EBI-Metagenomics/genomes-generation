#!/usr/bin/env python3
import argparse
import os
import shutil
import re

def main():
    # set arguments
    # arguments are passed to classes
    parser = argparse.ArgumentParser(
        description="Evaluate completeness and contamination of a MAG."
    )
    parser.add_argument("-o", "--output", dest="outdir", type=str, help="Output directory")
    parser.add_argument("--bins",
        dest="bins",
        nargs="+",
        help="Maxbin bins list"
    )
    parser.add_argument("-v", "--version", dest="version", type=str, help="Tool version")
    parser.add_argument("-a", "--accession", dest="accession", type=str, help="Run accession")

    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    for i in args.bins:
        # Define the regular expression pattern
        filename = os.path.basename(i)
        number = int(filename.split('.')[1])
        new_filename = f"{args.accession}_bin.{number}-maxbin2_v{args.version}.fa"
        shutil.move(i, os.path.join(args.outdir, new_filename))
if __name__ == "__main__":
    main()
