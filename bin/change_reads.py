#!/usr/bin/env python3

"""Changes ERR to ERZ in raw reads and dots to underscores by request"""

from Bio import SeqIO
import argparse
import gzip
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Change names in raw reads"
    )
    parser.add_argument("--reads",
                        dest="reads",
                        nargs="+",
                        help="Raw reads files"
                        )
    parser.add_argument("-f", "--from", dest="from_accession", type=str, help="Raw reads accession")
    parser.add_argument("-t", "--to", dest="to_accession", type=str, help="Assembly accession")
    parser.add_argument("--change_dots_to_underscores", dest="change_dots_to_underscores", action="store_true",
                        help="Specify this flag to change dots to underscores in reads names")
    args = parser.parse_args()

    for read_file in args.reads:
        with gzip.open("changed_" + os.path.basename(read_file), "wt") as out_handle, gzip.open(read_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                new_name = record.id.replace(args.from_accession, args.to_accession)
                if args.change_dots_to_underscores:
                    new_name = new_name.replace('.', '_')
                record.id = new_name
                SeqIO.write(record, out_handle, "fastq")


if __name__ == "__main__":
    main()