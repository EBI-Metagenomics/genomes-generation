#!/usr/bin/env python

import argparse
import csv
import sys
import re
import os
import fileinput


def input_args():
    """Multi fasta rename"""
    parser = argparse.ArgumentParser(
        description="Rename multi fasta"
    )
    parser.add_argument(
        "-i", "--input", help="indicate input FASTA file", required=False
    )
    parser.add_argument(
        "-m", "--map", help="map file for names", required=False, default="map.txt"
    )
    parser.add_argument(
        "-p", "--prefix", help="Prefix that would be included to header <prefix><digit>", required=False
    )
    parser.add_argument(
        "-o", "--output", help="indicate output FASTA file", required=True
    )
    args = parser.parse_args()
    return args

def rename(args):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv
    """
    print("Renaming " + args.input)
    with fileinput.hook_compressed(args.input, "r", encoding="utf-8") as fasta_in:
        with open(args.output, "w") as fasta_out, open(args.map, "w") as map_tsv:
            count = 0
            tsv_map = csv.writer(map_tsv, delimiter="\t")
            tsv_map.writerow(["original", "temporary", "short"])
            for line in fasta_in:
                if line.startswith(">"):
                    count += 1
                    fasta_out.write(f">{args.prefix}{count}\n")
                    name = line.strip()
                    short_name = name.split(' ')[0]
                    temporary_name = f"{args.prefix}{count}"
                    tsv_map.writerow([name, temporary_name, short_name])
                else:
                    fasta_out.write(line)
    print(f"Wrote {count} sequences to {args.output}.")


def main():
    args = input_args()
    rename(args)


if __name__ == "__main__":
    main()
