#!/usr/bin/env python3
#
# This file is part of the EukCC (https://github.com/openpaul/eukcc).
# Copyright (c) 2019 Paul Saary
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# provides all file operation functions
# used inthis package
#
# Modified by Alejandra Escobar (30 Jun 2023) for the genomes-generation 
# pipeline of MGnify

import os
import argparse
import subprocess
import logging
import tempfile
import gzip
import glob
import csv


# backup fasta handler, so we can use readonly directories
class fa_class:
    def __init__(self, seq, name, long_name):
        self.seq = seq
        self.name = name
        self.long_name = long_name

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)


def busco_score(path):
    d = {
        "BUSCO_C": None,
        "BUSCO_M": None,
        "BUSCO_D": None,
        "BUSCO_S": None,
        "BUSCO_F": None,
    }

    with open(path) as fin:
        for line in fin:
            if line.strip().startswith("C:"):
                elem = line.strip().replace("[", ",").replace("]", "").split(",")
                for e in elem:
                    x = e.split(":")
                    d["BUSCO_{}".format(x[0])] = x[1].replace("%", "")
                # this is the last info we need
                break
            if line.strip().startswith("# The lineage dataset is: "):
                d["BUSCO_lineage"] = line.strip().split()[5]
    return d

def genome_stats(path):
    with open(path) as fin:
        for line in fin:
            line = line.strip().lstrip()
            if 'Number of scaffolds' in line:
                total_fragments = line.split()[0]
            elif 'Total length' in line:
                genome_size = line.split()[0]
            elif 'Scaffold N50' in line:
                n50 = line.split()[0]
                units = line.split()[1]
                if units == 'MB':
                    n50 = int(n50)*1000000
                elif units == 'KB':
                    n50 = int(n50)*1000

    return {"N50": n50, "bp": genome_size, "contigs": total_fragments}

def eukcc_parser( concoct, metabat ):
    eukcc_data = {}
    with open(concoct, 'r') as conco_in, open(metabat, 'r') as meta_in:
        next(conco_in)
        for line in conco_in:
            bin_id,completeness,contamination,ncbi_lng = line.strip().split("\t")
            eukcc_data[bin_id] = {
                    'bin_id' : bin_id+'.fa',
                    'completeness': completeness,
                    'contamination': contamination}

        next(meta_in)
        for line in meta_in:
            bin_id,completeness,contamination,ncbi_lng = line.strip().split("\t")
            eukcc_data[bin_id] = {
                    'bin_id' : bin_id+'.fa',
                    'completeness': completeness,
                    'contamination': contamination}

    return eukcc_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", 
        help="path for the output table", 
        default="qc.csv", 
        type=str
    )
    parser.add_argument(
        "--eukcc_c", 
        help="Eukcc results for concoct binner", 
        type=str,
        required=True
    )
    parser.add_argument(
        "--eukcc_m",
        help="Eukcc results for metabat binner",
        type=str,
        required=True
    )
    parser.add_argument(
        "--busco_files",
        help="List of busco outputs",
        nargs="*",
        required=True
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="rerun even if output exists",
        default=False,
    )
    parser.add_argument(
        "--quiet", 
        action="store_true", 
        help="supress information", 
        default=False
    )
    parser.add_argument(
        "--debug", 
        action="store_true", 
        help="Make it more verbose", 
        default=False
    )
    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S: ", level=logLevel
    )

    # for each mag learn eukcc and BUSCO values
    fields = [
        "fasta",
        "completeness",
        "contamination",
        "BUSCO_C",
        "BUSCO_M",
        "BUSCO_D",
        "BUSCO_S",
        "BUSCO_F",
        "BUSCO_n",
        "BUSCO_lineage",
        "N50",
        "bp",
        "contigs",
    ]

    genomes_list = []
    for file_in in args.busco_files:
        prefix = file_in.replace('.short_summary.specific.txt','') 
        genomes_list.append(prefix)

    eukcc_data = eukcc_parser( args.eukcc_c, args.eukcc_m )

    with open(args.output, "w") as outfile:
        cout = csv.DictWriter(outfile, fieldnames=fields, extrasaction="ignore")
        cout.writeheader()

        for mag in genomes_list:
            #eukcc_p = mag+".eukcc.csv"
            busco_p = mag+".short_summary.specific.txt"
            stat = genome_stats(busco_p)

            eukcc = eukcc_data[mag+'.fa']
            stat = {**stat, **eukcc}

            stat["fasta"] = mag+'.fa'
            busco = busco_score(busco_p)
            stat = {**stat, **busco}
            cout.writerow(stat)


