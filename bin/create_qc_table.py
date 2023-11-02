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


def busco_stats(path):
    d = {
        "BUSCO_completeness": None,
        "BUSCO_contamination": None,
        "BUSCO_C": None,
        "BUSCO_M": None,
        "BUSCO_D": None,
        "BUSCO_S": None,
        "BUSCO_F": None,
        "N50": None,
        "bp": None,
        "contigs": None
    }

    with open(path) as fin:
        for line in fin:
            line = line.strip().lstrip()
            if 'Number of scaffolds' in line:
                total_fragments = line.split()[0]
                d["contigs"] = str(total_fragments)
            elif 'Total length' in line:
                genome_size = line.split()[0]
                d["bp"] = str(genome_size)
            elif 'Scaffold N50' in line:
                n50 = line.split()[0]
                units = line.split()[1]
                if units == 'MB':
                    n50 = int(n50)*1000000
                elif units == 'KB':
                    n50 = int(n50)*1000
                d["N50"] = str(n50)
            elif line.startswith("C:"):
                elem = line.strip().replace("[", ",").replace("]", "").split(",")
                for e in elem:
                    x = e.split(":")
                    d["BUSCO_{}".format(x[0])] = x[1].replace("%", "")
            elif line.startswith("# The lineage dataset is: "):
                d["BUSCO_lineage"] = str(line.strip().split()[5])
        d["BUSCO_completeness"] = str(100 - float(d["BUSCO_M"]))
        d["BUSCO_contamination"] = d["BUSCO_D"]
    return d


def eukcc_parser(eukcc_concat, genomes_list, output_file):
    eukcc_data = {}
    with open(eukcc_concat, 'r') as file_in:
        logging.info("Processing eukcc stats")
        for line in file_in:
            if "completeness" not in line:
                genome,completeness,contamination = line.strip().split(",")
                eukcc_data[genome] = {
                        'bin_id' : genome,
                        'completeness': completeness,
                        'contamination': contamination}
        logging.info(f"Eukcc stats has {len(eukcc_data)} records")
    with open(output_file, 'w') as file_out:
        logging.info("Writing eukcc stats")
        file_out.write("genome,completeness,contamination" + '\n')
        for item in eukcc_data:
            if item in genomes_list:
                file_out.write(','.join([eukcc_data[item]['bin_id'], eukcc_data[item]['completeness'], eukcc_data[item]['contamination']]) + '\n')
        logging.info("Writing eukcc stats DONE")
    return eukcc_data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", 
        help="path for the output table", 
        default="qc.csv", 
        type=str
    )
    parser.add_argument(
        "--output_busco",
        help="path for the busco output table",
        default="busco_final_qc.csv",
        type=str
    )
    parser.add_argument(
        "--output_eukcc",
        help="path for the eukccoutput table",
        default="eukcc_final_qc.csv",
        type=str
    )
    parser.add_argument(
        "--eukcc_concat", 
        help="Eukcc combined results for QC50 (euk_quality.csv)", 
        type=str,
        required=False
    )
    parser.add_argument(
        "--busco_files",
        help="List of busco outputs",
        nargs="*",
        required=False
    )
    parser.add_argument(
        "--genomes_list",
        help="List of genomes to output stats for",
        type=str,
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
    # check inputs
    if not(args.eukcc_concat or args.busco_files):
        logging.error("No necessary files provided")
        exit(1)

    # create list of genomes
    genomes_list = []
    with open(args.genomes_list, 'r') as file_in:
        for line in file_in:
            genomes_list.append(line.strip())
        logging.info(f"There are {len(genomes_list)} genomes")

    # parse eukcc tables
    eukcc_data = {}
    if args.eukcc_concat:
        eukcc_data = eukcc_parser(args.eukcc_concat, genomes_list, args.output_eukcc)

    busco_data = {}
    if args.busco_files:
        logging.info("Processing busco stats")
        with open(args.output_busco, 'w') as file_out:
            file_out.write("genome,completeness,contamination,BUSCO_C,BUSCO_M,BUSCO_D,BUSCO_S,BUSCO_F,BUSCO_n,BUSCO_lineage,N50,bp,contigs" + '\n')
            for input_file in args.busco_files:
                logging.info(f"Processing {input_file}")
                stats = busco_stats(input_file)
                prefix = input_file.replace('.short_summary.specific.txt', '') + '.fa'
                busco_data[prefix] = stats
                if prefix in genomes_list:
                    file_out.write(",".join([prefix, stats['BUSCO_completeness'], stats['BUSCO_contamination'],
                                             stats['BUSCO_C'], stats['BUSCO_M'], stats['BUSCO_D'], stats['BUSCO_S'],
                                             stats['BUSCO_F'], stats['BUSCO_n'], stats['BUSCO_lineage'], stats['N50'],
                                             stats['bp'], stats['contigs']]) + '\n')
        logging.info("Processing busco stats DONE")

    if args.busco_files and args.eukcc_concat:
        logging.info("Writing final output")
        with open(args.output, "w") as outfile:
            outfile.write("genome,eukcc_completeness,eukcc_contamination,BUSCO_completeness,BUSCO_contamination,BUSCO_C,BUSCO_M,BUSCO_D,BUSCO_S,BUSCO_F,BUSCO_n,BUSCO_lineage,N50,bp,contigs" + '\n')
            for genome in genomes_list:
                if genome in eukcc_data and genome in busco_data:
                    stats = busco_data[genome]
                    outfile.write(
                        ",".join([genome, eukcc_data[genome]['completeness'], eukcc_data[genome]['contamination'],
                                  stats['BUSCO_completeness'], stats['BUSCO_contamination'],
                                  stats['BUSCO_C'], stats['BUSCO_M'], stats['BUSCO_D'], stats['BUSCO_S'],
                                  stats['BUSCO_F'], stats['BUSCO_n'], stats['BUSCO_lineage'], stats['N50'],
                                  stats['bp'], stats['contigs']]) + '\n')

if __name__ == "__main__":
    main()