#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import sys


def main(assembly_dir, run_dir, rename_file, output):
    accessions = {}
    with open(rename_file, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')  # err,erz
            accessions[line[0]] = line[1]

    assemblies = {}
    for assembly_file in os.listdir(assembly_dir):
        accession = assembly_file.split('.')[0]
        assemblies[accession] = assembly_file

    runs = {}
    for run_file in os.listdir(run_dir):
        accession = run_file.split('.')[0].split('_')[0]
        runs.setdefault(accession, [])
        runs[accession].append(run_file)
    # double check for 2 fastq files or other strange cases
    for run_accession in runs:
        se, pe_forward, pe_reversed = False, False, False
        if len(runs[run_accession]) > 2:
            for raw_file in runs[run_accession]:
                if len(raw_file.split('_1')) >= 2:
                    pe_forward = raw_file
                elif len(raw_file.split('_2')) >= 2:
                    pe_reversed = raw_file
                else:
                    se = raw_file
            if se and pe_forward and pe_reversed:
                print(f"More than 3 raw files detected for {run_accession}, choosing PE reads")
                runs[run_accession] = [pe_forward, pe_reversed]

    with open(output, 'w') as file_out:
        file_out.write(','.join(['id', 'assembly', 'fastq_1', 'fastq_2', 'assembly_accession']) + '\n')
        for run_accession in accessions:
            list_values = []
            if run_accession not in accessions:
                print(f'No {run_accession} in fetched assemblies')
                continue
            assembly_accession = accessions[run_accession]
            list_values.extend([run_accession, os.path.join(assembly_dir, assemblies[run_accession])])
            if run_accession not in runs:
                print(f'No {run_accession} in fetched runs')
                continue
            if len(runs[run_accession]) == 2:
                if len(runs[run_accession][0].split('_1')) >= 2:
                    run1 = os.path.join(run_dir, runs[run_accession][0])
                    run2 = os.path.join(run_dir, runs[run_accession][1])
                else:
                    run1 = os.path.join(run_dir, runs[run_accession][1])
                    run2 = os.path.join(run_dir, runs[run_accession][0])
            else:
                run1 = os.path.join(run_dir, runs[run_accession][0])
                run2 = ""
            list_values.extend([run1, run2, assembly_accession])
            file_out.write(','.join(list_values) + '\n')

def parse_args():
    parser = argparse.ArgumentParser(description='The script creates a samplesheet input for genomes-generation pipeline')
    parser.add_argument('-a', '--assembly-dir', required=True,
                        help='Path to the directory containing downloaded assemblies (Assemblies/ACC/raw)')
    parser.add_argument('-r', '--run-dir', required=True,
                        help='Path to the directory containing downloaded raw reads (Raw_reads/ACC/raw)')
    parser.add_argument('-n', '--rename-file', required=True,
                        help='Path to the file with run-assembly correspondence information')
    parser.add_argument('-o', '--output-file', required=False,
                        help='Name of output file', default="samplesheet.tsv")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.assembly_dir, args.run_dir, args.rename_file, args.output_file)
