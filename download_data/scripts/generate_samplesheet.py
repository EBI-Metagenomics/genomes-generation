#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import sys

from ena_portal_api.ena_handler import EnaApiHandler
handler = EnaApiHandler()

def get_erz_list(input_dir):
    erz_list = set()
    for file in glob.glob(os.path.join(input_dir, 'ERZ*')):
        acc = file.strip().split('/')[-1].split('.')[0]
        erz_list.add(acc)
    return erz_list


def get_err_erz(erz_list, outfile_software):
    accessions = {}
    with open(outfile_software, 'w') as file_software:
        for erz_acc in erz_list:
            handler_request = handler.get_assembly(erz_acc)
            run_acc = handler_request["submitted_ftp"].strip().split('/')[-1].split('.')[0]
            assembly_software = handler_request["assembly_software"]
            if not run_acc.startswith(('ERR', 'DRR', 'SRR')):
                print('Invalid run name {} for assembly {}'.format(run_acc, erz_acc))
                sys.exit(1)
            file_software.write('\t'.join([run_acc, assembly_software]) + '\n')
            accessions[run_acc] = erz_acc
    return accessions


def main(assembly_dir, run_dir, output, output_assembly_software):
    erz_list = get_erz_list(assembly_dir)
    accessions = get_err_erz(erz_list, output_assembly_software)

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
            list_values.extend([run_accession, os.path.join(assembly_dir, assemblies[assembly_accession])])
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
    parser.add_argument('-o', '--output-file', required=False,
                        help='Name of output file', default="samplesheet.tsv")
    parser.add_argument('--assembly-software-filename', required=False,
                        help='Name of output file for assembly software information', default="assembly_software.tsv")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.assembly_dir, args.run_dir, args.output_file, args.assembly_software_filename)
