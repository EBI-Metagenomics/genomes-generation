#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import sys

from ena_portal_api.ena_handler import EnaApiHandler
handler = EnaApiHandler()

def parse_args():
    parser = argparse.ArgumentParser(description='The script filters list of assemblies')
    parser.add_argument('-a', '--assembly-study', required=True,
                        help='Assembly study accession')
    parser.add_argument('-r', '--raw-reads-study', required=True,
                        help='Raw reads study accession')
    parser.add_argument('-o', '--outdir', required=False, default='outdir',
                        help='Output directory')
    parser.add_argument('--keep-metat', action='store_true')
    parser.add_argument('--scientific-name', required=False, nargs="+",
                        help='Only runs with specified scientific_name(s) would be included')
    return parser.parse_args()


def main(args):
    allowed_library_source = ['METAGENOMIC']
    scientific_name = []
    # if argument came from bash script it is interpreted as list of one item
    if ',' in args.scientific_name[0]:
        scientific_name = args.scientific_name[0].split(',')
        scientific_name = [s.strip() for s in scientific_name]
    # if argument was used directly in script it is correct list of items
    elif args.scientific_name:
        scientific_name = args.scientific_name
    raw_reads_study = args.raw_reads_study
    runs = handler.get_study_runs(raw_reads_study, fields='scientific_name,library_source,library_strategy')
    list_runs = []
    for run in runs:
        if run['library_strategy'] == 'AMPLICON':
            continue
        if scientific_name:
            if run['scientific_name'] not in scientific_name:
                continue
        if args.keep_metat:
            allowed_library_source.append('METATRANSCRIPTOMIC')
        if run['library_source'] not in allowed_library_source:
            continue
        list_runs.append(run['run_accession'])

    assembly_study = args.assembly_study
    assemblies = handler.get_study_assemblies(assembly_study, fields='submitted_ftp')
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    with open(os.path.join(args.outdir, 'runs_assemblies.tsv'), 'w') as file_out:
        for assembly in assemblies:
            retrieved_run = assembly['submitted_ftp'].split('/')[-1].split('.')[0]
            if not retrieved_run.startswith(('ERR', 'DRR', 'SRR')):
                print('Invalid run name {} for assembly {}'.format(run_acc, erz_acc))
                sys.exit(1)
            if retrieved_run not in list_runs:
                continue
            file_out.write('\t'.join([retrieved_run, assembly['analysis_accession']]) + '\n')


if __name__ == '__main__':
    main(parse_args())