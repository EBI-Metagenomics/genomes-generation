#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import sys
import requests
import xmltodict
import json
import time
from typing import Tuple

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
    parser.add_argument('--output-samplesheet', required=False,
                        help='Name of output file', default="samplesheet.csv")
    parser.add_argument('-b', '--scientific-name', required=False,
                        help='Comma separated list of scientific_name(s) that would be included')
    parser.add_argument('-e', '--environment-biome', required=False,
                        help='Comma separated list of environment_biome(s) that would be included')
    parser.add_argument('--keep-metat', action='store_true')
    return parser.parse_args()


def transform_paths(ftp_path: str) -> Tuple[str, str]:
    """
    Transforms an FTP path according to privacy, it also returns if it's public or private. There are some limitations.
    Use https:// for public (-resume nextflow works for https:// but not for ftp and s3)
    Use FIRE S3 object key for private data (-resume is broken for S3 paths)
    """
    if ftp_path.startswith("ftp.sra.ebi.ac.uk/vol1/"):
        #s3_key = ftp_path.replace("ftp.sra.ebi.ac.uk/vol1/", "s3://era-public/")
        https_key = 'https://' + ftp_path
        print(f"Detected a public file for FTP path: {ftp_path}")
        return https_key, 'public'
    elif ftp_path.startswith("ftp.dcc-private.ebi.ac.uk/vol1/"):
        s3_key = ftp_path.replace("ftp.dcc-private.ebi.ac.uk/vol1/", "s3://era-private/")
        print(f"Detected a private file for FTP path: {ftp_path}")
        return s3_key, 'private'
    else:
        raise ValueError(
            f"Invalid FTP path: {ftp_path}. Must start with 'ftp.sra.ebi.ac.uk/vol1/' or 'ftp.dcc-private.ebi.ac.uk/vol1/'."
        )


def load_xml(assembly):
    retry_attempts = 3
    retry_delay_min = 15
    xml_url = "https://www.ebi.ac.uk/ena/browser/api/xml/{}".format(assembly)

    for attempt in range(1, retry_attempts + 1):
        r = requests.get(xml_url)
        if not r.ok:
            retry_delay = retry_delay_min * attempt
            if r.status_code == 500:
                print(f"Received HTTP 500 error for sample {assembly}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
                continue
            else:
                print(f"Unable to fetch xml for sample {assembly}. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
                continue
        if r.ok:
            data_dict = xmltodict.parse(r.content)
            json_dump = json.dumps(data_dict)
            json_data = json.loads(json_dump)
            return json_data
        else:
            print("Could not retrieve xml for accession", assembly)
            print(r.text)
            return None


def generate_software_input(erz_acc):
    handler_request = handler.get_assembly(erz_acc)
    assembly_software = handler_request["assembly_software"]
    if not assembly_software:
        try:
            json_analysis = load_xml(erz_acc)
            assembly_software = json_analysis["ANALYSIS_SET"]["ANALYSIS"]["ANALYSIS_TYPE"]["SEQUENCE_ASSEMBLY"][
                "PROGRAM"]
            assembly_software = assembly_software.replace(" ", "_" if "v" in assembly_software else "_v")
            if '_v' not in assembly_software:
                assembly_software = assembly_software.replace('v', '_v')
        except:
            print('Did not get software from XML. Leaving empty field')
    return assembly_software


def write_samplesheet_line(assembly, assembly_path, run, run_path, samplesheet):
    with open(samplesheet, 'a') as file_out:
        transformed_assembly = ""
        if assembly_path:
            # TODO: return s3 when nextflow will work with -resume
            transformed_assembly, _ = transform_paths(assembly_path)
        else:
            print(f'no assembly path {assembly}')
        chosen_run_path = []
        if run_path:
            se, pe_forward, pe_reversed = False, False, False
            if len(run_path) >= 2:
                for raw_file in run_path:
                    if len(raw_file.split('_1')) >= 2:
                        pe_forward = raw_file
                    elif len(raw_file.split('_2')) >= 2:
                        pe_reversed = raw_file
                    else:
                        se = raw_file
                if se and pe_forward and pe_reversed:
                    print(f"More than 3 raw files detected for {run}, choosing PE reads")
                    chosen_run_path = [pe_forward, pe_reversed]
                elif pe_forward and pe_reversed:
                    chosen_run_path = [pe_forward, pe_reversed]
                else:
                    chosen_run_path = [se]
            else:
                chosen_run_path = [run_path[0]]
        else:
            print(f'no runs path {assembly}')

        transformed_runs = []
        for item in chosen_run_path:
            # TODO: return when s3 will work with -resume
            transformed_path, _ = transform_paths(item)
            transformed_runs.append(transformed_path)

        line = f"{run},"
        if len(transformed_runs) == 2:
            line += f"{transformed_runs[0]},{transformed_runs[1]},"
        elif len(transformed_runs) == 1:
            line += f"{transformed_runs[0]},,"
        else:
            print("wrong length of runs path")
        assembler = generate_software_input(assembly)
        line += f"{assembly},{transformed_assembly},{assembler}"
        if transformed_assembly == '' or transformed_runs == []:
            print(f'Some data is missing for {run} - {assembly}. Skipping that raw in samplesheet.')
        else:
            file_out.write(line + '\n')


def generate_lists(raw_reads_study, assembly_study, outdir, input_scientific_name, input_env_biome, keep_metat, samplesheet_path):
    allowed_library_source = ['METAGENOMIC']
    run_fields = 'library_source,library_strategy,fastq_ftp'
    if input_scientific_name:
        input_scientific_name = [s.strip() for s in input_scientific_name.split(',')]
        run_fields += ',scientific_name'
    if input_env_biome:
        run_fields += ',environment_biome'
        input_input_env_biomes = [s.strip() for s in input_env_biome.split(',')]
    runs = handler.get_study_runs(raw_reads_study, fields=run_fields)
    print(f'Received {len(runs)} runs from {raw_reads_study}')
    list_runs = []
    run_paths = {}
    for run in runs:
        if run['library_strategy'] == 'AMPLICON':
            continue
        if input_scientific_name:
            if run['scientific_name'] not in input_scientific_name:
                continue
        if input_env_biome:
            if run['environment_biome'] not in input_input_env_biomes:
                continue
        if keep_metat:
            allowed_library_source.append('METATRANSCRIPTOMIC')
        if run['library_source'] not in allowed_library_source:
            continue
        list_runs.append(run['run_accession'])
        run_paths[run['run_accession']] = run['fastq_ftp']
    print(f'Received {len(list_runs)} runs after filtering')
    assembly_run = {}
    assemblies = handler.get_study_assemblies(assembly_study, fields='submitted_ftp,generated_ftp')
    print(f'Received {len(assemblies)} assemblies from {assembly_study}')
    with open(os.path.join(outdir, 'runs_assemblies.tsv'), 'w') as file_out:
        for assembly in assemblies:
            retrieved_run = assembly['submitted_ftp'].split('/')[-1].split('.')[0].split('_')[0]
            if not retrieved_run.startswith(('ERR', 'DRR', 'SRR')):
                print('Invalid run name {} for assembly {}'.format(run_acc, erz_acc))
                sys.exit(1)
            if retrieved_run not in list_runs:
                continue
            file_out.write('\t'.join([retrieved_run, assembly['analysis_accession']]) + '\n')
            assembly_run[assembly['analysis_accession']] = retrieved_run

            write_samplesheet_line(
                assembly=assembly['analysis_accession'],
                assembly_path=assembly['generated_ftp'],
                run=retrieved_run,
                run_path=run_paths[retrieved_run].split(';'),
                samplesheet=samplesheet_path
            )
    return assembly_run


def main(args):
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    # add header to input samplesheet
    samplesheet_path = os.path.join(args.outdir, args.output_samplesheet)
    with open(samplesheet_path, 'w') as file_out:
        file_out.write(','.join(['id', 'fastq_1', 'fastq_2', 'assembly_accession', 'assembly_fasta', 'assembler']) + '\n')

    # filter runs
    assembly_run_dict = generate_lists(
        raw_reads_study=args.raw_reads_study,
        assembly_study=args.assembly_study,
        outdir=args.outdir,
        input_scientific_name=args.scientific_name,
        input_env_biome=args.environment_biome,
        keep_metat=args.keep_metat,
        samplesheet_path=samplesheet_path
    )
    print(f'Samplesheet is done for {len(assembly_run_dict)} records')
    print('Done.')


if __name__ == '__main__':
    main(parse_args())