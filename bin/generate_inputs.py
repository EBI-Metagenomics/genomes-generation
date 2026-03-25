#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import re
import sys
import requests
import xmltodict
import json
import time
from typing import Tuple

from ena_portal_api.ena_handler import EnaApiHandler
handler = EnaApiHandler()

def parse_args():
    parser = argparse.ArgumentParser(description='The script generates input samplesheet for genomes-generation pipeline')
    parser.add_argument('-a', '--assembly-study', required=True,
                        help='Assembly study accession')
    parser.add_argument('-r', '--raw-reads-study', required=True,
                        help='Raw reads study accession')
    parser.add_argument('-o', '--outdir', required=False, default='outdir',
                        help='Output directory')
    parser.add_argument('--output-samplesheet', required=False,
                        help='Name of output file', default="samplesheet.csv")
    parser.add_argument('-b', '--scientific-name', required=False,
                        help='Comma separated list of scientific_name(s) to include. '
                             'For example: marine sediment metagenome,sediment metagenome')
    parser.add_argument('-e', '--environment-biome', required=False,
                        help='Comma separated list of environment_biome(s) to include. '
                             'For example: cold marine sediment biome,warm marine sediment biome')
    parser.add_argument('--keep-metat', action='store_true')
    return parser.parse_args()


def transform_paths(ftp_path: str) -> Tuple[str, str]:
    """
    Transforms an FTP path according to privacy, it also returns if it's public or private. There are some limitations.
    Use https:// for public (-resume nextflow works for https:// but not for ftp and s3)
    Use FIRE S3 object key for private data (-resume is broken for S3 paths)
    """
    if ftp_path.startswith("ftp.sra.ebi.ac.uk/vol1/"):
        # TODO: return s3 path when nextflow would be fixed
        # TODO s3_key = ftp_path.replace("ftp.sra.ebi.ac.uk/vol1/", "s3://era-public/")
        https_key = 'https://' + ftp_path
        print(f"Detected a public file for FTP path: {ftp_path}")
        return https_key, 'public'
    elif ftp_path.startswith("ftp.dcc-private.ebi.ac.uk/vol1/"):
        # TODO: return s3 path when nextflow would be fixed
        # TODO s3_key = ftp_path.replace("ftp.dcc-private.ebi.ac.uk/vol1/", "s3://era-private/")
        print(f"Detected a private file for FTP path: {ftp_path}")
        return ftp_path, 'private'
    else:
        raise ValueError(
            f"Invalid FTP path: {ftp_path}. Must start with 'ftp.sra.ebi.ac.uk/vol1/' or 'ftp.dcc-private.ebi.ac.uk/vol1/'."
        )


def load_xml(assembly, retry_attempts=1, retry_delay_min=5):
    xml_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{assembly}"

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
            print('Could not get the assembly software from ENA API. Leaving the field empty')
    return assembly_software


def get_assembly_data(assembly, assembly_path):
    # fetch assembly data
    assembler = generate_software_input(assembly)
    transformed_assembly = ""
    if assembly_path:
        # TODO: return s3 when nextflow will work with -resume
        transformed_assembly, _ = transform_paths(assembly_path)
    else:
        print(f'no assembly path {assembly}')
    return transformed_assembly, assembler


def get_run_data(run, run_path):
    # fetch runs data
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
        print(f'no runs path {run}')

    transformed_runs = []
    for item in chosen_run_path:
        # TODO: return when s3 will work with -resume
        transformed_path, _ = transform_paths(item)
        transformed_runs.append(transformed_path)
    runs_str = ","
    if len(transformed_runs) == 2:
        runs_str = f"{transformed_runs[0]},{transformed_runs[1]}"
    elif len(transformed_runs) == 1:
        runs_str = f"{transformed_runs[0]},"
    else:
        print(f"wrong length of run {run} paths: {','.join(transformed_runs)}")
    return runs_str


def write_samplesheet_line(assembly, assembly_path, run, run_path, samplesheet):
    transformed_assembly, assembler = get_assembly_data(assembly, assembly_path)
    runs_str = get_run_data(run, run_path)
    if not (transformed_assembly == "" or runs_str == ""):
        with open(samplesheet, 'a') as file_out:
            line = f"{run},{runs_str},{assembly},{transformed_assembly},{assembler}"
            file_out.write(line + '\n')
        return 1
    return 0


def get_run_from_assembly_xml(assembly):
    erz_acc = assembly['analysis_accession']
    json_analysis = load_xml(erz_acc)
    if not json_analysis:
        print(f"Could not load XML for assembly {erz_acc}")
        return None

    analysis = json_analysis.get("ANALYSIS_SET", {}).get("ANALYSIS", {})
    run_acc = None

    # 1. Try RUN_REF accession
    run_ref = analysis.get("RUN_REF")
    if run_ref:
        if isinstance(run_ref, list):
            run_acc = run_ref[0].get("@accession")
        elif isinstance(run_ref, dict):
            run_acc = run_ref.get("@accession")

    # 2. Try ANALYSIS_TYPE/SEQUENCE_ASSEMBLY/NAME
    if not run_acc:
        try:
            name = analysis["ANALYSIS_TYPE"]["SEQUENCE_ASSEMBLY"]["NAME"]
            match = re.search(r'[EDS]RR\d+', name)
            if match:
                run_acc = match.group(0)
        except (KeyError, TypeError):
            pass

    # 3. Try FILES/FILE @filename
    if not run_acc:
        try:
            files = analysis["FILES"]["FILE"]
            if isinstance(files, list):
                files = files[0]
            filename = files.get("@filename", "")
            match = re.search(r'[EDS]RR\d+', filename)
            if match:
                run_acc = match.group(0)
        except (KeyError, TypeError):
            pass

    if not run_acc or not run_acc.startswith(('ERR', 'DRR', 'SRR')):
        print(f"Invalid run name {run_acc} for assembly {erz_acc}")
        sys.exit(1)

    return run_acc


def fetch_data(raw_reads_study, assembly_study, outdir, input_scientific_name, input_env_biome, keep_metat):
    allowed_library_source = ['METAGENOMIC']
    run_fields = 'library_source,library_strategy,fastq_ftp'
    if input_scientific_name:
        input_scientific_name = [s.strip() for s in input_scientific_name.split(',')]
        run_fields += ',scientific_name'
    if input_env_biome:
        run_fields += ',environment_biome'
        input_env_biomes = [s.strip() for s in input_env_biome.split(',')]
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
            if run['environment_biome'] not in input_env_biomes:
                continue
        if keep_metat:
            allowed_library_source.append('METATRANSCRIPTOMIC')
        if run['library_source'] not in allowed_library_source:
            continue
        list_runs.append(run['run_accession'])
        run_paths[run['run_accession']] = run['fastq_ftp']
    print(f'Received {len(list_runs)} runs after filtering')
    assembly_run, assembly_path = {}, {}
    assemblies = handler.get_study_assemblies(assembly_study, fields='generated_ftp')
    print(f'Received {len(assemblies)} assemblies from {assembly_study}')

    assemblies_without_run = []
    with open(os.path.join(outdir, 'runs_assemblies.tsv'), 'w') as file_out:
        for assembly in assemblies:
            erz_acc = assembly['analysis_accession']
            retrieved_run = get_run_from_assembly_xml(assembly)
            if retrieved_run is None:
                assemblies_without_run.append(erz_acc)
                continue
            if retrieved_run not in list_runs:
                continue
            file_out.write('\t'.join([retrieved_run, erz_acc]) + '\n')
            assembly_run[erz_acc] = retrieved_run
            assembly_path[erz_acc] = assembly['generated_ftp']

    if assemblies_without_run:
        print(f"WARNING: {len(assemblies_without_run)} assemblies had no run detected: "
              f"{', '.join(assemblies_without_run)}")

    resolved = len(assembly_run)
    total = len(assemblies)
    print(f"Resolved runs for {resolved}/{total} assemblies")
    if resolved == 0:
        print("ERROR: no assemblies could be matched to a run. Aborting.")
        sys.exit(1)

    return assembly_run, assembly_path, run_paths


def generate_samplesheet(assembly_run, assembly_path, run_paths, samplesheet_path):
    written_lines = 0
    for assembly in assembly_run:
        written_lines += write_samplesheet_line(
            assembly=assembly,
            assembly_path=assembly_path[assembly],
            run=assembly_run[assembly],
            run_path=run_paths[assembly_run[assembly]].split(';'),
            samplesheet=samplesheet_path
        )
    return written_lines


def main(args):
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    # add header to input samplesheet
    samplesheet_path = os.path.join(args.outdir, args.output_samplesheet)
    with open(samplesheet_path, 'w') as file_out:
        file_out.write(','.join(['id', 'fastq_1', 'fastq_2', 'assembly_accession', 'assembly_fasta', 'assembler']) + '\n')

    # filter runs
    assembly_run_dict, assembly_path_dict, run_paths_dict = fetch_data(
        raw_reads_study=args.raw_reads_study,
        assembly_study=args.assembly_study,
        outdir=args.outdir,
        input_scientific_name=args.scientific_name,
        input_env_biome=args.environment_biome,
        keep_metat=args.keep_metat
    )
    written_lines = generate_samplesheet(assembly_run_dict, assembly_path_dict, run_paths_dict, samplesheet_path=samplesheet_path)
    print(f'Samplesheet is done for {written_lines} records')
    print('Done.')


if __name__ == '__main__':
    main(parse_args())