#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from Bio import SeqIO
import sys
import argparse
import os
import subprocess
import glob
import gzip

COVERAGE_FOLDER_NAME = 'coverage'

def process_genomes(genomes, name):
    ### Generating the contigs2bins.txt file
    ## Getting the list of bins

    files_list = os.listdir(genomes)
    bin_list = []
    for element in files_list:
        suffix = element.split('.')[-1]
        if suffix == 'fa':
            bin_list.append(element)

    print('Looking for your genomes...', file=sys.stdout)
    if len(bin_list) == 0:
        sys.exit("No genomes were found, please ensure that your assemblies are not compressed")

    with open(os.path.join(COVERAGE_FOLDER_NAME, name + "_contigs2bins.txt"), "w") as contigs2bins_out:
        bin_contigs = {}
        for bin_file in bin_list:
            bin_prefix = bin_file.replace('.fa', '')
            bin_contigs[bin_prefix] = []
            for record in SeqIO.parse(genomes + '/' + bin_file, "fasta"):
                bin_contigs[bin_prefix].append(str(record.id))
                contigs2bins_out.write(bin_prefix + '\t' + str(record.id) + '\n')
    return bin_list, bin_contigs

def process_metabat_depth(metabat_depth):
    ### Recovering the coverage files per run
    print('Looking for your metabat_depth.txt files...', file=sys.stdout)
    coverageFiles = metabat_depth
    if not coverageFiles:
        sys.exit(
            "No coverage files were found, please double check you have at least one metabat_depth.txt file in the provided path")
    ### Saving the contigs length and total average depth
    contig_cov = {}
    for cov_path in coverageFiles:
        with open(cov_path, 'r') as curr_cov:
            for line in curr_cov:
                if not line.startswith('contigName'):
                    contigName, contigLen, totalAvgDepth, bam, bamvar = line.strip().split('\t')
                    values = (contigLen, totalAvgDepth, bam, bamvar)
                    contig_cov[contigName] = values
    return contig_cov


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script reformat the coverage files generated for binning to be used in MAGs uploading"
    )
    parser.add_argument(
        "-g", "--genomes", dest="genomes", help="/full/path/to/dereplicated_genomes", required=True
    )
    parser.add_argument(
        "-n", "--name", dest="name", help="Run accession", required=True
    )
    parser.add_argument(
        "-m", "--metabat-depth", dest="metabat_depth", help="metabat_depth.txt", required=True, nargs='+'
    )
    args = parser.parse_args()
    if not os.path.exists(COVERAGE_FOLDER_NAME):
        os.mkdir(COVERAGE_FOLDER_NAME)

    bin_list, bin_contigs = process_genomes(args.genomes, args.name)
    contig_cov = process_metabat_depth(args.metabat_depth)

    ### Generating the output directories with their respective average coverage and coverage.tab files
    print('Creating coverage files per MAG...', file=sys.stdout)
    cov_tab_header = 'contigName\tcontigLen\ttotalAvgDepth\tsorted.bam\tsorted.bam-var'
    for final_bin in bin_list:
        assembled_pairs = 0
        assembly_length = 0
        primary_dir = os.path.join(COVERAGE_FOLDER_NAME, final_bin+'_coverage')
        if not os.path.exists(primary_dir):
            os.mkdir(primary_dir)
        cov_ave_out = os.path.join(primary_dir, final_bin+'_MAGcoverage.txt')
        cov_table = os.path.join(primary_dir, 'coverage.tab')

        with open(cov_table, "w") as large_out, open(cov_ave_out, "w") as short_out:
            large_out.write(cov_tab_header+'\n')

            for contig in bin_contigs[final_bin.replace('.fa','')]:
                large_out.write(contig+'\t'+'\t'.join(contig_cov[contig])+'\n')
                length = float(contig_cov[contig][0])
                read_depth = float(contig_cov[contig][1])
                assembled_pairs += (length * read_depth)
                assembly_length += length

            cov = assembled_pairs/assembly_length
            short_out.write(final_bin+'\t'+str(cov)+'\n')

    print('Done!', file=sys.stdout)