#!/usr/bin/env python3
import os
import argparse


##### This script generates the final output of taxonomy annotation for eukaryotic MAGs from BAT outputs
##### Alejandra Escobar, EMBL-EBI
##### July 4, 2023

def taxo_reporter( out_file, bat_files ):
    with open(out_file, 'w') as to_print, open('all_bin2classification.txt', 'w') as to_concat:
        to_print.write("\t".join(
            [
                'bin',
                'lineage_taxids',
                'lineage_names'
            ])+'\n')
        to_concat.write("\t".join(
            [
                '# bin',
                'classification',
                'reason',
                'lineage',
                'lineage scores'
            ])+'\n')
        for current_bin in bat_files:
            with open(current_bin, 'r') as file_in:
                next(file_in)
                for line in file_in:
                    line = line.rstrip().split("\t")
                    concat_list = []
                    genome = line.pop(0)
                    concat_list.append(genome)
                    classification = line.pop(0)
                    concat_list.append(classification)
                    reason = line.pop(0)
                    concat_list.append(reason)
                    lineage = line.pop(0)
                    concat_list.append(lineage)
                    lineage_scores = line.pop(0)
                    concat_list.append(lineage_scores)
                    to_concat.write("\t".join(concat_list)+'\n')
                    clean_lineage = []
                    for tax_rank in line:
                        name_rank = tax_rank.split(':')[0].split(' (')[0]
                        clean_lineage.append(name_rank)
                    clean_lineage = ";".join(clean_lineage)
                    to_print.write("\t".join([genome,lineage,clean_lineage])+'\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script generates the final output of taxonomy annotation for eukaryotic MAGs from BAT outputs"
    )
    parser.add_argument(
        "--output", 
        help="path for the output table", 
        default="euk_taxonomy.csv", 
        type=str
    )
    parser.add_argument(
        "--bat_names", 
        help="BAT names files (*.BAT_run.bin2classification.names.txt)", 
        nargs="*",
        required=True
    )
    args = parser.parse_args()

    taxo_reporter( args.output, args.bat_names )
