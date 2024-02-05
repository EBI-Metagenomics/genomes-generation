#!/usr/bin/env python3
# coding=utf-8

import argparse
import logging
import os
import sys
from Bio import SeqIO
from shutil import copy, rmtree
import numpy as np
import math
from itertools import combinations


def calculate_n50(seq_lens):
    seq_array = np.array(seq_lens)
    sorted_lens = seq_array[np.argsort(-seq_array)]
    csum = np.cumsum(sorted_lens)
    level = 50
    nx = int(int(np.sum(sorted_lens)) * (level / 100))
    csumn = min(csum[csum >= nx])
    l_level = int(np.where(csum == csumn)[0])
    n_level = int(sorted_lens[l_level])
    return n_level


def get_stats(binner, stats):
    """Get the stats file for a binner.
    The stats folder contains the stats for each binner, and also the
    combination of binners (such as binner1, binner12...)
    """
    # We need to add the prefix _ at the end of the binner name
    # to be sure we pick the right binner stat file
    # the binner stats files are named, i.e.: binner2_filtered_genomes.tsv binner23_filtered_genomes.tsv
    # for binner2 we may get binner23 if not prefixed with _
    binner = f"{binner}_"
    stats_dict = {}
    stats_file = None
    logging.debug(f"{binner}")
    for item in os.listdir(stats):
        if binner in item:
            stats_file = item
            break
    if not stats_file:
        logging.error("get_stats: No stats file. Exit")
        sys.exit(1)
    logging.debug(f"get_stats: Process {stats_file}")
    with open(os.path.join(stats, stats_file)) as file_in:
        next(file_in)
        for line in file_in:
            line = line.strip().split("\t")
            stats_dict[line[0] + ".fa"] = [float(line[1]), float(line[2])]
    return stats_dict


def get_bins_from_binner(binner_dir):
    binner_dict = {}
    for bin_name in os.listdir(binner_dir):
        bin_path = os.path.join(binner_dir, bin_name)
        binner_dict[bin_name] = {}
        logging.debug(f"get_bins_from_binner: Add {bin_path}")
        with open(bin_path) as file_in:
            for record in SeqIO.parse(file_in, "fasta"):
                binner_dict[bin_name][record.id] = len(record.seq)
    return binner_dict


def process_pair(binner1, binner2, stats_path, bins_2_stats):
    if not binner2:
        logging.info('-'.join(['-']*20) + f'> Add {binner1}')
        return get_bins_from_binner(binner1), get_stats(binner1, stats_path)
    else:
        logging.info('-'.join(['-']*20) + f'> Add {binner1}')
        bins1 = get_bins_from_binner(binner1)
        bins2 = binner2
        all_bin_pairs = {}
        for bin_1 in bins1:
            all_bin_pairs[bin_1] = {}
            for bin_2 in bins2:
                # find identical contigs between bin_1 and bin_2
                match_1_length, match_2_length, mismatch_1_length, mismatch_2_length = [0 for _ in range(4)]
                for contig in bins1[bin_1]:
                    if contig in bins2[bin_2]:
                        match_1_length += bins2[bin_2][contig]
                    else:
                        mismatch_1_length += bins1[bin_1][contig]
                for contig in bins2[bin_2]:
                    if contig in bins1[bin_1]:
                        match_2_length += bins1[bin_1][contig]
                    else:
                        mismatch_2_length += bins2[bin_2][contig]

                logging.debug(f'compare {bin_1} vs {bin_2}')
                # chose the highest % ID, dependinsh of which bin is a subset of the other
                ratio_1 = 100 * match_1_length / (match_1_length + mismatch_1_length)
                ratio_2 = 100 * match_2_length / (match_2_length + mismatch_2_length)
                logging.debug(f'ratio: {max([ratio_1, ratio_2])}')
                if max([ratio_1, ratio_2]) >= 80:
                    all_bin_pairs[bin_1][bin_2] = max([ratio_1, ratio_2])
            if all_bin_pairs[bin_1] == {}:
                all_bin_pairs.pop(bin_1)

        if not all_bin_pairs:
            logging.info('No ratios > 80')
            return bins2, bins_2_stats
        logging.debug(f'all_bin_pair ratios {all_bin_pairs}')

        logging.info("Choose best bin in pair")
        # choose bins
        best_bins, best_bins_stats, best_bins_contigs = {}, {}, {}
        bins_2_matches = []

        bins_1_stats = get_stats(binner1, stats_path)
        for bin_1 in all_bin_pairs:
            N50 = calculate_n50(list(bins1[bin_1].values()))
            score = bins_1_stats[bin_1][0] - bins_1_stats[bin_1][1] * 5 + 0.5 * math.log(N50)
            current_bin = bin_1
            current_contigs = bins1[bin_1]
            current_stats = bins_1_stats[bin_1]

            for bin_2 in all_bin_pairs[bin_1]:
                logging.debug(f'compare: >>> {bin_1} with {score}')
                # check for sufficient overlap (80% bin length)
                bins_2_matches.append(bin_2)
                # check if this bin is better than original
                N50 = calculate_n50(list(bins2[bin_2].values()))
                if (bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * 5 + 0.5 * math.log(N50)) > score:
                    current_bin = bin_2
                    current_contigs = bins2[bin_2]
                    current_stats = bins_2_stats[bin_2]
                    logging.debug(
                        f'compare: {bin_2} better then {bin_1} {bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * 5} > {score}')
                elif (bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * 5) == score:
                    logging.debug(
                        f'compare: {bin_2} the same as {bin_1} {bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * 5} = {score}')
                    if bins_1_stats[bin_1][0] > bins_2_stats[bin_2][0] and bins_1_stats[bin_1][1] < bins_2_stats[bin_2][1]:
                        logging.debug(
                            f'compare: {bin_1} better {bin_2} by completeness and contamination')
                    elif bins_1_stats[bin_1][0] == bins_2_stats[bin_2][0] and bins_1_stats[bin_1][1] == bins_2_stats[bin_2][1]:
                        logging.debug(
                            f'compare: {bin_1} [{len(bins1[bin_1])} contigs] the same as {bin_2} [{len(bins2[bin_2])} contigs] by completeness and contamination')
                    elif bins_1_stats[bin_1][0] < bins_2_stats[bin_2][0] and bins_1_stats[bin_1][1] > bins_2_stats[bin_2][1]:
                        logging.debug(
                            f'compare: {bin_2} better {bin_1} by completeness and contamination')
                    else:
                        logging.debug(
                            f'compare: difficult case {bin_1} [{len(bins1[bin_1])} contigs] ({bins_1_stats[bin_1][0]}, {bins_1_stats[bin_1][1]}) vs {bin_2} [{len(bins2[bin_2])} contigs] ({bins_2_stats[bin_2][0]}, {bins_2_stats[bin_2][1]})')
                else:
                    logging.debug(f'compare: {bin_1} wins')
            best_bins_contigs[current_bin] = current_contigs
            best_bins_stats[current_bin] = current_stats

        # retrieve bins from second group that were not found in first group
        untouched_bins2 = list(set(bins_2_stats.keys()).difference(set(bins_2_matches)))
        for bin_2 in untouched_bins2:
            best_bins_stats[bin_2] = bins_2_stats[bin_2]
            best_bins_contigs[bin_2] = bins2[bin_2]
        return best_bins_contigs, best_bins_stats


def parse_args():
    parser = argparse.ArgumentParser(description='The script creates a file that matches ERZ and read accessions')
    parser.add_argument('-i', '--input', required=True, nargs='+', help='Folders with bins')
    parser.add_argument('-s', '--stats', required=True, help='Folders with checkm stats')
    parser.add_argument('-v', '--verbose', action="store_true")
    return parser.parse_args()


def main(args):
    consolidated_bins = "consolidated_bins"
    dereplicated_bins = "dereplicated_bins"
    input_binners = args.input
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    binners = []
    for binner in input_binners:
        if os.path.exists(binner):
            if os.listdir(binner):
                binners.append(binner)
    if len(binners) < 2:
        logging.info('Number of binners is less then 2. No consolidation')
        if len(binners) == 1:
            logging.info(f'Dereplicated bins would be taken from {binners[0]}')
            if not os.path.exists(dereplicated_bins):
                os.mkdir(dereplicated_bins)
            with open("dereplicated_list.tsv", 'w') as file_out:
                for item in os.listdir(binners[0]):
                    copy(os.path.join(binners[0], item), os.path.join(dereplicated_bins, item))
                    file_out.write(item + '\n')
        sys.exit()
    else:
        logging.info('-'.join(['-']*20) + f'---> Processing {len(binners)} input binners')
        best_bins, best_stats = {}, {}
        for binner in binners:
            bins, stats = process_pair(binner, best_bins, args.stats, best_stats)
            best_bins = bins
            best_stats = stats
            logging.info(f'best set has {len(best_stats)} bins')
            for bin in bins:
                logging.debug(f'{bin} has {len(list(bins[bin].keys()))} contigs')
            logging.debug(f'stats {best_stats}')

    # consolidate folder and stats
    if not os.path.exists(consolidated_bins):
        os.mkdir(consolidated_bins)
    else:
        rmtree(consolidated_bins)
        os.mkdir(consolidated_bins)
    for bin in best_bins:
        for binner in binners:
            if bin in os.listdir(binner):
                copy(os.path.join(binner, bin), os.path.join(consolidated_bins, bin))
                break

    # write consolidated stats
    with open("consolidated_stats.tsv", 'w') as file_out:
        file_out.write("\t".join(['bin', 'completeness', 'contamination']) + '\n')
        for item in best_stats:
            file_out.write("\t".join([item, str(best_stats[item][0]), str(best_stats[item][1])]) + '\n')

    # generate dereplicated bins
    contig_mapping = {}
    for bin in os.listdir(consolidated_bins):
        bin_path = os.path.join(consolidated_bins, bin)
        with open(bin_path, 'r') as file_in:
            for record in SeqIO.parse(file_in, "fasta"):
                if record.id not in contig_mapping:
                    contig_mapping[record.id] = bin
                elif best_stats[bin] > best_stats[contig_mapping[record.id]]:
                    contig_mapping[record.id] = bin

    if not os.path.exists(dereplicated_bins):
        os.mkdir(dereplicated_bins)
    for bin in os.listdir(consolidated_bins):
        bin_path = os.path.join(consolidated_bins, bin)
        drep_bin_path = os.path.join(dereplicated_bins, bin)
        with open(bin_path, 'r') as file_in, open(drep_bin_path, 'w') as file_out:
            for record in SeqIO.parse(file_in, "fasta"):
                if contig_mapping[record.id] == bin:
                    SeqIO.write(record, file_out, "fasta")

    # remove empty files
    for item in os.listdir(dereplicated_bins):
        if not os.stat(os.path.join(dereplicated_bins, item)).st_size:
            os.remove(os.path.join(dereplicated_bins, item))

    # write dereplicated_list.tsv
    with open("dereplicated_list.tsv", 'w') as file_out:
        drep_bins = os.listdir(dereplicated_bins)
        file_out.write('\n'.join(drep_bins))

if __name__ == '__main__':
    args = parse_args()
    main(args)