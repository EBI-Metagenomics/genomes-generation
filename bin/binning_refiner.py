#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com
# t.thomas@unsw.edu.au

# Binning_refiner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Binning_refiner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# metaWRAP author notes:
# I thank the original creator of this script! This is a great idea! To make
# this script more usable as part of the metaWRAP binning pipeline, I
# removed unnecessary visual aaspects of the original Bin_refiner script
# and made it python2 compatible.

# Check out the original program: https://github.com/songweizhi/Binning_refiner
# And the publication: https://www.ncbi.nlm.nih.gov/pubmed/28186226


import os
import glob
import shutil
import argparse
from time import sleep
from sys import stdout
from Bio import SeqIO
import shutil

SEPARATOR = '__'

def check_folder(directory):
    if not directory:
        return None
    if not os.path.exists(directory):
        print(f"No directory {directory}")
        return None
    if len(os.listdir(directory)) == 0:
        print(f"No data in {directory}")
        return None
    else:
        print(f"There are {len(os.listdir(directory))} bins in {directory}")
    if directory[-1] == '/':
        return os.path.abspath(directory[:-1])
    else:
        return os.path.abspath(directory)

##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-1',
                    required=True,
                    help='first bin folder name')

parser.add_argument('-2',
                    required=True,
                    help='second bin folder name')

parser.add_argument('-3',
                    required=False,
                    help='third bin folder name')

parser.add_argument('-o',
                    required=True,
                    help='output folder name')

parser.add_argument('-n',
                    required=True,
                    help='name prefix')

parser.add_argument('-ms',
                    required=False,
                    default=524288,
                    type=int,
                    help='(optional) minimum size for refined bins, default = 524288 (0.5Mbp)')

args = vars(parser.parse_args())

wd = os.getcwd()
output_folder = os.path.join(wd, args.get('o'))
output_name = args.get('n')
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)
print(f'Output folder {output_folder}')
refined_output_folder = os.path.join(output_folder, 'refined')
if not os.path.exists(refined_output_folder):
    os.mkdir(refined_output_folder)
print(f'Refined output {refined_output_folder}')

# check given folders
input_bin_folder_1 = check_folder(args.get('1'))
input_bin_folder_2 = check_folder(args.get('2'))
input_bin_folder_3 = check_folder(args.get('3'))
# get input bin folder list
input_bin_folder_list = [i for i in [input_bin_folder_1, input_bin_folder_2, input_bin_folder_3] if i]
print(f"Given {len(input_bin_folder_list)} not empty bin folders")
if len(input_bin_folder_list) == 1:
    print("No refinement needed")
    for bin in os.listdir(input_bin_folder_list[0]):
        shutil.copy(os.path.join(input_bin_folder_list[0], bin), os.path.join(refined_output_folder, bin))
    exit()

bin_size_cutoff = args['ms']
bin_size_cutoff_MB = float("{0:.2f}".format(bin_size_cutoff / (1024 * 1024)))

################################################ Define folder/file name ###############################################
########################################################################################################################

# check input files
bin_folder_bins_ext_list_uniq = set()
for bin_folder in input_bin_folder_list:
    bins = os.listdir(bin_folder)
    if len(bins) == 0:
        print('No input bin detected from %s folder, please double-check!' % (bin_folder))
        exit()
    for binfile in bins:
        extension = os.path.splitext(binfile)[1]
        bin_folder_bins_ext_list_uniq.add(extension)

    # check whether bins in the same folder have same extension, exit if not
    if len(bin_folder_bins_ext_list_uniq) > 1:
        print('Different file extensions were found from %s bins, please use same extension (fa, fas or fasta) '
              'for all bins in the same folder.' % (bin_folder))
        exit()

contig_bin_dict, contig_length_dict = {}, {}
combined_all_bins_file = os.path.join(output_folder, 'combined_all_bins.fasta')
print(f'Combining all bins together to {combined_all_bins_file}')
with open(combined_all_bins_file, 'w') as combined_file:
    for each_folder in input_bin_folder_list:
        each_folder_basename = os.path.basename(each_folder)
        sleep(1)
        print(f'Add folder/bin name to contig name for {each_folder} bins')
        for each_bin in os.listdir(each_folder):
            bin_file_name, bin_file_ext = os.path.splitext(each_bin)
            for each_contig in SeqIO.parse(os.path.join(each_folder, each_bin), 'fasta'):
                each_contig_new_id = SEPARATOR.join([each_folder_basename, bin_file_name, each_contig.id])
                if each_contig.id not in contig_bin_dict:
                    contig_bin_dict[each_contig.id] = [SEPARATOR.join([each_folder_basename, bin_file_name])]
                    contig_length_dict[each_contig.id] = len(each_contig.seq)
                elif each_contig.id in contig_bin_dict:
                    contig_bin_dict[each_contig.id].append(SEPARATOR.join([each_folder_basename, bin_file_name]))
                each_contig.id = each_contig_new_id
                each_contig.description = ''
                SeqIO.write(each_contig, combined_file, 'fasta')

contig_assignments_file = os.path.join(output_folder, 'contig_assignments.txt')
with open(contig_assignments_file, 'w') as contig_assignments:
    for each in contig_bin_dict:
        if len(contig_bin_dict[each]) == len(input_bin_folder_list):
            contig_assignments.write('\t'.join(contig_bin_dict[each] + [each, str(contig_length_dict[each])]) + '\n')

contig_assignments_file_sorted = os.path.join(output_folder, 'contig_assignments_sorted.txt')
contig_assignments_file_sorted_one_line = os.path.join(output_folder, 'contig_assignments_sorted_one_line.txt')
os.system(f'cat {contig_assignments_file} | sort > {contig_assignments_file_sorted}')

with open(contig_assignments_file_sorted, 'r') as contig_assignments_sorted, open(contig_assignments_file_sorted_one_line, 'w') as contig_assignments_sorted_one_line:
    current_match = ''
    current_match_contigs = []
    current_length_total = 0
    n = 1
    for each in contig_assignments_sorted:
        each_split = each.strip().split('\t')
        current_contig = each_split[-2]
        current_length = int(each_split[-1])
        matched_bins = '\t'.join(each_split[:-2])
        if current_match == '':
            current_match = matched_bins
            current_match_contigs.append(current_contig)
            current_length_total += current_length
        elif current_match == matched_bins:
            current_match_contigs.append(current_contig)
            current_length_total += current_length
        elif current_match != matched_bins:
            refined_bin_name = 'refined_bin%s' % n
            if current_length_total >= bin_size_cutoff:
                contig_assignments_sorted_one_line.write('Refined_%s\t%s\t%sbp\t%s\n' % (n, current_match, current_length_total,'\t'.join(current_match_contigs)))
                n += 1
            current_match = matched_bins
            current_match_contigs = []
            current_match_contigs.append(current_contig)
            current_length_total = 0
            current_length_total += current_length
    if current_length_total >= bin_size_cutoff:
        contig_assignments_sorted_one_line.write('Refined_%s\t%s\t%sbp\t%s\n' % (n, current_match, current_length_total,'\t'.join(current_match_contigs)))
    else:
        n -= 1

sleep(1)
print('The number of refined bins: %s' % n)

refined_bins = open(contig_assignments_file_sorted_one_line)

for each_refined_bin in refined_bins:
    each_refined_bin_split = each_refined_bin.strip().split('\t')
    each_refined_bin_name = output_name + '_' + each_refined_bin_split[0]
    if len(input_bin_folder_list) == 2:
        each_refined_bin_contig = each_refined_bin_split[4:]
    if len(input_bin_folder_list) == 3:
        each_refined_bin_contig = each_refined_bin_split[5:]

    added_contigs = []
    stdout.write(f'Extracting refined bin: {each_refined_bin_name}.fa')
    refined_bin_file = os.path.join(refined_output_folder, f'{each_refined_bin_name}.fa')
    with open(refined_bin_file, 'w') as refined_bin_handle:
        input_contigs = SeqIO.parse(combined_all_bins_file, 'fasta')
        for each_input_contig in input_contigs:
            each_input_contig_id = each_input_contig.id.split(SEPARATOR)[-1]
            if each_input_contig_id in each_refined_bin_contig:
                if each_input_contig_id not in added_contigs:
                    each_input_contig.id = each_input_contig_id
                    each_input_contig.description = ''
                    SeqIO.write(each_input_contig, refined_bin_handle, 'fasta')
                    added_contigs.append(each_input_contig_id)


