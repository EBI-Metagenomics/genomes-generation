#!/usr/bin/env python3
# coding=utf-8

import argparse
import glob
import os
import sys

from ena_portal_api.ena_handler import EnaApiHandler
handler = EnaApiHandler()


def main(input_dir, outfile):
    erz_list = get_erz_list(input_dir)
    print('Obtained a list of ERZ accessions, generating renaming file...')
    with open(outfile, 'w') as file_out:
        for erz_acc in erz_list:
            ftp_loc = handler.get_assembly(erz_acc)["submitted_ftp"]
            run_acc = ftp_loc.strip().split('/')[-1].split('.')[0]
            if not run_acc.startswith(('ERR','DRR','SRR')):
                print('Invalid run name {} for assembly {}'.format(run_acc, erz_acc))
                sys.exit(1)
            file_out.write(','.join([run_acc, erz_acc]) + '\n')


def get_erz_list(input_dir):
    erz_list = set()
    for file in glob.glob(os.path.join(input_dir, 'ERZ*')):
        acc = file.strip().split('/')[-1].split('.')[0]
        erz_list.add(acc)
    return erz_list


def parse_args():
    parser = argparse.ArgumentParser(description='The script creates a file that matches ERZ and read accessions')
    parser.add_argument('-d', '--input_dir', required=True,
                        help='Path to the directory containing downloaded assemblies')
    parser.add_argument('-o', '--outfile', required=True,
                        help='Path to outfile')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input_dir, args.outfile)