#!/usr/bin/env python3
# coding=utf-8

import argparse

STEPS = {
    'GGP:EUK_MAGS_GENERATION:EUKCC_CONCOCT': 0,
     'GGP:EUK_MAGS_GENERATION:EUKCC_METABAT': 1,
     'GGP:EUK_MAGS_GENERATION:FILTER_QUALITY': 2,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:RENAME_AND_CHECK_SIZE_BINS_BINNER1': 3,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:RENAME_AND_CHECK_SIZE_BINS_BINNER2': 4,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:RENAME_AND_CHECK_SIZE_BINS_BINNER3': 5,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:CHECKM2_BINNER1': 6,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:CHECKM2_BINNER2': 7,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:CHECKM2_BINNER3': 8,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE12:BINNING_REFINER': 9,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE23:BINNING_REFINER': 10,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE13:BINNING_REFINER': 11,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE123:BINNING_REFINER': 12,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE12:CHECKM2_REFINE': 13,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE23:CHECKM2_REFINE': 14,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE13:CHECKM2_REFINE': 15,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:REFINE123:CHECKM2_REFINE': 16,
     'GGP:PROK_MAGS_GENERATION:BIN_REFINEMENT:CONSOLIDATE_BINS': 17,
     'GGP:PROK_MAGS_GENERATION:CHECKM2': 18,
     'GGP:PROK_MAGS_GENERATION:DREP': 19
}

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collect stats of processes.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Pipeline logging file")
    parser.add_argument('-o', '--output', type=str, required=False, help="Output file", default="output_logging.txt")
    return parser.parse_args()

def main():
    data = {}
    args = parse_args()
    with open(args.input, 'r') as file_in:
        run_accesison = ''
        for line in file_in:
            line = line.strip()
            if 'GGP' in line:
                run_accesison, step = line.split('\t')[:2]
                data.setdefault(run_accesison, ['' for _ in range(len(STEPS))])
            else:
                data[run_accesison][STEPS[step]] = line

    with open(args.output, 'w') as file_out:
        for run in data:
            file_out.write(run + '\n')
            for id in range(len(data[run])):
                if data[run][id]:
                    step = list(STEPS.keys())[id]
                    file_out.write(step + ' -- ' + data[run][id] + '\n')
    print('Completed')


if __name__ == "__main__":
    main()
