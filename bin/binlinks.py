#!/usr/bin/env python3

import pysam
from Bio import SeqIO
from collections import defaultdict
import os
import argparse
import logging
import csv

BIN_DEFAULT_SEPARATOR = '_'
BAM_DEFAULT_SEPARATOR = '_'

def is_in(read, contig_map, within=1000):
    if read.reference_name not in contig_map.keys():
        return False
    if read.reference_start <= within or read.reference_end <= within:
        return True
    elif read.reference_start > (
        contig_map[read.reference_name] - within
    ) or read.reference_end > (contig_map[read.reference_name] - within):
        return True
    else:
        return False


def keep_read(read, contig_map, within=1000, min_ANI=98, min_cov=0):
    ani = (
        (read.query_alignment_length - read.get_tag("NM"))
        / float(read.query_alignment_length)
        * 100
    )
    cov = read.query_alignment_length / float(read.query_length) * 100
    if ani >= min_ANI and cov >= min_cov and is_in(read, contig_map, within) is True:
        return True
    else:
        return False


def contig_map(bindir, bin_sep, suffix=".fa"):
    logging.debug("Map contigs")
    m = {}
    bins = [item for item in os.listdir(bindir) if 'unbinned' not in item]
    for f in bins:
        logging.debug(f"Processing {f}")
        if f.endswith(suffix) is False:
            continue
        path = os.path.join(bindir, f)
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if bin_sep != BIN_DEFAULT_SEPARATOR:
                    m[record.name.replace(bin_sep, BIN_DEFAULT_SEPARATOR)] = len(record.seq)
                else:
                    m[record.name] = len(record.seq)
    logging.debug(f"Map contigs {len(m)}")
    return m


def bin_map(bindir, bin_sep, suffix=".fa"):
    logging.debug("Bin map")
    contigs = defaultdict(str)
    contigs_per_bin = defaultdict(int)
    bins = [item for item in os.listdir(bindir) if 'unbinned' not in item]
    for f in bins:
        if f.endswith(suffix) is False:
            continue
        path = os.path.join(bindir, f)
        binname = os.path.basename(f)
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contigs[record.name] = binname
                contigs_per_bin[binname] += 1
    logging.debug(f"Bin map contigs:{len(contigs)}, contigs_per_bin:{len(contigs_per_bin)}")
    logging.debug(f'Change separator in bin dictionary if {bin_sep} != {BIN_DEFAULT_SEPARATOR}')
    return_contigs = defaultdict(str)
    if bin_sep != BIN_DEFAULT_SEPARATOR:
        for i in contigs:
            return_contigs[i.replace(bin_sep, BIN_DEFAULT_SEPARATOR)] = contigs[i]
        return return_contigs, contigs_per_bin
    else:
        return contigs, contigs_per_bin


def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    From: https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch():
        if not read.is_paired or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def main():
    # set arguments
    # arguments are passed to classes
    parser = argparse.ArgumentParser(
        description="Evaluate completeness and contamination of a MAG."
    )
    parser.add_argument("--bindir", dest="bindir", type=str, help="Run script on these bins")
    parser.add_argument("--bam",
        dest="bam",
        type=str,
        help="Bam with allr eads aligned against all contigs making up the bins",
    )
    parser.add_argument("--bam-separator",
        dest="bam_sep",
        type=str,
        help="Separator in BAM names, ex ERZ.1 (sep = .) or ERZ_1 (sep = _)",
        default=BAM_DEFAULT_SEPARATOR,
    )
    parser.add_argument("--bin-separator",
        dest="bin_sep",
        type=str,
        help="Separator in BINS.fa names, ex ERZ.1 (sep = .) or ERZ_1 (sep = _)",
        default=BIN_DEFAULT_SEPARATOR,
    )
    parser.add_argument(
        "--out",
        "-o",
        type=str,
        required=False,
        help="Path to output table (Default: links.csv)",
        default="links.csv",
    )
    parser.add_argument(
        "--ANI", type=float, required=False, help="ANI of matching read", default=99
    )
    parser.add_argument(
        "--within",
        type=int,
        required=False,
        help="Within this many bp we need the read to map",
        default=1000,
    )
    parser.add_argument(
        "--contigs",
        "-c",
        action="store_true",
        default=False,
        help="Instead of bins print contigs",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        dest="quiet",
        action="store_true",
        default=False,
        help="Silcence most output",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="Debug and thus ignore safety",
    )
    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S: ",
        level=logLevel,
    )

    samfile = pysam.AlignmentFile(args.bam, "rb")

    cm = contig_map(args.bindir, args.bin_sep)
    bm, contigs_per_bin = bin_map(args.bindir, args.bin_sep)

    link_table = defaultdict(lambda: defaultdict(int))
    bin_table = defaultdict(lambda: defaultdict(int))

    # generate link table
    logging.info("Parsing Bam file. This can take a few moments")
    for read, mate in read_pair_generator(samfile):
        if keep_read(read, cm, args.within, min_ANI=args.ANI) and keep_read(
            mate, cm, args.within, min_ANI=args.ANI
        ):
            # fill in the table
            link_table[read.reference_name][mate.reference_name] += 1
            if read.reference_name != mate.reference_name:
                link_table[mate.reference_name][read.reference_name] += 1
    logging.debug(f"Link table has {len(link_table)} records")

    # generate bin table
    for contig_1, dic in link_table.items():
        for contig_2, links in dic.items():
            bin_table[bm[contig_1]][bm[contig_2]] += 1
    out_data = []
    logging.debug("Constructing output dict")
    if args.contigs:
        logging.debug("contig mode")
        logging.debug(f"Link table has {len(link_table)} records")
        for contig_1, linked in link_table.items():
            for contig_2, links in linked.items():
                out_data.append(
                    {
                        "bin_1": bm[contig_1],
                        "bin_2": bm[contig_2],
                        "contig_1": contig_1,
                        "contig_2": contig_2,
                        "links": links,
                        "bin_1_contigs": contigs_per_bin[bm[contig_1]],
                        "bin_2_contigs": contigs_per_bin[bm[contig_2]],
                    }
                )
    else:
        logging.debug("bins mode")
        logging.debug(f"Bin table has {len(bin_table)} records")
        for bin_1, dic in bin_table.items():
            for bin_2, links in dic.items():
                out_data.append({"bin_1": bin_1, "bin_2": bin_2, "links": links})
    # results
    logging.info("Writing output")
    with open(args.out, "w") as fout:
        cout = csv.DictWriter(fout, fieldnames=list(out_data[0].keys()))
        cout.writeheader()
        for row in out_data:
            cout.writerow(row)


if __name__ == "__main__":
    main()
