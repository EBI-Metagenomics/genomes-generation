#!/usr/bin/env python3

import sys
import argparse

def get_contamination(cat_file, tax):
    total_bp = 0
    total_cont = 0
    cont_tax = []
    main_tax = "Not classified"
    linen = 0
    with open(cat_file, "r") as f:
        for line in f:
            linen += 1
            if linen == 1:
                total_bp = int(line.split()[-2].replace(',', ''))
            else:
                line = line.rstrip()
                cols = line.split("\t")
                if cols[0] == tax:
                    length = int(cols[-1])
                    taxon = cols[1]
                    if taxon != "not classified" and taxon != "NA" and main_tax == "Not classified":
                        main_tax = taxon
                    elif taxon != "not classified" and taxon != "NA":
                        cont_tax.append(taxon)
                        total_cont += length
    cont = float(total_cont)/total_bp*100
    return main_tax, cont_tax, cont


def get_contigs(contig_stats, taxa):
    excl_contigs = set()
    with open(contig_stats, "r") as f:
        for line in f:
            if line[0] != "#":
                line = line.rstrip()
                cols = line.split("\t")
                if cols[1] == "classified":
                    contig = cols[0]
                    kingdom = cols[6].split(":")[0]
                    phylum = cols[7].split(":")[0]
                    if kingdom in taxa or phylum in taxa:
                        excl_contigs.add(contig)
    return excl_contigs

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Detect decontamination on CAT results"
    )
    parser.add_argument(
        "-s", "--summary", dest="summary", help="[genome_prefix].summary.txt", required=True
    )
    parser.add_argument(
        "-n", "--names", dest="names", help="[genome_prefix].contig2classification.official_names.txt", required=False
    )
    parser.add_argument(
        "-i", "--sample-name", dest="sample_name", help="genome prefix", required=False
    )
    args = parser.parse_args()

    res_kingdom = get_contamination(args.summary, "superkingdom")
    res_phylum = get_contamination(args.summary, "phylum")
    with open(args.sample_name+".cont-stats.tsv", "w") as outstats:
        outstats.write("Lineage\t%s\t%s\n" % (res_kingdom[0], res_phylum[0]))
        outstats.write("Contamination\t%.2f%% (%s)\t%.2f%% (%s)\n" % (res_kingdom[2], ",".join(res_kingdom[1]),
                                                                  res_phylum[2], ",".join(res_phylum[1])))
    with open(args.sample_name+".cont-contigs.txt", "w") as outcontigs:
        if res_kingdom[1] or res_phylum[1]:
            excl_tax = res_kingdom[1]+res_phylum[1]
            excl_contigs = get_contigs(args.names, excl_tax)
            for c in excl_contigs:
                outcontigs.write("%s\n" % c)
