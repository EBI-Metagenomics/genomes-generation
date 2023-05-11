#!/usr/bin/env python3

import sys
import os

if len(sys.argv) < 3:
    print("ERROR! usage: python script.py bins_qa.tab bins_taxonomy.tab")
    sys.exit()

bins = {}
with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.strip("\n")
        cols = line.split("\t")
        if cols[0] != "Bin Id":
            name = cols[0]
            complet = float(cols[11])
            cont = float(cols[12])
            heter = float(cols[13])
            bins[name] = [complet, cont, heter]

with open(sys.argv[2], "r") as f:
    for line in f:
        line = line.strip("\n")
        cols = line.split("\t")
        if cols[0] != "Bin Id":
            name = cols[0]
            taxa = cols[3]
            try:
                bins[name].append(taxa)
            except:
                continue

for ele in bins.keys():
    print("%s\t%.2f\t%.2f\t%.2f\t%s" % (ele, bins[ele][0], bins[ele][1], bins[ele][2], bins[ele][3]))