#!/usr/bin/env python3
#
# This file is part of the EukCC (https://github.com/openpaul/eukcc).
# Copyright (c) 2019 Paul Saary
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# provides all file operation functions
# used inthis package
import os
import argparse
import subprocess
import logging
import tempfile
import gzip
import glob
import csv


# backup fasta handler, so we can use readonly directories
class fa_class:
    def __init__(self, seq, name, long_name):
        self.seq = seq
        self.name = name
        self.long_name = long_name

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)


def Fasta(path):
    """
    Iterator for fasta files
    """
    entry = False

    with open(path) as fin:
        for line in fin:
            if line.startswith(">"):
                if entry is not False:
                    entry.seq = "".join(entry.seq)
                    yield entry
                # define new entry
                long_name = line.strip()[1:]
                name = long_name.split()[0]
                entry = fa_class([], name, long_name)
            else:
                entry.seq.append(line.strip())
        # yield last one
        entry.seq = "".join(entry.seq)
        yield entry


def N50(l, x=0.5):
    l = sorted(l)
    l.reverse()
    t = sum(l) * x
    for i, x in enumerate(l):
        if i == 0:
            cs = x
        else:
            cs = cs + x
        if cs >= t:
            return x


def genome_stats(path):
    lengths = [len(x) for x in Fasta(path)]
    return {"N50": N50(lengths), "bp": sum(lengths), "contigs": len(lengths)}


def eukcc_score(path):
    with open(path) as fin:
        for row in csv.DictReader(fin, delimiter="\t"):
            return row


def busco_score(path):
    d = {
        "BUSCO_C": None,
        "BUSCO_M": None,
        "BUSCO_D": None,
        "BUSCO_S": None,
        "BUSCO_F": None,
    }

    with open(path) as fin:
        for line in fin:
            if line.strip().startswith("C:"):
                elem = line.strip().replace("[", ",").replace("]", "").split(",")
                for e in elem:
                    x = e.split(":")
                    d["BUSCO_{}".format(x[0])] = x[1].replace("%", "")
                # this is the last info we need
                break
            if line.strip().startswith("# The lineage dataset is: "):
                d["BUSCO_lineage"] = line.strip().split()[5]
    return d


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", help="path for the output table", default="qc.csv", type=str
    )
    parser.add_argument("--mags", help="path to mag folder", type=str)
    parser.add_argument("--qc_dir", help="path to qc folder", type=str)
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="rerun even if output exists",
        default=False,
    )
    parser.add_argument(
        "--quiet", action="store_true", help="supress information", default=False
    )
    parser.add_argument(
        "--debug", action="store_true", help="Make it more verbose", default=False
    )

    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S: ", level=logLevel
    )

    # find mags
    mags = glob.glob("{}/*.fa".format(args.mags))

    # for each mag learn eukcc and BUSCO values
    # also determine N50
    #
    fields = [
        "fasta",
        "completeness",
        "contamination",
        "BUSCO_C",
        "BUSCO_M",
        "BUSCO_D",
        "BUSCO_S",
        "BUSCO_F",
        "BUSCO_n",
        "BUSCO_lineage",
        "N50",
        "bp",
        "contigs",
    ]

    with open(args.output, "w") as outfile:
        cout = csv.DictWriter(outfile, fieldnames=fields, extrasaction="ignore")
        cout.writeheader()

        for mag in mags:
            name = os.path.basename(mag).replace(".fa", "")
            # learn eukcc score
            eukcc_p = os.path.join(args.qc_dir, "eukcc", name, "eukcc.csv")
            busco_p = os.path.join(
                args.qc_dir, "busco", name, "short_summary.specific.txt"
            )
            stat = genome_stats(mag)
            eukcc = eukcc_score(eukcc_p)
            stat = {**stat, **eukcc}
            stat["fasta"] = os.path.basename(mag)
            busco = busco_score(busco_p)
            stat = {**stat, **busco}
            cout.writerow(stat)
