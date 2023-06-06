#!/usr/bin/env python3

import argparse
import sys
import random
import re
bases = {"R": ["A", "G"],
         "Y": ["C", "T"],
         "K": ["G", "T"],
         "M": ["A", "C"],
         "S": ["C", "G"],
         "W": ["A", "T"],
         "B": ["C", "G", "T"],
         "D": ["A", "G", "T"],
         "H": ["A", "C", "T"],
         "V": ["A", "C", "G"],
         "N": ["A", "G", "C", "T"]
         }

dna = re.compile("[^AGCT]{1}")

vbases = ["A","C","G","T"]

def ren_fasta(args):
    n = 0
    for line in open(args.fasta_file, "r"):
        if line[0] == ">":
            name = line.strip("\n").replace(">","")
            n += 1
            if args.clean and args.prefix is not None:
                print(">%s_%i" % (args.prefix, n))
            elif args.clean and args.prefix is None:
                print(">%s" % (name.split()[0]))
            else:
                print(">%s_%i\t%s" % (args.prefix, n, name))
        else:
            print(line.strip("\n"))

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Rename multifasta file')
        parser.add_argument('-f', dest='fasta_file', help='Input FASTA file')
        parser.add_argument('-p', dest='prefix', help='Header prefix', default=None)
        parser.add_argument('--clean', help="Retain only new header (default: False)", default=False, action="store_true")
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit()
        else:
            args = parser.parse_args()
            ren_fasta(args) 

