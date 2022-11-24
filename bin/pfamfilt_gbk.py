#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO

def main():
    with open(sys.argv[2], "w") as o:
        for rec in SeqIO.parse(sys.argv[1], "genbank"):
            pf=False
            for feat in (rec.features):
                if feat.type == 'PFAM_domain':
                    pf=True
            if pf:
                SeqIO.write(rec, o, "genbank")

if __name__ == "__main__":
    main()
