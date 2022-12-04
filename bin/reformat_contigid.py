#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import argparse
import os
def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--fna", required=True, metavar='FILE')
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    with open(args.name+".fna","w") as fnaout:
        records=SeqIO.parse(args.fna, "fasta")
        for record in records:
            record.id=record.id.split(" ")[0]
            record.id= re.sub('\.[0-9]+$', '', record.id)
            fnaout.write(">%s\n%s\n" %(record.id, record.seq) )

if __name__ == "__main__":
    main()
