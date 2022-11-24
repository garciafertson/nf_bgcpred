#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--gbk", required=True, metavar='FILE')
    parser.add_argument("--tsv", required=True, metavar='FILE')
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    #outbed=".".join([args.name,"bed"])
    #print(args.name)
    with open(args.tsv) as f, open(args.gbk, "r") as gb, open(args.name+"_dp.bed") as b, open(args.name+"_dp.fna") as fnaout:
        record_dict=SeqIO.index(gb, "genbank")
        next(f)
        for line in f:
            line=line.rstrip()
            contig_name=line.split("\t")[0]
            start=line.split("\t")[2]
            end=line.split("\t")[3]
            feature=line.split("\t")[1]
            record=record_dict[feature]
            b.write("%s\t%s\t%\t%s\n" %(contig_name,start,end,feature))
            #extract sequence from prediction
            SeqIO.write(record,fnaout,"fasta")

if __name__ == "__main__":
    main()
