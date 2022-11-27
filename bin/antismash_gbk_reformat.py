#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--gbk", required=True, metavar='FILE')
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    #outbed=".".join([args.name,"bed"])
    #print(args.name)
    with open(args.gbk, "r") as f, open(args.name+"_as.bed", "w") as b, open(args.name+"_as.fna", "w") as fn:
        i=1
        for record in SeqIO.parse(f,"genbank"):
            contig=record.description
            start=record.annotations["structured_comment"]['antiSMASH-Data']['Orig. start']
            end=record.annotations["structured_comment"]['antiSMASH-Data']['Orig. end']
            record.id="as|" + args.name+"-"+record.id+"_"+str(i)
            record.description=""
            i+=1
            SeqIO.write(record, fn, "fasta")
            b.write("%s\t%s\t%s\t%s\n" %(contig, start, end, record.id))

if __name__ == "__main__":
    main()
