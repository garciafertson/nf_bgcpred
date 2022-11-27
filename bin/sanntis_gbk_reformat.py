#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, metavar='FILE')
    parser.add_argument("--contigs", required=True, metavar='FILE')
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    #outbed=".".join([args.name,"bed"])
    #print(args.name)

    with open(args.gff, "r") as f, open(args.name+"_sn.bed", "w") as b, open(args.name+"_sn.fna", "w") as fnaout :
        record_dict=SeqIO.index(args.contigs, "fasta")
        i=1
        next(f)
        for line in f:
            line=line.rstrip()
            contig_name=line.split("\t")[0]
            start=line.split("\t")[3]
            end=line.split("\t")[4]
            feature="sn|"+args.name+"_clstr_"+str(i)
            i+=1
            b.write("%s\t%s\t%s\t%s\n" %(contig_name,start,end,feature))
            #extract sequence from prediction
            record=record_dict[contig_name]
            record.id=feature
            record.seq=record.seq[(int(start)-1):int(end)]
            SeqIO.write(record,fnaout,"fasta")

if __name__ == "__main__":
    main()
