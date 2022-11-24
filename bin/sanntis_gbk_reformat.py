#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, metavar='FILE')
    parser.add_argument("--contigs")
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    #outbed=".".join([args.name,"bed"])
    #print(args.name)

    with open(args.gff) as f, open(args.contigs, "r") as fn, open(args.name+"_sn.bed") as b, open(args.name+"_sn.fna") as fnaout :
        record_dict=SeqIO.index(fn, "fasta")
        i=1
        next(f)
        for line in f:
            line=line.rstrip()
            contig_name=line.split("\t")[0]
            start=line.split("\t")[3]
            end=line.split("\t")[4]
            feature=args.name+"_sanntis_"+str(i)
            i+=1
            b.write("%s\t%s\t%\t%s\n" %(contig_name,start,end,feature))
            #extract sequence from prediction
            record=record_dict[contig_name]
            record.id=feature
            record.seq=record.seq[int(start-1):int(end)]
            SeqIO.write(record,fnaout,"fasta")

if __name__ == "__main__":
    main()
