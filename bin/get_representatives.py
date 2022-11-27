#!/usr/bin/env python3
import pandas as pd
import re
import sys, os
import argparse
import pandas as pd
from Bio import SeqIO

def file_is_not_empty(path):
    return os.stat(path).st_size!=0

def main():

#read input files, bedfile, cluster file and fasta files
    parser=argparse.ArgumentParser()
    parser.add_argument("--bedfile", required=True, metavar='FILE')
    parser.add_argument("--clusters", required=True, metavar='FILE')
    parser.add_argument("--fna", required=True, metavar='FILE')
    args=parser.parse_args()

    all_bgc=set()
    representative_bgc=set()
    incluster_bgc=set()

    #Readbedfile to get length information
    #with open(args.bedfile, "r") as f:
    bed=pd.read_csv(args.bedfile, index_col=3, header=None, sep='\t')
    bed["len"]=bed[2]-bed[1]
    all_bgc=set(bed.index)

    #Read cluster file, for each line get the corresponding length of bgc from bed file
    #Keep those BGC absent in cluster file
    with open(args.clusters, "r") as f:
        for line in f:
            line=line.rstrip()
            bgcincl=set(line.split("\t"))
            incluster_bgc.update(bgcincl) #incluster_bgc.update(bgcincl)
            rep=bed.loc[bgcincl,"len"].idxmax()
            representative_bgc.add(rep)

    #get set of representatives
    representative_bgc.update(all_bgc.difference(incluster_bgc))
    bglist=list(representative_bgc)
    #output bed file with repesentative bgcs
    bedrep=bed.loc[bglist, [0,1,2]]
    bedrep=bedrep.reset_index()
    bedrep=bedrep.loc[:,[0,1,2,3]]
    bedrep.to_csv("bgc_representative.bed", header=None, index=False, sep="\t")
    with open(args.fna, "r") as f, open("bgc_representative.fna", "w") as o:
        records=SeqIO.parse(f, "fasta")
        for rec in records:
            if rec.id in representative_bgc:
                SeqIO.write(rec, o, "fasta")


    #open fasta file and collect sequnces in representative bgc set
    #df.to_csv("predicted_representative.tsv", index=True,sep="\t")
if __name__ == "__main__":
    main()
