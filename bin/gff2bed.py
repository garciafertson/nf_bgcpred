#!/usr/bin/env python3
import argparse

def main( ):
    parser=argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, metavar='FILE')
    parser.add_argument("--name", required=True, metavar='FILE')
    args=parser.parse_args()

    #outbed=".".join([args.name,"bed"])
    #print(args.name)
    with open(args.input_gbk, "r") as f, open(args.name+".bed", "w") as o:
        for line in f:
            line=line.rstrip()
            if line.startswith("#"):
                next()
            else:
                bgcid=line.split("\t")[0]
                start=line.split("\t")[2]
                end=line.split("\t")[3]
                strand=line.split("\t")[5]
                geneid="_".join([bgcid,start,end,strand])
                o.write("%s\t%s\t%s\t%s\n" %(bgcid, start, end, geneid))

if __name__ == "__main__":
    main()
