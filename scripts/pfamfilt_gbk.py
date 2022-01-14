from Bio import SeqIO
import sys

#for rec in SeqIO.parse(sys.argv[1], "genbank"):
#   SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")

with open(sys.argv[2], "w") as o:
    for rec in SeqIO.parse(sys.argv[1], "genbank"):
        pf=False
        for feat in (rec.features):
            if feat.type == 'PFAM_domain':
                pf=True
        if pf:
            SeqIO.write(rec, o, "genbank")
