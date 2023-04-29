#!/usr/bin/env python3

from Bio import AlignIO
from Bio import SeqIO
import sys
import argparse

######
if __name__=='__main__':
    ##### parse argument
    parser = argparse.ArgumentParser(description='Filter out species in the fasta if all base are N (missing).')
    parser.add_argument('--in_fasta',help='input fasta')
    parser.add_argument('--out_fasta',help='output fasta')

    args = parser.parse_args()

    ##### process
    data=SeqIO.parse(args.in_fasta,'fasta')
    new_fasta=[]
    for rec in data:
        if not len([i for i in rec.seq if not i=='N'])==0:
            new_fasta.append(rec)
    SeqIO.write(new_fasta,args.out_fasta,'fasta')


