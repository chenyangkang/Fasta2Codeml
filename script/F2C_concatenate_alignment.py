#!/usr/bin/env python3

from Bio import AlignIO
from Bio import Seq, SeqIO
import os
import numpy as np
from collections import Counter
import sys
import argparse
from collections import defaultdict

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Concatenate multiple alignment based on entry name')
    parser.add_argument('--alignment_list_file',help='File that specify the alignments. One file name each line.')
    parser.add_argument('--out_path',help='Output path (path+name)')

    args = parser.parse_args()
    

    ####
    aln_dict =defaultdict(dict)

    #### append
    with open(args.alignment_list_file,'r') as f:
        path_list=[i.strip() for i in f.readlines() if not i.strip()=='']
        
    
    for path in path_list:
#         base_name='.'.join(path.split('.fasta'))[:-1]
        try:
            aln = AlignIO.read(path,'fasta')
            if not len(aln[0])%3==0:
                print('Your alignment column count is not multiple of 3. Terminating.')
                raise

            for sp_aln in aln:
                aln_dict[sp_aln.name][path]=sp_aln
        except:
            continue


    #### unique cds idx
    all_idx_list = [list(aln_dict[i].keys()) for i in aln_dict.keys()]
    unique_idx=set()
    for idx in all_idx_list:
        unique_idx.update(idx)
    unique_idx = sorted(list(unique_idx))


    def get_one(aln_dict, cds_idx):
        for sp in aln_dict:
            if cds_idx in aln_dict[sp]:
                return len(aln_dict[sp][cds_idx])

    out=''

    for sp in aln_dict:
        seq=''
        for idx in unique_idx:
            if idx in aln_dict[sp]:
                seq+=str(aln_dict[sp][idx].seq)
            else:
                seq+='-'*get_one(aln_dict, idx)

        out+=f'>{sp}\n'
        out+=seq
        out+='\n'


    with open(args.out_path,'w') as f:
        f.write(out)
















