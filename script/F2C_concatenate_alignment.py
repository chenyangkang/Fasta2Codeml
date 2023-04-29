#!/usr/bin/env python3

from Bio import AlignIO
from Bio import Seq, SeqIO
import os
import numpy as np
from collections import Counter
import sys
import argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Concatenate multiple alignment based on entry name')
    parser.add_argument('--out_dir',help='Output directory')


    gene_name=sys.argv[1]

    cds_path_list = [i for i in os.listdir(f'codon_NT_alignment') if i.startswith(gene_name)]
    cds_idx_list = sorted([int(i.split('.')[1].split('cds')[1]) for i in cds_path_list])

    from collections import defaultdict
    aln_dict =defaultdict(dict)


    #### append
    for idx in cds_idx_list:
        try:
            path = f'{gene_name}.cds{idx}.aln.NT.Trimed.fasta'
            aln = AlignIO.read(f'codon_NT_alignment/{path}','fasta')
            if not len(aln[0])%3==0:
                continue

            for sp_aln in aln:
                aln_dict[sp_aln.name][idx]=sp_aln
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


    with open(f'./codon_NT_alignment_gene/{gene_name}.fasta','w') as f:
        f.write(out)

    out = AlignIO.read(f'./codon_NT_alignment_gene/{gene_name}.fasta','fasta')
    new_out = AlignIO.read(f'./codon_NT_alignment_gene/{gene_name}.fasta','fasta')

    ##### remove bad columns
    column_id_to_remove=[]
    already_deleated = 0
    for column in range(int(len(out[0])/3)):
        counted = Counter(out[:,column*3]+out[:,column*3+1]+out[:,column*3+2])
        nonsense = 0
        if '-' in counted:
            nonsense+=counted['-']
        if '!' in counted:
            nonsense+=counted['!']
        if 'N' in counted:
            nonsense+=counted['N']
        if nonsense/np.sum(list(counted.values()))>0.7:
            if column==0:
                new_out=new_out[:,column*3+3:]
                already_deleated+=1
            elif column==int(len(out[0])/3)-1:
                new_out=new_out[:,:-3]
                already_deleated+=1
            else:
                new_out=new_out[:,0:(column-already_deleated)*3] + new_out[:,(column+1-already_deleated)*3:]
                already_deleated+=1

            if not len(new_out[0].seq)%3==0:
                raise


    #### remove overlacking species
    new_new_out = ''
    for sp in new_out:
        counted = Counter(str(sp.seq))
        nonsense = 0
        if '-' in counted:
            nonsense+=counted['-']
        if '!' in counted:
            nonsense+=counted['!']
        if 'N' in counted:
            nonsense+=counted['N']
        if nonsense/np.sum(list(counted.values()))>0.7:
            continue
        else:
            new_new_out+=f'>{sp.name}\n'
            new_new_out+=str(sp.seq)
            new_new_out+='\n'


    with open(f'./codon_NT_alignment_gene/{gene_name}.fasta','w') as f:
        f.write(new_new_out)














