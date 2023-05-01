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
    parser = argparse.ArgumentParser(description='Filter codons/species by fractions missing.')
    parser.add_argument('--input_path',help='Input fasta file')
    parser.add_argument('--output_path',help='Output path (path+name)')
    parser.add_argument('--codon_frac',help='remove codons with species miss for more than this fraction', default=0.5)
    parser.add_argument('--sp_frac',help='remove species with sites miss for more than this fraction', default=0.5)

    args = parser.parse_args()
    
    ####
    out = AlignIO.read(args.input_path,'fasta')
    new_out = AlignIO.read(args.input_path,'fasta')
    
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
        if nonsense/np.sum(list(counted.values()))>float(args.codon_frac):
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
        if nonsense/np.sum(list(counted.values()))>float(args.sp_frac):
            continue
        else:
            new_new_out+=f'>{sp.name}\n'
            new_new_out+=str(sp.seq)
            new_new_out+='\n'


    with open(args.output_path,'w') as f:
        f.write(new_new_out)
    

    