#!/usr/bin/env python3

import sys
from Bio import AlignIO
from Bio import SeqIO
import os
import argparse


######
if __name__=='__main__':

    ######
    parser = argparse.ArgumentParser(description='Co-filter alignment and tree.')
    parser.add_argument('--foreground_file',help='Foreground file')
    parser.add_argument('--in_fasta',help='Input fasta alignment')
    parser.add_argument('--in_tree',help='Input tree')
    parser.add_argument('--out_phy',help='Output phylip alignment')
    parser.add_argument('--out_tree',help='Output tree')
    

    args = parser.parse_args()


    ######

    foreground_list_path = args.foreground_file

    input_alignment_path = args.in_fasta
    output_alignment_path = args.out_phy

    input_species_tree_path=args.in_tree
    ouput_tree_path = args.out_tree


    #### read foreground
    with open(foreground_list_path,'r') as f:
        foreground_list = f.readlines()
        foreground_list = [i.strip() for i in foreground_list if not i.strip()=='']


    #### get alignment entry names
    aln = SeqIO.parse(input_alignment_path,'fasta')
    aln_name_list= [i.name for i in aln]


    #### get tree leaf names
    with open(input_species_tree_path,'r') as f:
        tree = f.read()
    import newick
    tree1 = newick.loads(tree)[0]


    #### get concensus names
    final_leaf_names = set(list(tree1.get_leaf_names())) & set(aln_name_list)



    #### get trimed alignment
    aln = SeqIO.parse(input_alignment_path,'fasta')
    a = [i for i in aln if i.name in final_leaf_names]

    #### get trimed tree
    tree1.prune_by_names(final_leaf_names, inverse=True)
    tree1.remove_redundant_nodes(keep_leaf_name=True) 



    ######
    newick_tree = tree1.newick


    ##### checking point: both foreground and background shoule be in the tree
    if len(set(foreground_list) & set(tree1.get_leaf_names()))==0:
        print('no foreground remained after filtering!')
        raise
    elif len(set(tree1.get_leaf_names())-set(foreground_list))==0:
        print('no background remained after filtering!')
        raise

    #### annotaiton: foreground and background
    for foreground_item in foreground_list:
        newick_tree = newick_tree.replace(foreground_item, foreground_item+f' #1')


    #### checkpoint: tree file taxa equal to alignment taxa
    assert len(a) == len(set(tree1.get_leaf_names()))

    ##### writing tree file
    with open(ouput_tree_path,'w') as f:
        f.write(str(len(set(tree1.get_leaf_names()))) +' 1\n\n')
        f.write(newick_tree.replace(',',', '))
        f.write(';\n\n')

    ###### writing alignment
    with open(output_alignment_path, 'w') as f:
        f.write(str(len(a))+' '+str(len(a[0]))+'\n')
        for seq in a:
            f.write(seq.name + '  ')
            f.write(str(seq.seq)+'\n')



