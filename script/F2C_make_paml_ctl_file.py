#!/usr/bin/env python3


import sys
import argparse

if __name__=='__main__':
    ######
    parser = argparse.ArgumentParser(description='Generate codeml null and alternative model input file (branch-site).')
    parser.add_argument('--input_seq_path',help='Input phylip alignment sequence file')
    parser.add_argument('--input_tree_path',help='Input tree file')
    parser.add_argument('--ctl_path_null',help='Output name for null model configuration')
    parser.add_argument('--ctl_path_alt',help='Output name for alternative model configuration')

    parser.add_argument('--paml_output_null',help='Output path for null model running result')
    parser.add_argument('--paml_output_alt',help='Output path for alternative model running result')
    

    args = parser.parse_args()


    input_seq_path = args.input_seq_path
    input_tree_path = args.input_tree_path
    ctl_path_null = args.ctl_path_null
    ctl_path_alt = args.ctl_path_alt

    paml_output_null = args.paml_output_null
    paml_output_alt = args.paml_output_alt



    null_model=f'''
          seqfile = {input_seq_path}            * Path to the alignment file
         treefile = {input_tree_path}           * Path to the tree file
          outfile = {paml_output_null}            * Path to the output file

            noisy = 3              * How much rubbish on the screen
          verbose = 1              * More or less detailed report

          seqtype = 1              * Data type
            ndata = 1           * Number of data sets or loci
            icode = 0              * Genetic code 
        cleandata = 0              * Remove sites with ambiguity data?

            model = 2         * Models for ω varying across lineages
          NSsites = 2          * Models for ω varying across sites
        CodonFreq = 7        * Codon frequencies
          estFreq = 0        * Use observed freqs or estimate freqs by ML
            clock = 0          * Clock model
        fix_omega = 1         * Estimate or fix omega
            omega = 1        * Initial or fixed omega

    '''

    with open(ctl_path_null,'w') as f:
        f.write(null_model)



    alternative_model=f'''
          seqfile = {input_seq_path}            * Path to the alignment file
         treefile = {input_tree_path}           * Path to the tree file
          outfile = {paml_output_alt}            * Path to the output file

            noisy = 3              * How much rubbish on the screen
          verbose = 1              * More or less detailed report

          seqtype = 1              * Data type
            ndata = 1           * Number of data sets or loci
            icode = 0              * Genetic code 
        cleandata = 0              * Remove sites with ambiguity data?

            model = 2         * Models for ω varying across lineages
          NSsites = 2          * Models for ω varying across sites
        CodonFreq = 7        * Codon frequencies
          estFreq = 0        * Use observed freqs or estimate freqs by ML
            clock = 0          * Clock model
        fix_omega = 0         * Estimate or fix omega
            omega = 0.5        * Initial or fixed omega

    '''


    with open(ctl_path_alt,'w') as f:
        f.write(alternative_model)




