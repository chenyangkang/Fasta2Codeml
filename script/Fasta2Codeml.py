#!/usr/bin/env python3


#####
import argparse
import os
import subprocess

###### Arguments
parser = argparse.ArgumentParser(description='A codeml (PAML package) wrapper to make life easier. Dummy input unaligned multi-species fasta file (a single gene), and output codeml result.')

parser.add_argument('--out_dir',help='Output directory. Must be full path (!).')
parser.add_argument('--project_name',help='Project name. You can use the gene name as project name. Should contain NO space anywhere.')
parser.add_argument('--foreground_file',help='text format foreground species list, with names corresponding to the headers in the fasta sequence file. Should contain NO space anywhere. Must be full path (!).')
parser.add_argument('--Multi_file', action='store_true', default=False, help='Add this tag ff you use multiple fasta file (for example, multiple cds files for a single gene)')
parser.add_argument('--Multi_file_list',help='If --Multi_file is set, specify a text file with each line as one entry of the bunch of fasta files.')
parser.add_argument('--fasta',help='If --Multi_file set as False (default), specify the fasta file path (the only file).')
parser.add_argument('--muscle',help='Path to muscle aligner')
parser.add_argument('--macse',help='Path to macse program (.jar file)')
parser.add_argument('--raxml',help='Path to raxml program')
parser.add_argument('--codeml',help='Path to codeml program')

parser.add_argument('--G',help='Mem use during macse. Use a single number. default is 5G', default=5)
parser.add_argument('--T',help='Thread use while building the tree, default 1', default=1)
parser.add_argument('--boostrap',help='Boostrap while building the tree', default=100)

# parser.add_argument('--macse',help='Path to macse program')



args = parser.parse_args()
print(args.Multi_file)

###
project_path = os.path.join(args.out_dir,args.project_name)

### make project folder
if not os.path.exists(project_path):
    os.mkdir(project_path)

### chdir
os.chdir(project_path)

### define linux command wrapper
def process_command(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        buff = p.stdout.readline()
        buff = bytes.decode(buff)
        print(buff)
        if buff == '' and p.poll() != None:
            break

    p.wait()
    

## decide multi-or single processing
if args.Multi_file:
    ### multi file process
    pass
    
else:
    ### single file process
    ## step 0: get fasta base name
    base_name = args.fasta.split('/')[-1]
    if base_name.endswith('.fasta'):
        base_name = '.'.join(base_name.split('.')[:-1])
    else:
        pass
        
        
    ### step 1: filterout_full_N_seq
    ## run step 1
    print('-------------------- Runing Step 1: filterout_full_N_seq -----------------------')
    command = f"F2C_filterout_full_N_seq.py --in_fasta {args.fasta} --out_fasta {project_path}/{base_name}.step1.fasta"
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 1 finished! -----------------------------------\n')
    
    ### Step 2: run muscle alignment
    print('-------------------- Runing Step 2: run muscle alignment -----------------------')
    command = f'{args.muscle} -in {project_path}/{base_name}.step1.fasta -out {project_path}/{base_name}.step2.fasta -maxiters 5 -diags'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2 finished! -----------------------------------\n')
    
    ### Step 3: MACSE to refine alignment
    print('-------------------- Runing Step 3: run muscle alignment -----------------------')
    command=f'java -jar -Xmx{args.G}G {args.macse} -prog refineAlignment -align {project_path}/{base_name}.step2.fasta -out_AA {project_path}/{base_name}.step3.AA.fasta -out_NT {project_path}/{base_name}.step3.NT.fasta'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 3 finished! -----------------------------------\n')
    
    ### Step 4: Trim alignment using MACSE
    print('-------------------- Runing Step 4: Trim alignment using MACSE -----------------------')
    command=f'java -jar -Xmx{args.G}G {args.macse} -prog exportAlignment -codonForInternalFS NNN -codonForInternalStop NNN -charForRemainingFS N -codonForExternalFS NNN -codonForFinalStop NNN -align {project_path}/{base_name}.step3.NT.fasta -out_NT {project_path}/{base_name}.step4.fasta'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 4 finished! -----------------------------------\n')
    
    ### Step 5: Build tree
    print('---------------------------- Runing Step 5: Build tree --------------------------')
    command=f'F2C_convert.sh {project_path}/{base_name}.step2.fasta > {project_path}/{base_name}.step5.phy'
    print('Command: ',command)
    process_command(command)
    
    
    command=f'{args.raxml} -f a -x 42 -p 42 -# {args.boostrap} -m GTRGAMMA -s {project_path}/{base_name}.step5.phy -n {base_name} -T {args.T} -w {project_path}'
    print('Command: ',command)
    process_command(command)
    
    print('--------------------------- Step 5 finished! -----------------------------------\n')
    
    ### Step 6: Co-filter alignment and fasta
    print('-------------------- Runing Step 6: Co-filter alignment and fasta -----------------------')
    
    command=f'F2C_co_filter_aln_and_tree.py --foreground_file {args.foreground_file} --in_fasta {project_path}/{base_name}.step4.fasta --in_tree {project_path}/RAxML_bestTree.{base_name} --out_phy {project_path}/{base_name}.step6.phy --out_tree {project_path}/{base_name}.step6.tre'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 6 finished! -----------------------------------\n')
    
    ### Step 7: Make codeml configure files
    print('-------------------- Runing Step 7: Make codeml configure files -----------------------')
    command=f'F2C_make_paml_ctl_file.py --input_seq_path {project_path}/{base_name}.step6.phy --input_tree_path {project_path}/{base_name}.step6.tre --ctl_path_null {project_path}/{base_name}.null.ctl --ctl_path_alt {project_path}/{base_name}.alternative.ctl --paml_output_null {project_path}/{base_name}.null.res.txt --paml_output_alt {project_path}/{base_name}.alt.res.txt'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 7 finished! -----------------------------------\n')
    
    ### Step 8: Run codeml
    print('-------------------- Runing Step 8: Run codeml -----------------------')
    print('----------------- Runing Alternative Model --------------------')
    command=f'codeml {project_path}/{base_name}.alternative.ctl > {project_path}/{base_name}.alt.log'
    print('Command: ',command)
    process_command(command)
    print('----------------- Alternative Model Done --------------------')
    print('----------------- Runing Null Model --------------------')
    command=f'codeml {project_path}/{base_name}.null.ctl > {project_path}/{base_name}.null.log'
    print('Command: ',command)
    process_command(command)
    print('----------------- Null Model Done --------------------')
    print('--------------------------- Step 8 finished! -----------------------------------\n')
    
    ### Step 9: Stats codeml results
    print('-------------------- Runing Step 9: Stats codeml results -----------------------')
    command=f'F2C_Parse_paml_output.py --null_file {project_path}/{base_name}.null.res.txt --alt_file {project_path}/{base_name}.alt.res.txt --out_file {project_path}/final_stats.txt --df 1'
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 9 finished! -----------------------------------\n')
    
    print('\nFINISH.')
    
    
    
    
    
    
    
    