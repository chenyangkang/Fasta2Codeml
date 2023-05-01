#!/usr/bin/env python3


#####
import argparse
import os
import subprocess
import time

###### Arguments
parser = argparse.ArgumentParser(description='A codeml (PAML package) wrapper to make life easier. Dummy input unaligned multi-species fasta file (a single gene), and output codeml result.')

parser.add_argument('--out_dir',help='Output directory. Must be full path (!).')
parser.add_argument('--project_name',help='Project name. You can use the gene name as project name. Should contain NO space anywhere.')
parser.add_argument('--foreground_file',help='text format foreground species list, with names corresponding to the headers in the fasta sequence file. Should contain NO space anywhere. Must be full path (!).')
parser.add_argument('--multi_file', action='store_true', default=False, help='Add this tag ff you use multiple fasta file (for example, multiple cds files for a single gene).')
parser.add_argument('--multi_file_list',help='If --multi_file is set, specify a text file with each line as one entry of the bunch of fasta files. CDS files in the path must be full path (!).')
parser.add_argument('--fasta',help='If --Multi_file set as False (default), specify the fasta file path (the only file).')
parser.add_argument('--muscle',help='Path to muscle aligner')
parser.add_argument('--macse',help='Path to macse program (.jar file)')
parser.add_argument('--raxml',help='Path to raxml program')
parser.add_argument('--codeml',help='Path to codeml program')

parser.add_argument('--G',help='Mem use during macse. Use a single number. default is 5G', default=5)
parser.add_argument('--T',help='Thread use while building the tree, default 1', default=1)
parser.add_argument('--boostrap',help='Boostrap while building the tree, Default 100', default=100)
parser.add_argument('--codon_frac',help='remove codons with species miss for more than this fraction. Default 0.5', default=0.5)
parser.add_argument('--sp_frac',help='remove species with sites miss for more than this fraction. Default 0.5', default=0.5)

# parser.add_argument('--macse',help='Path to macse program')



args = parser.parse_args()
print('Multi-file mode: ',args.multi_file)

###
project_path = os.path.join(args.out_dir,args.project_name)

### make project folder
if not os.path.exists(project_path):
    os.mkdir(project_path)


### define linux command wrapper
def process_command(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        buff = p.stdout.readline()
#         err = p.stderr.readlines()
#         if not len(err)==0:
#             print(err)
#             raise
        buff = bytes.decode(buff)
        print(buff)
        if buff == '' and p.poll() != None:
            break

    p.wait()
    
def path_check(path_list):
    for prereq in path_list:
        if not os.path.exists(prereq):
            print(f'{prereq} not exists')
            raise
            
            
## decide multi-or single processing
if args.multi_file:
    ### multi file process
    with open(args.multi_file_list,'r') as f:
        cds_path_list = [i.strip() for i in f.readlines() if not i.strip()=='']
    
    ### chdir
    os.chdir(project_path)
    
    ###
    alignment_list = []
    for cds_path in cds_path_list:
        try:
            ## step 0: get fasta base name
            print(f'============= Start processing {cds_path} =============')
            base_name = cds_path.split('/')[-1]
            if base_name.endswith('.fasta'):
                base_name = '.'.join(base_name.split('.')[:-1])
            else:
                pass
            ### step 1.1: filterout_full_N_seq
            ## run step 1
            print('-------------------- Runing Step 1.1: filterout_full_N_seq -----------------------')
            command = f"F2C_filterout_full_N_seq.py --in_fasta {cds_path} --out_fasta {project_path}/{base_name}.step1.fasta"

            print('Command: ',command)
            process_command(command)
            print('--------------------------- Step 1.1 finished! -----------------------------------\n')

            ### Step 1.2: run muscle alignment
            print('-------------------- Runing Step 1.2: run muscle alignment -----------------------')
            command = f'{args.muscle} -in {project_path}/{base_name}.step1.fasta -out {project_path}/{base_name}.step2.fasta -maxiters 5 -diags'

            print('Command: ',command)
            process_command(command)
            print('--------------------------- Step 1.2 finished! -----------------------------------\n')

            ### Step 1.3: MACSE to refine alignment
            print('-------------------- Runing Step 1.3: run muscle alignment -----------------------')
            command=f'java -jar -Xmx{args.G}G {args.macse} -prog refineAlignment -align {project_path}/{base_name}.step2.fasta -out_AA {project_path}/{base_name}.step3.AA.fasta -out_NT {project_path}/{base_name}.step3.NT.fasta'

            print('Command: ',command)
            process_command(command)
            print('--------------------------- Step 1.3 finished! -----------------------------------\n')

            ### Step 1.4: Trim alignment using MACSE
            print('-------------------- Runing Step 4: Trim alignment using MACSE -----------------------')
            command=f'java -jar -Xmx{args.G}G {args.macse} -prog exportAlignment -codonForInternalFS NNN -codonForInternalStop NNN -charForRemainingFS N -codonForExternalFS NNN -codonForFinalStop NNN -align {project_path}/{base_name}.step3.NT.fasta -out_NT {project_path}/{base_name}.step4.fasta'

            print('Command: ',command)
            process_command(command)
            print('--------------------------- Step 1.4 finished! -----------------------------------\n\n\n')
            
            alignment_list.append(f'{project_path}/{base_name}.step4.fasta')
            
        except Exception as e:
            print(f'While processing {cds_path}, error occured: ')
            print(e)
            print(f'Continue processing other files...')
            
    ### Step 2.1: Concatenate files
    with open(f'{project_path}/alignment_file_list.txt','w') as f:
        for line in alignment_list:
            f.write(line+'\n')
    time.sleep(3)
            
    print('-------------------- Runing Step 2.1: Concatenate files -----------------------')
    command=f'F2C_concatenate_alignment.py --alignment_list_file {project_path}/alignment_file_list.txt --out_path {project_path}/{args.project_name}.step5.fasta'

    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2.1 finished! -----------------------------------\n')

    ### Step 2.1.5 Trim alignment by codons missing fraction and species missing fraction
    print('-------------------- Step 2.1.5 Trim alignment by codons missing fraction -----------------------')
    command=f'F2C_filter_codons_by_fraction.py --input_path {project_path}/{args.project_name}.step5.fasta --output_path {project_path}/{args.project_name}.step5.5.fasta --codon_frac {args.codon_frac} --sp_frac {args.sp_frac}'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2.1.5 finished! -----------------------------------\n')
    
    
    ### Step 2.2: Build tree
    print('---------------------------- Runing Step 2.2: Build tree --------------------------')
    path_check([f'{project_path}/{args.project_name}.step5.5.fasta'])
    command=f'F2C_convert.sh {project_path}/{args.project_name}.step5.5.fasta > {project_path}/{args.project_name}.step6.phy'
    print('Command: ',command)
    process_command(command)
    
    
    command=f'{args.raxml} -f a -x 42 -p 42 -# {args.boostrap} -m GTRGAMMA -s {project_path}/{args.project_name}.step6.phy -n {args.project_name} -T {args.T} -w {project_path}'
    print('Command: ',command)
    process_command(command)
    
    print('--------------------------- Step 2.2 finished! -----------------------------------\n')
    
    
    ### Step 2.3: Co-filter alignment and fasta
    print('-------------------- Runing Step 2.3: Co-filter alignment and fasta -----------------------')
    path_check([args.foreground_file,f'{project_path}/{args.project_name}.step5.5.fasta',f'{project_path}/RAxML_bestTree.{args.project_name}'])
    
    command=f'F2C_co_filter_aln_and_tree.py --foreground_file {args.foreground_file} --in_fasta {project_path}/{args.project_name}.step5.5.fasta --in_tree {project_path}/RAxML_bestTree.{args.project_name} --out_phy {project_path}/{args.project_name}.step7.phy --out_tree {project_path}/{args.project_name}.step7.tre'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2.3 finished! -----------------------------------\n')
    
    ### Step 2.4: Make codeml configure files
    print('-------------------- Runing Step 2.4: Make codeml configure files -----------------------')
    path_check([f'{args.project_name}.step7.phy', f'{args.project_name}.step7.tre'])
    
    command=f'F2C_make_paml_ctl_file.py --input_seq_path {args.project_name}.step7.phy --input_tree_path {args.project_name}.step7.tre --ctl_path_null {args.project_name}.null.ctl --ctl_path_alt {args.project_name}.alternative.ctl --paml_output_null {args.project_name}.null.res.txt --paml_output_alt {args.project_name}.alt.res.txt'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2.4 finished! -----------------------------------\n')
    
    ### Step 2.5: Run codeml
    print('-------------------- Runing Step 2.5: Run codeml -----------------------')
    print('----------------- Runing Alternative Model --------------------')
    command=f'codeml {args.project_name}.alternative.ctl > {args.project_name}.alt.log'
    print('Command: ',command)
    process_command(command)
    print('----------------- Alternative Model Done --------------------')
    print('----------------- Runing Null Model --------------------')
    command=f'codeml {args.project_name}.null.ctl > {args.project_name}.null.log'
    print('Command: ',command)
    process_command(command)
    print('----------------- Null Model Done --------------------')
    print('--------------------------- Step 2.5 finished! -----------------------------------\n')
    
    ### Step 2.6: Stats codeml results
    print('-------------------- Runing Step 2.6: Stats codeml results -----------------------')
    command=f'F2C_Parse_paml_output.py --null_file {project_path}/{args.project_name}.null.res.txt --alt_file {project_path}/{args.project_name}.alt.res.txt --out_file {project_path}/final_stats.txt --df 1'
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 2.6 finished! -----------------------------------\n')
    
    print('\nFINISH.')
    

    pass
    
else:
    ### single file process
    ### chdir
    os.chdir(project_path)

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
    
    
    ### Step 4.5 Trim alignment by codons missing fraction and species missing fraction
    print('-------------------- Step 4.5 Trim alignment by codons missing fraction -----------------------')
    command=f'F2C_filter_codons_by_fraction.py --input_path {project_path}/{base_name}.step4.fasta --output_path {project_path}/{base_name}.step4.5.fasta --codon_frac {args.codon_frac} --sp_frac {args.sp_frac}'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 4.5 finished! -----------------------------------\n')
    
    
    ### Step 5: Build tree
    print('---------------------------- Runing Step 5: Build tree --------------------------')
    command=f'F2C_convert.sh {project_path}/{base_name}.step4.5.fasta > {project_path}/{base_name}.step5.phy'
    print('Command: ',command)
    process_command(command)
    
    
    command=f'{args.raxml} -f a -x 42 -p 42 -# {args.boostrap} -m GTRGAMMA -s {project_path}/{base_name}.step5.phy -n {base_name} -T {args.T} -w {project_path}'
    print('Command: ',command)
    process_command(command)
    
    print('--------------------------- Step 5 finished! -----------------------------------\n')
    
    ### Step 6: Co-filter alignment and fasta
    print('-------------------- Runing Step 6: Co-filter alignment and fasta -----------------------')
    path_check([args.foreground_file, f'{project_path}/{base_name}.step4.5.fasta', f'{project_path}/RAxML_bestTree.{base_name}'])
    command=f'F2C_co_filter_aln_and_tree.py --foreground_file {args.foreground_file} --in_fasta {project_path}/{base_name}.step4.5.fasta --in_tree {project_path}/RAxML_bestTree.{base_name} --out_phy {project_path}/{base_name}.step6.phy --out_tree {project_path}/{base_name}.step6.tre'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 6 finished! -----------------------------------\n')
    
    ### Step 7: Make codeml configure files
    print('-------------------- Runing Step 7: Make codeml configure files -----------------------')
    path_check([f'{base_name}.step6.phy', f'{base_name}.step6.tre'])
    command=f'F2C_make_paml_ctl_file.py --input_seq_path {base_name}.step6.phy --input_tree_path {base_name}.step6.tre --ctl_path_null {base_name}.null.ctl --ctl_path_alt {base_name}.alternative.ctl --paml_output_null {base_name}.null.res.txt --paml_output_alt {base_name}.alt.res.txt'
    
    print('Command: ',command)
    process_command(command)
    print('--------------------------- Step 7 finished! -----------------------------------\n')
    
    ### Step 8: Run codeml
    print('-------------------- Runing Step 8: Run codeml -----------------------')
    print('----------------- Runing Alternative Model --------------------')
    command=f'codeml {base_name}.alternative.ctl > {base_name}.alt.log'
    print('Command: ',command)
    process_command(command)
    print('----------------- Alternative Model Done --------------------')
    print('----------------- Runing Null Model --------------------')
    command=f'codeml {base_name}.null.ctl > {base_name}.null.log'
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
    
    
    
    
    
    
    
    