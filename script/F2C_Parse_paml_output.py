#!/usr/bin/env python3

import sys
import argparse


######
parser = argparse.ArgumentParser(description='Parse and Stats for CODEML output.')
parser.add_argument('--null_file',help='Null model output file')
parser.add_argument('--alt_file',help='Alternative model output file')
parser.add_argument('--out_file',help='Output file path')
parser.add_argument('--df',help='Degree of freedom',default=1)

args = parser.parse_args()



#######
input_null_model_file_path=args.null_file
input_alternative_model_file_path=args.alt_file
out_file=args.out_file
df=args.df


#######
def extract_likelihood(paml_res_path):
    with open(paml_res_path,'r') as f:
        data = f.read()
    import re
    res = re.findall(
        r'lnL\(.*\):\s*(-\d*\.\d*)\s*\+*',
        data
    )
    lnL=float(res[0])
    return lnL

#######
null_l = extract_likelihood(input_null_model_file_path)
alt_l = extract_likelihood(input_alternative_model_file_path)
d_delt_l = 2*(alt_l-null_l)
import scipy
p=scipy.stats.chi2.pdf(d_delt_l,df=1)


#######
with open(out_file,'w') as f:
    f.write(f'{null_l},{alt_l},{d_delt_l},{p}\n')

print(f'''
Null model log-likelihood: {null_l},
Alternative model log-likelihood: {alt_l},
2_delta_l: {d_delt_l},
p_values: {p}
''')









