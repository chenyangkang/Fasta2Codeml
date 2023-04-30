# Fasta2Codeml
A codeml (PAML package) to make life easier. Dummy input unaligned multi-species fasta file (a single gene), and output codeml result.

## Prerequisites
1. Codeml (PAML version 4.10.6)
2. MACSE (.jar form)
3. MUSCLE
4. RAXML

5. biopython (v1.81, python package)
6. newick (v1.9.0, python package)

must be installed beforehand

## Input data
1. A single gene fasta sequence file (multi-species, not aligned).
2. A text file which indicate the foreground species. On species each line.


## Installation
Simply add ./script to your environment

## Example
cd to example/test_space/

type:
```
Fasta2Codeml.py \
--out_dir /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/test_space \
--project_name Simple_test \
--foreground_file /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/foreground.txt \
--fasta ../single_gene/CLOCK.fasta \
--muscle /beegfs/store4/chenyangkang/software/ParaAT2.0/muscle \
--macse /beegfs/store4/chenyangkang/software/macse_v2.07.jar \
--raxml /beegfs/store4/chenyangkang/software/standard-RAxML/raxml \
--codeml /beegfs/store4/chenyangkang/miniconda3/bin/codeml \
--boostrap 2

```

## To do
Develope multiple alignment mode. For example, porcess multiple cds files and concatenate them, then run codeml on the whole gene.
