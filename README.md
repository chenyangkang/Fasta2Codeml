# Fasta2Codeml
A codeml (PAML package) wrapper to make life easier. Dummy input unaligned multi-species fasta file (a single gene), and output codeml result.

## Prerequisites
1. Codeml (PAML version 4.10.6)
2. MACSE (.jar form)
3. MUSCLE
4. RAXML

5. biopython (v1.81, python package)
6. newick (v1.9.0, python package)

must be installed beforehand

# Installation
Simply add ./script to your environment


# Input data
1. A single gene fasta sequence file (multi-species, not aligned).
2. A text file which indicate the foreground species. On species each line.


# Example
cd to example/test_space/

**Change the absolute path in the command lines below to to your path.**

##For single fasta file mode
type:
```
Fasta2Codeml.py \
--out_dir /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/test_space_single_gene \
--project_name Simple_test \
--foreground_file /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/foreground.txt \
--fasta /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/single_gene/CLOCK.fasta \
--muscle /beegfs/store4/chenyangkang/software/ParaAT2.0/muscle \
--macse /beegfs/store4/chenyangkang/software/macse_v2.07.jar \
--raxml /beegfs/store4/chenyangkang/software/standard-RAxML/raxml \
--codeml /beegfs/store4/chenyangkang/miniconda3/bin/codeml \
--boostrap 10
```

## For multi-file model (multi-cds file mode):
```
Fasta2Codeml.py \
--out_dir /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/test_space_multi_cds \
--project_name Simple_multi_test \
--foreground_file /beegfs/store4/chenyangkang/DEV/Fasta2Codeml/example/foreground.txt \
--multi_file \
--multi_file_list cds_list.txt \
--muscle /beegfs/store4/chenyangkang/software/ParaAT2.0/muscle \
--macse /beegfs/store4/chenyangkang/software/macse_v2.07.jar \
--raxml /beegfs/store4/chenyangkang/software/standard-RAxML/raxml \
--codeml /beegfs/store4/chenyangkang/miniconda3/bin/codeml \
--boostrap 10
```

# Workflow underneath
1. Remove species that contain only "N"s.
2. Run muscle alignment with 5 iterations.
3. Refine alignment using MACSE.
4. Replace frameshift(!) and stop codon with NNN using MACSE.
5. Concatenate files (if in multi-file mode).
6. Build tree with raxml `-f a -x 42 -p 42 -m GTRGAMMA `.
7. Co-filter fasta file and tree file. Remove codon columns with more than 70% species missed, and remove species with more than 70% codons as "NNN". Trim and annotate tree with the foreground information provided. Output alignment as phylip format.
8. Generate codeml configuration files for both branch-site null model (omega=1) and alternative model.
9. Run both codeml model.
10. Generate p values and other statistics using scipy.


# Disclaimer
This package is not liable for any ouput results and scientific/commercial interpretation (see the MIT license).
