#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=5:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/BPW_cinereus_stacks

# run ustacks
/users/PHS0328/bw120818/bin/populations \
-P ./stacks.denovo/stacks.M11 \
-M /fs/scratch/PHS0328/BPW_cinereus_stacks/info/popmap_155.tsv \
-R 0.5 \
-t 40 \
--structure \
--write-single-snp

