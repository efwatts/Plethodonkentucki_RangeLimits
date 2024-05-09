#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=5:00:00
#SBATCH --ntasks=40


# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW

# set variables
sample=DH_52791
index=1
M=12

# run ustacks
~/bin/ustacks \
-t gzfastq \
-f ./raw_demultiplexed/${sample}.1.fq.gz \
-o ./stacks.denovo/stacks.M12/stacks.M${M} \
-i ${index} \
--name ${sample} \
-M ${M} \
-p 40