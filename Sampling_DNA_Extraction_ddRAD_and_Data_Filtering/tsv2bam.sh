#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=12:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW

# run ustacks
~/bin/tsv2bam \
-P /fs/scratch/PHS0328/kentucki_EFW/stacks.denovo/stacks.M12/ \
-M /fs/scratch/PHS0328/kentucki_EFW/scripts/EFW_popmap_all.tsv \
-R /fs/scratch/PHS0328/kentucki_EFW/raw_demultiplexed/ \
-t 8