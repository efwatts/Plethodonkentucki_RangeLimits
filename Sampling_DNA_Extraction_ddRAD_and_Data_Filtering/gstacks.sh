#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=48:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW

# run ustacks
~/bin/gstacks \
-P /fs/scratch/PHS0328/kentucki_EFW/stacks.denovo/stacks.M12/ \
-M /fs/scratch/PHS0328/kentucki_EFW/scripts/EFW_popmap_all.tsv \
-O /fs/scratch/PHS0328/kentucki_EFW/stacks.denovo/stacks.M12/ \
-t 40