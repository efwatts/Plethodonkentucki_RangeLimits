#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=5:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW

# run sstacks
~/bin/sstacks \
-c ./stacks.denovo/stacks.M12/ \
-s ./stacks.denovo/stacks.M12/$1 \
-o ./stacks.denovo/stacks.M12/ \
-p 40
