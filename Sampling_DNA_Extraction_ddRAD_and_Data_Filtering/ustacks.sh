#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=10:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW

# run ustacks
~/bin/ustacks \
-t gzfastq \
-f ./raw_demultiplexed/$1.1.fq.gz \
-o ./stacks.denovo/stacks.M12/stacks.M$3 \
-i $2 \
--name $1 \
-M $3 \
-p 40
