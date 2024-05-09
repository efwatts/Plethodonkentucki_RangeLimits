#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=1:00:00
#SBATCH --ntasks=40


# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW/library1

# absolute path to working direcrtory
work=/fs/scratch/PHS0328/kentucki_EFW/library1

# set file for popmap, folder for demultiplexed reads, output folder, and log file
popmap=$work/info/EFW_popmap.tsv

# Output directory with results of denovo.pl until gstacks step
out=$work/tests.denovo/stacks.M9

# Move into the directory
cd $out

~/bin/populations \
-P $out \
-M $popmap \
-t 40 \
-r 0.8
