#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=96:00:00
#SBATCH --ntasks=40


# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW/library1/stacks.denovo

# absolute path to working direcrtory
work=/fs/scratch/PHS0328/kentucki_EFW/library1

# set file for popmap, folder for demultiplexed reads, output folder, and log file
popmap=$work/info/EFW_popmap.tsv
reads_dir=$work/cleaned_barcode.dist2

# This creates a new output directory per Stacks run
out=$work/stacks.denovo/stacks.M$1
mkdir -p $out
# Move into the new directory
cd $out

~/bin/denovo_map.pl \
--samples $reads_dir \
--popmap $popmap \
--out-path $out \
-M $1 \
-n $1 \
-T 40 \
--min-samples-per-pop 0.8

#gstacks \
#-P $out \
#-M $popmap \
#-t 40
