#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=4:00:00
#SBATCH --ntasks=40

# change to denovo folder 
cd /fs/scratch/PHS0328/cinereus_2022

# absolute path to working direcrtory
work=/fs/scratch/PHS0328/kentucki_EFW/library1

# set file for popmap, folder for demultiplexed reads, output folder, and log file
popmap=$work/info/EFW_popmap.tsv

# Output directory with results of denovo.pl until gstacks step
out=$work/stacks.denovo/stacks.M9

# Move into the directory
cd $out

# tsv2bam
#tsv2bam -P $out -M $popmap -t 40

# gstacks 
#gstacks -P $out -M $popmap -t 40


# populations
#mkdir -p $out/populations

~/bin/populations \
-P $out \
-M $popmap \
-t 40 \
-R 0.8 \
--write-single-snp \
--genepop \
--plink \
--structure
