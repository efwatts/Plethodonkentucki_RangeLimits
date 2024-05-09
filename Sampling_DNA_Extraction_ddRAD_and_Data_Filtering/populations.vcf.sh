#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=1:00:00
#SBATCH --ntasks=40


# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW/stacks.denovo/stacks.M12

# absolute path to working direcrtory
work=/fs/scratch/PHS0328/kentucki_EFW

# set file for popmap, folder for demultiplexed reads, output folder, and log file
popmap=$work/scripts/popmap_35ind.tsv

# Output directory with results of denovo.pl until gstacks step
out=$work/stacks.denovo/stacks.M12

# Move into the directory
cd $out

~/bin/populations \
-P $out \
-M $popmap \
-t 40 \
-O /fs/scratch/PHS0328/kentucki_EFW/stacks.denovo/DroppedMissing/M12_R1_vcf \
--vcf
