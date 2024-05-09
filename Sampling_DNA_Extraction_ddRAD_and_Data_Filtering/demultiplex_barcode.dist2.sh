#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=10:00:00
#SBATCH --ntasks=40

# setwd
cd /fs/scratch/PHS0328/kentucki_EFW/library1/

# demultiplex
~/bin/process_radtags -f ./ddRAD_rawfiles/R1_R2.stitched.fq --inline_index -b ./info/EFW_barcodes.txt -o ./cleaned_barcode.dist2 --renz_1 sphI --renz_2 ecoRI -c -q -D -w 0.15 -s 20 -r --barcode_dist_1 2
