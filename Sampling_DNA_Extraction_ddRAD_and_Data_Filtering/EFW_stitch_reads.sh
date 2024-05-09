#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=10:00:00
#SBATCH --ntasks=40

# setwd
cd /fs/scratch/PHS0328/kentucki_EFW/library1/ddRAD_rawfiles

# reduce
sed 's/NTCACGAT+NGATCTCG/ATCACG/g' R1.fq > R1_x.fq
sed 's/NTCACGAT+NGATCTCG/ATCACG/g' R2.fq > R2_x.fq

# reverse complement R2

awk 'NR%4==2' R2_x.fq | awk '{ print "\n\n\n"$1;}' | sed -e '1,2d' | rev | tr ACGT TGCA > R2.seq.rc

# reverse the quality scores for R2

awk 'NR%4==0' R2_x.fq | awk '{ print "\n\n\n"$1;}' | rev > R2.qual

# stitch R1 to rev. R2

paste R1_x.fq R2.seq.rc R2.qual | awk '{gsub("\t","",$0); print $0;}' | awk '{gsub(" 1:N:0:"," stitched:N:0:",$0); print $0;}' > R1_R2.stitched.fq


