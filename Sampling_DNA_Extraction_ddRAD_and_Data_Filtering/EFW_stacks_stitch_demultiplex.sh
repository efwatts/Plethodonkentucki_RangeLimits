## Performing analyses on Ohio OSC HPC ##


## Pipeline for PE ddRAD:
# 1. Check read quality with fastQC
# 2. Stitch PE reads 1 + 2 into single 300bp locus
# 3. Demultiplex with stacks
# 4  Optimize stacks flitering parameters (R80)
# 5. Filter all samples 


## ---- 1. Quality check with fastQC: ----

# Log on to Ohio HPC through the GUI
# fastQC already downloaded there, including a GUI

# ** I ended up running fastQC on my own PC, with files stored on hard drive **


## ---- 2. Stitch reads ----

# navigate to scratch folder on OSC HPC
cd /fs/scratch/PHS0328/kentucki_EFW/library1

# for stacks, first make directory structure (Rochette & Catchen 2017)
mkdir ./alignments/
mkdir ./cleaned/
mkdir ./genome/
mkdir ./info/
mkdir ./raw/
mkdir ./stacks.denovo/
mkdir ./stacks.ref/
mkdir ./tests.denovo/
mkdir ./tests.ref/

## ---- 2. Gather raw data & Stitch PE reads ----

# concatenation (stitch) of PE reads - borrowed from Hime et al. 2018; Murphy et al. 2019

# copy all raw read files to common directory (adjust paths as necessary)
# navigate to raw reads directory
cd ./ddRAD_rawfiles/

# unzip new library
gunzip -c EFW_lib1_CKDL220004510-1a_HHJNFBBXX_L2_1.gz > EFW_lib1_R1
gunzip -c EFW_lib1_CKDL220004510-1a_HHJNFBBXX_L2_2.gz > EFW_lib1_R2

# concatenate all Read 1 files from all libraries; repeat with Read 2
#Don't need this EFW
cat EFW_lib1_R1 > R1
cat EFW_lib1_R2 > R2

# submit job that will fix indexes, reverse stitch R1 to R2
# index fixing removes the dual index and replaces with the correct 6bp sequence
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/library1/scripts/EFW_stitch_reads.sh


## ---- 3. Demultiplex with stacks ----
# set working directory
cd /fs/scratch/PHS0328/kentucki_EFW/library1

# barcode dist1
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/library1/scripts/demultiplex_barcode.dist1.sh ../ddRAD_rawfiles/EFW_lib1_CKDL220004510-1a_HHJNFBBXX_L2_1.fq.gz ../ddRAD_rawfiles/EFW_lib1_CKDL220004510-1a_HHJNFBBXX_L2_2.fq.gz ../info/EFW_barcodes.txt

#after this step, look at the process_radtags log file to see which ones succeeded vs. failed (success is >900,000)
#remove failed from EFW_popmap_test (and other popmap files)

## ---- 4. Process with stacks: optimize parameters ----

## Using a subset of data to optimize filtering parameters (R80 method)

########################ALL SAMPLES#########################################


##Doing populations on librarys 1 and 2 
cd /fs/scratch/PHS0328/kentucki_EFW

sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/all.scripts/populations_structure.sh


## ---- 4. Process with stacks ----

## Using a subset of data to optimize filtering parameters (R80 method)

# set working directory
cd /fs/scratch/PHS0328/kentucki_EFW

for M in {1..12}; do
mkdir -p ./stacks.denovo/stacks.M${M}
done

for run in {1..12}; do
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/denovo_n20.sh ${run}
done

#do this to determine which M is the best to use
## count number of r80 loci
for i in {1..12}; do
cat ./stacks.denovo/stacks.M${i}/populations.sumstats.tsv | \
    grep -v '^#' |
    cut -f 1 |
    sort -n -u |
    wc -l
done


###just copied the below from Brian's

## ---- 5. Process complete data set ----

## Using M = n = 12, run pipeliline:
# n=9 was the best fit for library 1, but n=12 was for library 2....and together, library 3,
#n = 3 and 9 resulted in very few snps

# 1. ustacks - build loci within individuals
# 2. cstacks - assembles catalog of loci across individuals
# 3. sstacks - match samples to catalog
# 4. tsv2bam - convert files
# 5. gstacks - genotype individuals
# 6. populations - compute sumstats; export files


## USTACKS

## test on a single sample
###CHANGE everything to M12 (M2 is shown all together to be the best, but when I tried M3 and M9, M9 was better)

# testing on single sample 
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/ustacks_individual.sh

# success! Now, loop through sample names; first one at index 0 (DH_52791) already done above

# par1 = individual samples
# par2 = sample index
# par3 = M=n

# array of individual sample names
cd /fs/scratch/PHS0328/kentucki_EFW
sample=(DH_52791 DH_65590 DH_57064 EFW_0006 EFW_0016 EFW_0025 SRK_3200 SRK_3186 DH_65668 EFW_0007 EFW_0017 EFW_0026 SRK_3204 DH_60223 DH_69838 EFW_0009 EFW_0018 SRK_2975 DH_60235 DH_77587 EFW_0010 EFW_0020 SRK_3025 DH_60245 DH_77595 EFW_0011 EFW_0021 SRK_3027 DH_60255 DH_78196 EFW_0012 EFW_0022 DH_60257 DH_78217 EFW_0014 EFW_0024 SRK_3198)
for i in {1..37}
do
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/ustacks.sh ${sample[${i}-1]} $i 12
done

## CSTACKS
cd /fs/scratch/PHS0328/kentucki_EFW

sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/cstacks.sh


## SSTACKS
cd /fs/scratch/PHS0328/kentucki_EFW/
sample=(DH_52791 DH_65590 DH_57064 EFW_0006 EFW_0016 EFW_0025 SRK_3200 SRK_3186 DH_65668 EFW_0007 EFW_0017 EFW_0026 SRK_3204 DH_60223 DH_69838 EFW_0009 EFW_0018 SRK_2975 DH_60235 DH_77587 EFW_0010 EFW_0020 SRK_3025 DH_60245 DH_77595 EFW_0011 EFW_0021 SRK_3027 DH_60255 DH_78196 EFW_0012 EFW_0022 DH_60257 DH_78217 EFW_0014 EFW_0024 SRK_3198)
for i in {1..37}
do
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/sstacks.sh ${sample[${i}-1]}
done


## TSV2BAM

cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/tsv2bam.sh


## GSTACKS
cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/gstacks.sh


## POPULATIONS STRUCTURE
cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/populations.structure.sh

## POPULATIONS FASTA
cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/populations.fasta.sh

##POPULATIONS VCF
cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/populations.vcf.sh

##POPULATIONS FASTA of SNPs
cd /fs/scratch/PHS0328/kentucki_EFW/
sbatch --account=PHS0328 /fs/scratch/PHS0328/kentucki_EFW/scripts/populations.fasta.snp.sh






# (** if desired, can continue from here with ipyrad instead **)


# Can run Denovo_map.pl – this will run following programs: (I did this in line 92)

# 1. ustacks - build loci
# 2. cstacks - assembles catalog of loci
# 3. sstacks - match samples to catalog
# 4. tsv2bam - convert files
# 5. gstacks - genotype individuals
# 6. populations - compute sumstats; export files
        ##populations can output FASTA and treemix and vcf and genepop and plink and structure and phylip and more
