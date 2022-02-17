#!/bin/bash 
#$ -l h_rt=12:00:00 
#$ -N filter_clear 
#$ -l h_vmem=8G
#$ -cwd
#$ -j yes
# Initialize the environment modules
. /etc/profile.d/modules.sh

# The needs the input files as arguments. It also needs the chromosomes file in 
# the same directory as it is located in.

# Load modules
module load igmm/apps/bcftools/1.13
module load roslin/plink/1.90p 

# Get the input filename and use the basename for all intermediate files
input_file="$1"
name=$(basename "$input_file" ".vcf.gz")

## Filter 1
bcftools view  -m2 -M2 -v snps "-iQUAL>30 && FMT/DP>5 & FMT/DP<(AVG(FMT/DP)*2) & FMT/GQ>30" "$input_file"  -O z -o "${name}_filter1.vcf.gz"  

## Filter 2
bcftools view -q "0.01:minor"  -Ou "${name}_filter1.vcf.gz" | bcftools view -i "F_MISSING < 0.10" -Ou | bcftools view -i "FMT/DP>3" -o  "${name}_filter2.vcf.gz"  

# THIS IS WHERE THE ERROR ECCOURS MOST LIKELY -- SUDDENLY 20 MILLION VARIANTS AGAIN
# Renaming the chromosomes
# chromoosmes.txt are the names of chromosomes paired with corresponding number
bcftools annotate --rename-chrs "chromosomes.txt" -O z -o "${name}_filter2_renamed.vcf.gz" "${name}_filter2.vcf.gz" 

## Index the file
bcftools index "${name}_filter2_renamed.vcf.gz"

## Subset to only variant on chromosomes
bcftools view -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,MT" -O z -o "${name}_filter2_renamed_subset.vcf.gz" "${name}_filter2_renamed.vcf.gz"

# Filter 3; filter using plink just using recode to get the proper file format
plink --vcf  "${name}_filter2_renamed_subset.vcf.gz" \
--allow-extra-chr \
--double-id \
--recode \
--out "${name}_filter_plink1" 

# Filter with plink
plink --file "${name}_filter_plink1" \
--allow-extra-chr \
--make-bed \
--mind "0.1" \
--geno "0.1" \
--out "${name}_filter_plink2"

# Recode to generate map file needed to conver to vcf or something
plink --bfile  "${name}_filter_plink2" \
--allow-extra-chr \
--recode \
--out "${name}_filter_plink3"

# Make the vcf
plink --file "${name}_filter_plink3" \
--recode "vcf" \
--vcf-iid \
-allow-extra-chr \
--out "${name}_pass"
