#!/bin/sh
# Grid Engine options (lines prefixed with #$
#$ -N vcf2geneticMap
#$ -cwd                  
#  These options are:
#  job name: -N
#  use the current working directory: -cwd

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load modules
module load roslin/R/4.0.0
module load igmm/apps/vcftools/0.1.13

# Path to R script to generate genetic maps
SCRIPT=vcf2geneticMap.R
# VCF file is ancestral.imputed.vcf.gz
VCF=
# GEN is output prefix for genetic map file
GEN=

while getopts ":g:v:" opt; do
  case $opt in
    g) GEN=${OPTARG};;    
    v) VCF=${OPTARG};;
  esac
done

if [[ -z ${GEN} ]] | [[ -z ${VCF} ]]
then
  echo "Error"
  exit 1
fi

# Convert VCF to Plink format files
vcftools --gzvcf ${VCF} --chr ${SGE_TASK_ID} --plink --out ${VCF%.gz}_${SGE_TASK_ID}
rm ${VCF%.gz}_${SGE_TASK_ID}.ped
rm ${VCF%.gz}_${SGE_TASK_ID}.log

# Extract columns 1 and 4 to file
cut -f1,4 ${VCF%.gz}_${SGE_TASK_ID}.map > ${VCF%.gz}_${SGE_TASK_ID}.tmp
mv ${VCF%.gz}_${SGE_TASK_ID}.tmp ${VCF%.gz}_${SGE_TASK_ID}.map

# Pass file to R script to generate map
R CMD ${SCRIPT} ${VCF%.gz}_${SGE_TASK_ID}.map ${GEN}_LG${SGE_TASK_ID}.txt
