#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N RelateChr_CHR_
#$ -cwd       
#$ -l h_vmem=80G
#$ -pe sharedmem 1
#$ -o ./Relatechr_CHR_.txt
#$ -l h_rt=10:00:00
#$ -j yes
#$ -P roslin_HighlanderLab
#$ -l hostname='node3*.ecdf.ed.ac.uk'

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load intel/2020u4
module load roslin/R/4.1.0
module load GCC/10.3.0 
 
# Run the program
PATH_TO_RELATE="/home/v1jobste/jobsteter/Honeybees/Relate/relate_v1.1.8_x86_64_dynamic/"
GENMAPS="/home/v1jobste/jobsteter/Honeybees/Relate/RelateGenMaps/"

${PATH_TO_RELATE}/bin/Relate --mode All \
        --haps ./AncestralChr_CHR__haploid.haps \
        --sample ./Chr_CHR__missing.sample \
        --map $GENMAPS/Chr_CHR_.gmap \
        -N 200000 \
        -m 3.4e-9 \
        -o Chr_CHR_ \
        --seed 1 &

./cpumemlog.sh $! > cpumemlogChr_CHR_
