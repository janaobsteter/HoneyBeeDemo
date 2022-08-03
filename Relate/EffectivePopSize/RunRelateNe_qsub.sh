#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N RelateChr_CHR_Ne
#$ -cwd       
#$ -l h_vmem=80G
#$ -pe sharedmem 1
#$ -o ./Relatechr_CHR_Ne.txt
#$ -l h_rt=25:00:00
#$ -j yes
#$ -P roslin_HighlanderLab
#$ -l hostname='node3*.ecdf.ed.ac.uk'

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load intel/2020u4
module load roslin/R/4.1.0
module load GCC/10.3.0 
module load anaconda
source activate jana_tsinfer
 
# Run the program
PATH_TO_RELATE="/home/v1jobste/jobsteter/Honeybees/Relate/relate_v1.1.8_x86_64_dynamic/"
GENMAPS="/home/v1jobste/jobsteter/Honeybees/Relate/RelateGenMaps/"


${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh   \
	      -i output/Chr_CHR_ \
              -m 3.4e-9 \
              --poplabels SamplePopGroup.csv \
              --seed 1 \
	      --years_per_gen 1 \
              -o Chr_CHR__Ne &

./cpumemlog.sh $! Cpumemlog_CHR__Ne.txt
