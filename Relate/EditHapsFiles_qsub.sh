#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N EditHapsFile_CHR_
#$ -cwd
#$ -l h_vmem=20G
#$ -o EditHapsFile_CHR_.txt
#$ -l h_rt=02:00:00
#$ -j yes
#$  -P roslin_HighlanderLab
#  These options are:
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem



# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Anaconda
module load anaconda
module load igmm/apps/bcftools/1.9

# Activate environment
source activate jana_tsinfer
export LD_LIBRARY_PATH=/exports/cmvm/eddie/eb/groups/HighlanderLab/anaconda/envs/jana_tsinfer/lib:$LD_LIBRARY_PATH

# Edit the haps file
python EditHapsFile.py _CHR_  
