#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N CreateInputForEstsfs_0
#$ -cwd
#$ -l h_vmem=4G
#$ -o CreateInputForEstsfs_0.txt
#$ -l h_rt=10:00:00
#$ -j yes
#$  -P roslin_HighlanderLab
#  These options are:
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem



# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Anaconda
module load anaconda
source activate jana_tsinfer

# Run the program
python CreateInputForEstsfs.py 0
