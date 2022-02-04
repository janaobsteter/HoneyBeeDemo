#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N RunEstsfs_0
#$ -cwd
#$ -l h_vmem=16G
#$ -o RunEstsfsProgram_0.txt
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
source activate jana_tsinfer

# Run the program
./est-sfs config-rate6.txt  EstSfs_Dict0E.csv  seedfile.txt EstsfsOutput_rate6/outputEstsfs0.txt EstsfsOutput_rate6/output-0-pvalues.txt
