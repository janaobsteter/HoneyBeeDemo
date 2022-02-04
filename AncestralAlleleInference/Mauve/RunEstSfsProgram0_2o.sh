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
./est-sfs config-kimura.txt  EstSfs_Dict0E.csv  est-sfs-release-2.03/seedfile.txt EstsfsOutput_2o/outputEstsfs0.txt EstsfsOutput_2o/output-0-pvalues.txt
