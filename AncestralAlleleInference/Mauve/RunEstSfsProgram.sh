#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N RunEstsfs__CYCLE_
#$ -cwd
#$ -l h_vmem=16G
#$ -o RunEstsfsProgram__CYCLE_.txt
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
./est-sfs config-kimura.txt  EstSfs_Dict_CYCLE_E.csv  est-sfs-release-2.03/seedfile.txt EstsfsOutput2o/outputEstsfs_CYCLE_.txt EstsfsOutput2o/output-_CYCLE_-pvalues.txt
