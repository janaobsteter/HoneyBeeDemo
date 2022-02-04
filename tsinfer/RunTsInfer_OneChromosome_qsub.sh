#!/bin/bash
# Grid Engine options (lines prefixed with #$)
#$ -N RunTsinfer_CHROMOSOME_
#$ -cwd
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o RunTsinfer_CHROMOSOME_.txt
#$ -l h_rt=30:00:00
#$ -j yes
#$ -P roslin_HighlanderLab
#  These options are:
#  job name: -Cut_GENOTYPES
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
module load anaconda
source activate jana_tsinfer

bash RunTsInfer_OneChromosome.sh _CHROMOSOME_ &
./cpumemlog.sh $! > cpumemlog_CHROMOSOME_.txt
