#!/bin/bash
#$ -N shapeit4
#$ -l h_rt=02:00:00
#$ -l h_vmem=16G
#$ -cwd
#$ -j yes
#$ -pe sharedmem 4
#$ -R y
# Initialize the environment modules
. /etc/profile.d/modules.sh

# Parse the input files and options
while getopts "i:g:r:o:c:" flag
do
	case "${flag}" in
		i) input_file=${OPTARG};;
		g) gmap_file=${OPTARG};;
		r) reference_file=${OPTARG};;
		o) output_file=${OPTARG};;
		c) region=${OPTARG};;
	esac
done

# Load anaconda
module load anaconda

# Load the shapeit4am envirionment
source activate shapeit4am

# Run the shapeit4 tool on the input file
shapeit4 \
--sequencing \
--thread 4 \
--region "$region" \
--input "$input_file" \
--map "$gmap_file"  \
--reference "$reference_file" \
--output "output/${output_file}" \
--log "${input_file}_shapeit4.log"
