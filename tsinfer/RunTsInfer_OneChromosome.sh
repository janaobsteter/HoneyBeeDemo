chromosome=$1
vcf=/exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/Filtering_vcf/data/pipeline/hi_quality_removed_het/HiFiVCF/Chr${chromosome}.vcf
meta=/exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/Honeybees/tsinfer/Meta.csv
anc=/exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/Honeybees/tsinfer/AncestralAlleles_HiFiVcf.csv

python TsinferDroneInference_Eddie.py $chromosome $vcf $meta $anc
