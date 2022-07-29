module load anaconda
source activate jana_tsinfer

python ComputeWeightedStatistics.py "/home/v1jobste/jobsteter/Honeybees/tsinfer/OnlyAncestralInference/"  "/exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/Honeybees/GenomicData/Wragg/MetaData/Meta.csv"
