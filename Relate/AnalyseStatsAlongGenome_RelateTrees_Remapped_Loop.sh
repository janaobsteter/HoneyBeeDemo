for chr in $(seq 1 16)
do
	python AnalyseStatsAlongGenome_RelateTrees_Remapped.py ${chr}
done
