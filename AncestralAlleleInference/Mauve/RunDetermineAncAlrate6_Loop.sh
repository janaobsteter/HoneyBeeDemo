for cycle in $(seq 0 199)
do
	echo $cycle
	python ExtractAncestralAlleleFromEstsfsOutput_rate6.py ${cycle} 
done
