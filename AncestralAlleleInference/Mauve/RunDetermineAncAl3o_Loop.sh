for cycle in $(seq 0 199)
do
	echo $cycle
	python ExtractAncestralAlleleFromEstsfsOutput_3o.py ${cycle} 
done
