for cycle in $(seq 0 199)
do
	echo $cycle
	Rscript DetermineAncAllel_Pos.R $cycle
done
