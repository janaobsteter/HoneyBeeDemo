for cycle in $(seq 0 199)
do
	sed "s/_CYCLE_/${cycle}/g" RunEstSfsProgram_3o.sh > RunEstSfsProgram${cycle}_3o.sh
	qsub RunEstSfsProgram${cycle}_3o.sh
done	
