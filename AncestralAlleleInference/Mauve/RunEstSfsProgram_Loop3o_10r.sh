for cycle in $(seq 0 199)
do
	sed "s/_CYCLE_/${cycle}/g" RunEstSfsProgram_3o_10r.sh > RunEstSfsProgram${cycle}_3o_10r.sh
	qsub RunEstSfsProgram${cycle}_3o_10r.sh
done	
