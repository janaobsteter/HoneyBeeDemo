for cycle in $(seq 0 199)
do
	sed "s/_CYCLE_/${cycle}/g" RunEstSfsProgram_rate6.sh > RunEstSfsProgram${cycle}_rate6.sh
	qsub RunEstSfsProgram${cycle}_rate6.sh
done	
