for cycle in $(seq 0 199)
do
	sed "s/_CYCLE_/${cycle}/g" RunEstSfsProgram_2o.sh > RunEstSfsProgram${cycle}_2o.sh
	qsub RunEstSfsProgram${cycle}_2o.sh
done	
