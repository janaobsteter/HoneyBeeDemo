for cycle in $(seq 0 199)
do
	sed "s/_CYCLE_/${cycle}/g" CreateInputForEstsfs.sh > CreateInputForEstsfs${cycle}.sh
	qsub CreateInputForEstsfs${cycle}.sh
done
