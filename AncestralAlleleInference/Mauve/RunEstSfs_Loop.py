for cycle in $(seq 1 200)
do
	sed "s/_CYCLE_/${cycle}/g" RunEstSfs.sh > RunEstSfs${cycle}.sh
	qsub RunEstSfs${cycle}.sh
done
