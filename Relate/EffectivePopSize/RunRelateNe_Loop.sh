for chr in $(seq 1 1)
do
	sed "s/_CHR_/${chr}/g" RunRelateNe_qsub.sh > RunRelate${chr}Ne_qsub.sh
	qsub RunRelate${chr}Ne_qsub.sh
done
