for chr in  $(seq 1 1)
do
	sed "s/_CHR_/${chr}/g" RunRelate_qsub.sh > RunRelate${chr}_qsub.sh
	qsub RunRelate${chr}_qsub.sh
done
