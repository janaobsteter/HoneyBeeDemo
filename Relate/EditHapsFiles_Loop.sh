for chr in $(seq 1 16)
do
	sed "s/_CHR_/${chr}/g" EditHapsFiles_qsub.sh > EditHapsFiles${chr}_qsub.sh
	qsub EditHapsFiles${chr}_qsub.sh
done
