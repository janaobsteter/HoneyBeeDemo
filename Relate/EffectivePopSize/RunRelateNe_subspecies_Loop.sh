subspecies=$1
for chr in $(seq 1 1)
do
	sed "s/_SUBSPECIES_/${subspecies}/g" RunRelateNe_subspecies_qsub.sh > RunRelate${chr}Ne_${subspecies}_qsub.sh
	sed -i "s/_CHR_/${chr}/g" RunRelate${chr}Ne_${subspecies}_qsub.sh 
	qsub RunRelate${chr}Ne_${subspecies}_qsub.sh
done
