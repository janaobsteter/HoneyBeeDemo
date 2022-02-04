for chr in 7
do 
	sed "s/_CHROMOSOME_/${chr}/g" RunTsInfer_OneChromosome_qsub.sh > RunTsInfer_Chromosome${chr}_qsub.sh
	qsub RunTsInfer_Chromosome${chr}_qsub.sh
done
