popType=$1
pop=$2
foldType=$3
cutOff=$4

python GetSFS_Singer_whole_genome_KeepInds.py ${popType} ${pop} ${foldType} ${cutOff}
bash column_to_row.sh Genome_pop${pop}_unfold_count.sfs 
noPop=$(cat NoSamples_Pop${pop}.csv)
length=$(grep "^141516," /sietch_colab/jobste/Data/FilteringData/NewLength_Chromosomes.txt | cut -d"," -f2)
sed "s/ThisIsWhereSFSGoes/$(cat Genome_pop${pop}_unfold_count_one_row.txt)/"  honeybee_unfold_template.blueprint  > honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POPNO_/$(cat NoSamples_Pop${pop}.csv)/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POPNO1_/$(cat NoSamples_Pop1${pop}.csv)/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POPNO2_/$(cat NoSamples_Pop2${pop}.csv)/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POPNO3_/$(cat NoSamples_Pop3${pop}.csv)/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POPNO4_/$(cat NoSamples_Pop4${pop}.csv)/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_LENGTH_/${length}/" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_CHR_/Genome/g" honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i 's/SFS\\n/SFS/g'  honeybee_unfold_Genome_Pop${pop}.blueprint
sed -i "s/_POP_/${pop}/g"  honeybee_unfold_Genome_Pop${pop}.blueprint


java -cp stairway_plot_es Stairbuilder honeybee_unfold_Genome_Pop${pop}.blueprint
bash honeybee_unfold_Genome_Pop${pop}.blueprint.sh
