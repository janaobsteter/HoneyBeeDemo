
chr=$1
cut -f5 aMel_chr${chr}.wga.bed | cut -f1,2,3,4 -d"," |  grep -v -w "'Amel'" | awk 'NF'> OutspeciesNames${chr}.txt
cut -f8 aMel_chr${chr}.wga.bed | cut -f1,2,3,4 -d"," |  grep -v -w "'Amel'" | awk 'NF' > OutspeciesAlleles${chr}.txt
paste FocalChrPos${chr}.txt OutspeciesNames${chr}.txt OutspeciesAlleles${chr}.txt > OutspeciesInfo${chr}.txt
