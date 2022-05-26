chr=$1

awk '{print $1, $2+1}' OFS="\t"  aMel_chr${chr}.wga.bed > FocalChrPos${chr}.txt
grep -v Amel FocalChrPos${chr}.txt > tmp
grep -v Acer tmp > tmp1
grep -v Ador tmp1 > tmp2
grep -v Aflo tmp2 > tmp3
mv tmp3 FocalChrPos${chr}.txt
rm tmp tmp1 tmp2
