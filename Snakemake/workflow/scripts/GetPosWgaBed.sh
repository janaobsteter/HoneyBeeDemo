chr=$1

awk '{print $1, $2+1}' OFS="\t"  aMel_chr${chr}.wga.bed > FocalChrPos${chr}.txt
grep -v Amel FocalChrPos${chr}.txt |
grep -v Acer tmp |
grep -v Ador tmp1 |
grep -v Aflo tmp2 |
mv tmp3 FocalChrPos${chr}.txt
