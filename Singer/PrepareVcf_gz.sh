vcfDir=$1
vcf=$2
bcftools view -H ${vcfDir}/${vcf}.vcf.gz | cut -f10- > Genotypes_${vcf}.txt
sed -i "s/0/0|0/g" Genotypes_${vcf}.txt
sed -i "s/1/1|1/g" Genotypes_${vcf}.txt
bcftools view -H ${vcfDir}/${vcf}.vcf.gz | cut -f1-9 > Info_${vcf}.txt
bcftools view -h ${vcfDir}/${vcf}.vcf.gz > Header_${vcf}.txt
paste Info_${vcf}.txt Genotypes_${vcf}.txt > Body_${vcf}.txt
cat Header_${vcf}.txt Body_${vcf}.txt > ${vcf}_DH.vcf
bcftools view -c 1  ${vcf}_DH.vcf > ${vcf}_DH_polymorphic.vcf 

rm  Info_${vcf}.txt Body_${vcf}.txt Header_${vcf}.txt Genotypes_${vcf}.txt
