
vcfDir=$1
vcf=$2

bcftools view -H ${vcfDir}/${vcf}.vcf | cut -f10- > Genotypes_${vcf}.txt
bash JoinCol.sh Genotypes_${vcf}.txt > Genotypes_diploid_${vcf}.txt

numSamples=$(head -n 1 Genotypes_diploid_${vcf}.txt | awk '{print NF}')
colSamples=$(($numSamples + 9))

echo "Cutting info"
bcftools view -H ${vcfDir}/${vcf}.vcf | cut -f1-9 > Info_${vcf}.txt
echo "Cutting header"
bcftools view -h ${vcfDir}/${vcf}.vcf | head -n -1 > Header_${vcf}.txt
bcftools view -h ${vcfDir}/${vcf}.vcf | tail -n1 | cut -f1-${colSamples} > Names_${vcf}.txt
paste Info_${vcf}.txt Genotypes_diploid_${vcf}.txt > Body_${vcf}.txt
awk 'BEGIN {OFS = FS = "\t"} NF {NF -= 1} 1' <Names_${vcf}.txt  > Names_${vcf}_noLastColumn.txt
awk 'BEGIN {OFS = FS = "\t"} NF {NF -= 1} 1' <Body_${vcf}.txt > Body_${vcf}_noLastColumn.txt
cat Header_${vcf}.txt Names_${vcf}_noLastColumn.txt > FullHeader_${vcf}.txt
cat  FullHeader_${vcf}.txt Body_${vcf}_noLastColumn.txt  > ${vcf}_diploid.vcf
echo "Excluding non polymorphic"
bcftools view -c 1  ${vcf}_diploid.vcf > ${vcf}_diploid_polymorphic.vcf

rm  Info_${vcf}.txt Body_${vcf}.txt Header_${vcf}.txt Genotypes_${vcf}.txt ${vcf}_diploid.vcf
