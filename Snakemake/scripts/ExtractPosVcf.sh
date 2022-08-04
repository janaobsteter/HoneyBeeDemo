vcf=$1
alignedPos=$2
output=$3
outputPos=$4

bcftools 	--gzvcf $vcf
	        --out $3
	        --positions $alignedPos
	        --get-INFO AC
	        --get-INFO AF
