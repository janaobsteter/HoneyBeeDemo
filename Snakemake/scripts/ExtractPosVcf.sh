vcf=$1
alignedPos=$2
output=$3

bcftools 	--gzvcf $vcf
	        --out VcfInfo
	        --positions $alignedPos
	        --get-INFO AC
	        --get-INFO AF

awk '{print $1 "_" $2}' VcfInfo.INFO > $3
