include: ancestral_inference.smk
include: multiple_genome_alignment.smk
include: qc_vcf.smk

configfile: "config/tsinfer.yaml"

rule match_ancestral_vcf:


awk 'NR==FNR{a[NR] = $1; next} FNR<29{print}; FNR==5{printf "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"}; FNR>29{$8=a[FNR-29]; print}' OFS="\t" All_AAInfo.txt test.vcf awk '{FNR>28{print a[FNR-28]}}' All_AAInfo.txt test.vcf
