from os import listdir
from pathlib import Path
configfile: "../config/tsinfer.yaml"

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]


rule all:
    input:
        vcf=expand("Tsinfer/Chr{chromosome}_ancestral.vcf.gz", chromosome = range(1, config['noChromosomes'] + 1)),
        info=expand("Tsinfer/Info{chromosome}.INFO", chromosome = range(1, config['noChromosomes'] + 1)),

workdir: config['workdir']

VCFs =[x for x in list_full_paths(config['vcfDir']) if x.endswith(".vcf.gz")]
noVCFs = len(VCFs)
#if noVCFs == config['noChromosomes']:
    #print("Your VCFs are split.")
#
#if (noVCFs != config['noChromosomes']) and noVCFs != 1:
    #print("Check you VCFs.")

if noVCFs == 1:
    #print("Splitting VCFs.")
    rule split_vcfs: #Input VCFs need to be compressed and indexes
        input:
            vcf=VCFs[0]
        output:
            vcf=[config['vcfDir'] + x for x in expand("Chr{{chromosome}}.vcf.gz")]
        params:
            vcfDir=config['vcfDir']
        #envmodules:
        #    config['bcftoolsModule']
        shell:
            """
            bcftools view -r {wildcards.chromosome} -O z {input} > {output.vcf}
            bcftools index {output.vcf}
            """

rule get_major:
    input:
        config['vcfDir'] + "Chr{chromosome}.vcf.gz"
    output:
        infoVcf="Tsinfer/tmp.gz",
        info="Tsinfer/Info{chromosome}.INFO"
    params:
        prefix="Tsinfer/Info{chromosome}"
    shell:
        """
        bcftools +fill-tags {input} -Oz -o {output.infoVcf} -- -t AC,AF
        vcftools --gzvcf {output.infoVcf} --out {params.prefix} --get-INFO AC --get-INFO AF --get-INFO REF
        """

rule decompress:
    input:
        config['vcfDir'] + "Chr{chromosome}.vcf.gz" #This will take
    output:
        config['vcfDir'] + "Chr{chromosome}.vcf"
    shell:
        """
        gunzip {input}
        """


rule extract_vcf_pos:
    input:
        rules.decompress.output
    output:
        "Tsinfer/VcfPos{chromosome}.txt"
    #conda:
    #    config['tsinferEnv']
    #envmodules:
    #    config['bcftoolsModule']
    params:
        vcfDir=config['vcfDir']
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > tmp
        awk '{{print $1"_"$2}}' tmp > {output}
        """

rule match_ancestral_vcf:
    input:
        vcfPos=rules.extract_vcf_pos.output, # This has more lines,
        ancestralAllele=config['ancestralAllele'] # The ancestral file has to have chr_pos and AA, split with a tab
    output:
        "Tsinfer/AncestralVcfMatch{chromosome}.txt" # THis files needs to contain all the variants from the vcf (blank space)
    shell:
        """
        for line in $(cat {input.vcfPos});
        do
          grep $line {input.ancestralAllele} || echo "";
        done > tmp
        awk -F"," '{{print $1"\t"$2}}' tmp > {output}
        """

rule change_infoAA_vcf:
    input:
        vcf=rules.decompress.output,
        ancestralAllele=rules.match_ancestral_vcf.output
    output:
        "Tsinfer/Chr{chromosome}_ancestral.vcf"
    shell:
        """
        HEADERNUM=$(( $(grep "##" {input.vcf} | wc -l) + 1 ))
        INFOLINE=$(( $(grep -Fn "INFO" {input.vcf} | cut --delimiter=":" --fields=1  | head -n1) ))
        awk -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<=HEADER{{{{print}}}}; \
        FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\\n"}}}}; \
        FNR>HEADER{{{{$8=a[FNR-HEADER]; print}}}}' OFS="\t" {input.ancestralAllele} {input.vcf} > {output}
        """
rule compress_vcf:
    input:
        rules.change_infoAA_vcf.output
    output:
        "Tsinfer/Chr{chromosome}_ancestral.vcf.gz"
    shell:
        """
        bgzip {input}
        bcftools index {output}
        """
