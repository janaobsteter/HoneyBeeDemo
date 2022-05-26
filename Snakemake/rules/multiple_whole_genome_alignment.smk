configfile: "config/ancestral.yaml"

rule_all:
    input:
        "OutspeciesInfo_All_aligned.txt"

rule multiple_genome_alignment:
    input:
        config=config["config"]
    output:
        "evolverApis.hal"
    shell:
        "scripts/RunCactus.sh" #cactus apisAlignment {input.config} {output}
#These files are on Andraz's computer > bin> cactus >Honeybees

rule hal_to_maf:
    input:
        "evolverApis.hal"
    output:
        "apisCactus.maf.gz"
    shell:
        "docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.0.3 hal2maf "
        "{input}  --refGenome Amel {output}"


rule maf_to_bed:
    input:
        maf="apisCactus.maf.gz"
        target="Amel"
    output:
        expand("aMel_chr{chromosome}.wga.bed", chromosome = config['chromosomes'])
    conda:
        "envs/py27.yaml"
    script:
        "scripts/maf_to_bed.py -i {input.maf} -r {input.target} -c {wildcards.chromosome} | "
        "sort -k1,1 -k2,2n -u > aMel_chr${wildcards.chromosome}.wga.bed -"

rule extract_focal_positions_aligned:
    input:
        expand("aMel_chr{chromosome}.wga.bed", chromosome = config['chromosomes'])
    output:
        "FocalChrPos{chromosome}.txt"
    script:
        "bash scripts/GetPosWgaBed.sh {wildcards.chromosome}"

# rule combine_focal_positions_aligned:
#     input:
#         expand("FocalChrPos{chromosome}.txt", chromosome = config['chromosomes'])
#     output:
#         "FocalChrPos_All.txt"
#     script:
#         "bash scripts/CombinePos.sh config['chromosomes']"


rule extract_alleles_bed:
    input:
        expand("aMel_chr{chromosome}.wga.bed", chromosome = config['chromosomes'])
    output:
        "OutspeciesInfo{chromosome}.txt"
    script:
        "scripts/ExtractAllelesBed.sh {wildcards.chromosome}"

rule combine_alleles:
    input:
        expand("OutspeciesInfo{chromosome}.txt", chromosome = config['chromosomes'])
    output:
        alleles="OutspeciesInfo_All_aligned.txt"
        positions="FocalChrPos_All_aligned.txt"
    script:
        "scripts/CombineAlleles.sh config['chromosomes']"
