configfile: "config/multiple_whole_genome_alignment.yaml"
chromosomeDict = {'1': 'NC_037638.1',
                  '2': 'NC_037639.1',
                  '3': 'NC_037640.1',
                  '4': 'NC_037641.1',
                  '5': 'NC_037642.1',
                  '6': 'NC_037643.1',
                  '7': 'NC_037644.1',
                  '8': 'NC_037645.1',
                  '8': 'NC_037646.1',
                  '10': 'NC_037647.1',
                  '11': 'NC_037648.1',
                  '12': 'NC_037649.1',
                  '13': 'NC_037650.1',
                  '14': 'NC_037651.1',
                  '15': 'NC_037652.1',
                  '16': 'NC_037653.1',
                  'MT': 'NC_001566.1',
}

def mapChromosomeName(wildcards):
    return(chromosomeDict[wildcards.chromosome])

rule all:
    input:
        "OutspeciesInfo_All_aligned.txt"

rule multiple_genome_alignment:
    input:
        config=config["cactusConfig"]
    output:
        "MultipleGenomeAlignment/evolverApis.hal"
    shell:
        "docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.0.3 "
        "cactus apisAlignment {input.config} {output}"
#These files are on Andraz's computer > bin> cactus >Honeybees

rule hal_to_maf:
    input:
        "MultipleGenomeAlignment/evolverApis.hal"
    output:
        "MultipleGenomeAlignment/apisCactus.maf.gz"
    shell:
        "docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.0.3 "
        "hal2maf {input}  --refGenome Amel {output}"


rule maf_to_bed:
    input:
        maf="MultipleGenomeAlignment/apisCactus.maf.gz",
        target="Amel",
        chrName=mapChromosomeName
    output:
        expand("MultipleGenomeAlignment/aMel_chr{chromosome}.wga.bed", chromosome = config['chromosomes'])
    conda:
        "envs/py27.yaml"
    script:
        "scripts/maf_to_bed.py -i {input.maf} -r {input.target} -c {input.chrName} | "
        "sort -k1,1 -k2,2n -u > aMel_chr${wildcards.chromosome}.wga.bed -"

rule extract_focal_positions_aligned:
    input:
        "MultipleGenomeAlignment/aMel_chr{chromosome}.wga.bed"
    output:
        "MultipleGenomeAlignment/FocalChrPos{chromosome}.txt"
    script:
        "bash scripts/GetPosWgaBed.sh {wildcards.chromosome}"

# rule combine_focal_positions_aligned:
#     input:
#         expand("FocalChrPos{chromosome}.txt", chromosome = config['chromosomes'])
#     output:
#         "FocalChrPos_All.txt"
#     script:
#         "bash scripts/CombinePos.sh config['chromosomes']"


rule extract_alleles_bed: #IS IT OK TO HAVE EXPAND HERE?
    input:
        "MultipleGenomeAlignment/aMel_chr{chromosome}.wga.bed"
    output:
        "MultipleGenomeAlignment/OutspeciesInfo{chromosome}.txt"
    script:
        "scripts/ExtractAllelesBed.sh {wildcards.chromosome}"

rule rename_outspecies_info:
    input:
        "MultipleGenomeAlignment/OutspeciesInfo{chromosome}.txt"
    output:
        "MultipleGenomeAlignment/OutspeciesInfo{chromosome}_renamed.txt"
    shell:
        "awk -v NUM={wildcards.chromosome} -v OFS='\t' '{print NUM,$2,$3,$4}' {input} > {output}"

rule combine_alleles:
    input:
        expand("MultipleGenomeAlignment/OutspeciesInfo{chromosome}_renamed.txt", chromosome = config['chromosomes'])
    output:
        alleles="MultipleGenomeAlignment/OutspeciesInfo_All_aligned.txt",
        positions="MultipleGenomeAlignment/FocalChrPos_All_aligned.txt"
    script:
        "scripts/CombineAlleles.sh config['chromosomes']"
