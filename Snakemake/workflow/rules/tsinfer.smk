#include: ancestral_inference.smk
#include: multiple_genome_alignment.smk
#include: qc_vcf.smk
include prepare_files_for_tsinfer.smk

configfile: "../config/tsinfer.yaml"
workdir: config['workdir']

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1))

rule prepare_sample_file:
    input:
        vcf="Tsinfer/Chr{chromchromosomeosome}.vcf",
        meta=config['meta']
    output:
        "Tsinfer/Chr{chromosome}.samples"
    script:
        "scripts/PrepareTsinferSampleFile.py"

rule infer:
    input:
        "Tsinfer/Chr{chromosome}.samples"
    output:
        "Tsinfer/Chr{chromosome}.trees"
    script:
        "scripts/InferTrees.py"
