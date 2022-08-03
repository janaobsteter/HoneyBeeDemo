This is the pipeline to prepare the vcf files (by chromosome) to be suitable to run with Relate
1) Run PrepareRelateInputFiles.py
- This uses bcftools to transform .vcf to .haps, read in the file and add the ancestral allele information
2) Run EditHapsFile.py
- This computed haploid .haps from diploid .haps file, computes missingness and writes .samples
