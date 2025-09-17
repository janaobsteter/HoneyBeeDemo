
chr=$1
PATH_TO_RELATE=/home/jobste/Relate/relate_lib/
${PATH_TO_RELATE}/bin/Convert \
                 --mode ConvertToTreeSequence \
                 --compress  \
                 --anc Masked_snps/Chr${chr}.anc \
                 --mut Masked_snps/Chr${chr}.mut \
                 -o Chr${chr}
