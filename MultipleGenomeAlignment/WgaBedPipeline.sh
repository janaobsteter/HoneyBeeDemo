#Get all alignments in one file
cat AlignmentChr[1-9].wga.bed > AlignmentAllChrs.wga.bed
cat AlignmentChr1[0-2].wga.bed >> AlignmentAllChrs.wga.bed
#Get aligned positions to filter the vcf
awk ‘{print $1, $2+1}’ OFS=“\t” AlignmentAllChrs.wga.bed > FocalChrPos.txt #wga is 0-based and the vcf is 1-based that is why I add 1 to position in wga
#awk ‘{print $0, $1"_“$2}’ OFS=“\t” AlignmentAllChrs.wga.bed > AlignmentAllChrsPos.wga.bed
### Link to vcfs with info
ln -s ~/HighlanderLab/1_VCFsToUse/AllFinalAncestralINFO.vcf
#Get vcf info of postions that are present in the alignment
vcftools --vcf AllFinalAncestralINFO.vcf --positions FocalChrPos.txt --get-INFO AC --get-INFO AF --out All_VCFSnps_Focal
#Filter wga.bed file no to have so large files
awk ‘NR==FNR{a[$1"_“$2-1];next} $1”_“$2 in a {print $0, $1"_“$2+1}’ OFS=“\t” All_VCFSnps_Focal.INFO AlignmentAllChrs.wga.bed > All_AlignmentAllChrsPos.wga.bed
./CreateInputForEstsfs_fromWGAbed.py All_AlignmentAllChrsPos.wga.bed All_VCFSnps_Focal.INFO All_est-sfs_inputFile.txt
#Run est-sfs
~/Programas/est-sfs-release-2.03/est-sfs config-file.txt All_est-sfs_inputFile.txt seed-file.txt output-All-sfs.txt output-All-pvalues.txt
#Get vcf to add est-sfs ancestral allele states
ln -s ~/HighlanderLab/1_VCFsToUse/Geno_2021_TODAS.vcf
#Run script to get est-sfs ancestral allele states
./ExtractAncestralAlleleFromEstsfsOutput.py output-All-pvalues.txt All_estsfsAncestralStates.txt
