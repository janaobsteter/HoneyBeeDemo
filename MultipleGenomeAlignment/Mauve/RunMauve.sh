MAUVEDIR="/home/v1jobste/jobsteter/Honeybees/mauve_snapshot_2015-02-13/linux-x64/"
HOMEDIR="/home/v1jobste/jobsteter/Honeybees/MultipleGenomeAlignment"
REFDIR="/home/v1jobste/jobsteter/Honeybees/MultipleGenomeAlignment/ReferenceGenomes/"


$MAUVEDIR/progressiveMauve --output multiBees_long.xmfa  --output-guide-tree=multiBees_long.tree --backbone-output=multiBees_long.backbone \
		 $REFDIR/Amel/GCF_003254395.2_Amel_HAv3.1_genomic.fna  \
		 $REFDIR/Acer/GCF_001442555.1_ACSNU-2.0_genomic.fna \
		 $REFDIR/Ador/GCF_000469605.1_Apis_dorsata_1.3_genomic.fna \
	         $REFDIR/Aflo/GCF_000184785.3_Aflo_1.1_genomic.fna  
##\
##		 $REFDIR/Alab/GCA_014066325.1_ASM1406632v1_genomic.fna

