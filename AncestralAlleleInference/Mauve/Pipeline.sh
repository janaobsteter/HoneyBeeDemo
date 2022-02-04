# This is the workflow to obtain the ancestral alleles for the Apis mellifera

1) Create the input for the est-sfs from the .vcf file and the alignment file. We used Mauve to alignt the genomes of Apis mellifera, Apis cerana, Apis dorsata, Apis florea with default parameters. After that we've used the java class SnpExporter to transform the .xmfa to .snps file.
To create the dictionary for the est-sfs I've created a python script that creates the dictionary in chunks - for now it creates 100 chunks!
python CreateInputForEstsfs.py

Output files are names Estsfs_Dict$N.csv

2) Edit the output files to extract the dictionary (they also contain the names of the snps!)
bash EditEstFile.sh
The output files are called EstSfs_Dict$NE.csv

3) Extract the names of the SNPs with 
bash GetNamesEstSfsDict.sh
The output files are called EstSfsNames$N.csv

4) Run the est-sfs with kimura setting file:
mkdir EstsfsOutput/
bash RunEstSfs_Loop.sh
This will submit 200 jobs for the est-sfs

5) After this you will have all the output in the EstsfsOutput/ directory.
You then need to run
