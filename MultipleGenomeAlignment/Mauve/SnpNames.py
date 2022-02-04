from __future__ import division
import pandas as pd
import os
from collections import defaultdict
from math import ceil
import sys

# This is a script to prepare the input file for est-sfs for the inferrence of the ancestral allele
# The target specie is Apis mellifera and outgroups are Apis cerana, Apis dorsata and Apis florea
# The input files are the output from the multiple whole genome alignment with Mauve and a vcf file of Apis mellifera samples
# A auxiliary input file is the chromosome-name conversion file holding different names for Apis mellifera chromosomes
# 1) to convert .xmfa to .snps file with Mauve and its class SnpExporter
# 2) extract the INFO (AF and AC) for the SNPs in the .snps file from the vcf file
# 3) create a dictionary to hold the est-sfs coded SNPs for Apis mellifera samples from the vcf and outfroups from .snps

noCycle= int(sys.argv[1])

# Set home directory
homeDir = "/home/v1jobste/jobsteter/Honeybees/MultipleGenomeAlignment/"
#upDir = "/home/x/HoneybeeGeno/"

snpsFile = homeDir + "multiBees_long.snps"
chrConvFile = homeDir + "ChromosomeConversion.csv"
vcfInfo = homeDir + "VCFSnps_AmelChr.INFO"
#vcFile = upDir + "HarpurPNAS.vcf"

######### --- 1 ---###########
# Transform .xmfa to .snps
#os.system("java -Xmx20g -cp ~/bin/mauve_snapshot_2015-02-13/Mauve.jar org.gel.mauve.analysis.SnpExporter -f multiBees_long.xmfa -o multiBees_long.snps")

# Read in the chromosome-name conversion file
print("Reading in chromosome conversion data")
chrConv = pd.read_csv(chrConvFile)
chrConv.loc[:, "VCFName"] = "carnica." + chrConv.Name
# Create a dictionary of RefSeq:VCFName
#chrConvD = {rs:vcf for rs, vcf in zip(chrConv.RefSeq, chrConv.VCFName)}


######### --- 2 ---###########
# First extract the first three columns from .snps file (the snp pattern, Apis mellifera config and Apis mellifera position)
# Also remove lines with NULL for Apis mellifera from the .snps file
# print("Extracting columns from the .snps file and removing the null lines")
# os.system("cut -f1,2,3 " + snpsFile + " > MASnps.txt")
# os.system("""awk '{if (($2 == "null")) next; else print $0 }' MASnps.txt > MASnps_Amel.txt""")
# # Keep only the SNPs on Apis mellifera chromosomes
# os.system("""cut -f3 -d"," ChromosomeConversion.csv  > Chromosomes_Amel.txt""")
# os.system("grep -Fwf Chromosomes_Amel.txt MASnps_Amel.txt > MASnps_AmelChr.txt")
# # Sort the .snps file
# os.system("sort --version-sort -k2,3 MASnps_AmelChr.txt > MASnps_AmelChrSorted.txt")

# Set the name of the modified .snps file
print("Reading in the .snps file.")
# Read in the snps file
snps = pd.read_csv("MASnps_AmelChrSortedVcf.txt", sep=" ", header=None) #MASnps_AmelChrSorted.txt
snps.columns =["SNP_pattern", "seq1Contig", "seq1Pos", "Pos"]
nRowCycle = ceil(len(snps) / 200)
snps = snps[nRowCycle * noCycle: nRowCycle * (noCycle+1)]
pd.DataFrame({"ID": snps.Pos}).to_csv("SNPNames" + str(noCycle) + ".csv", index=None, header=None)


