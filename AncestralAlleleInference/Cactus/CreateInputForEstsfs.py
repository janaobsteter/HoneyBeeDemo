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
homeDir = "/home/v1jobste/jobsteter/Honeybees/MultipleGenomeAlignment/EstSfs/"
#upDir = "/home/x/HoneybeeGeno/"

# Set input file names
snpsFile = homeDir + "multiBees_long.snps"
chrConvFile = homeDir + "ChromosomeConversion.csv"
vcfInfo = homeDir + "VCFSnps_AmelChr1.INFO"
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

if os.path.exists(homeDir + "EstSfs_Dict.csv"):
    os.remove(homeDir + "EstSfs_Dict.csv")

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
# # Replace RefSeq chromosome name by the name that appears in the vcf file (could be different to what is set now!)
# for rs, vcf in zip(chrConv.RefSeq, chrConv.VCFName):
#     os.system("sed -i s/" + rs + "/" + vcf + "/g MASnps_AmelChrSorted.txt")
#
#
# # Extract only the VCF SNPs from the alignment file
# # Create a Position columns
# os.system("""awk '{$4=$2"_"$3} 1' MASnps_AmelChrSorted.txt > tmp && mv tmp MASnps_AmelChrSorted.txt""")
# os.system("""awk '{print $1"_"$2}' """ + vcfInfo + """ > VcfPos.txt """)
# os.system("grep -Fwf VcfPos.txt MASnps_AmelChrSorted.txt > MASnps_AmelChrSortedVcf.txt")
#

# Set the name of the modified .snps file
print("Reading in the .snps file.")
# Read in the snps file
snps = pd.read_csv("MASnps_AmelChrSortedVcf.txt", sep=" ", header=None) #MASnps_AmelChrSorted.txt
snps.columns =["SNP_pattern", "seq1Contig", "seq1Pos", "Pos"]
nRowCycle = ceil(len(snps) / 200)
snps = snps[nRowCycle * noCycle: nRowCycle * (noCycle+1)]


# Prepare a list of SNPs to extract from the vcf file
print("Extracting SNP info from the vcf file")
# Cut only the positions from the file
os.system("cut -f2,3 MASnps_AmelChrSortedVcf.txt > MASnps_AmelPosChr.txt")
# Use vcftools to extract the alternate allele count (AC) and allele frequency (AF) for the given set of SNPs from the vcf file
#os.system("vcftools --vcf " + vcFile + " --positions MASnps_AmelPosChr.txt --get-INFO AC --get-INFO AF --out VCFSnps_AmelChr")
# Read in the vcf info
snpsInfo = pd.read_csv(vcfInfo, sep="\t")
snpsInfo.loc[:, "Pos"] = snpsInfo['CHROM'].str.cat(snpsInfo['POS'].apply(str),sep="_")



######### --- 3 ---###########
# Create representation of the four SNPs for the est-sfs
snpDict = {"A": (1,0,0,0), "C": (0,1,0,0), "G":(0,0,1,0), "T":(0,0,0,1), "-":(0,0,0,0), "N":(0,0,0,0)}

# Create a list of species following the order on the .snps file
print("Creating a dictionary with est-sfs coded snps")
species = ["Amel", "Acer", "Ador", "Aflo"]

# Create a dictionary holding est-sfs coded info for each snp
specieValues = defaultdict()


# Create a loop to write the snp info to dictionary
print("Creating an est-sfs dictionary of alelle counts for each specie.")
for snpPos, snpPattern in iter(zip(snps.Pos, snps.SNP_pattern)):
    # Check whether the SNP is also in the vcf file
    if snpPos in list(snpsInfo.Pos):
        snpSpecDict = {}
        snpLine = snpsInfo.query("Pos == @snpPos")
        amelAllele = list(snpPattern)[0].upper()
        # Extract count and frequency for the alternate allele and compute the total allele count
        # This is now an iterator!!! Iterate only once!!!! Even with a sum!
        altFreq = {altAl:float(altValue) for altAl, altValue in zip(list(snpLine.ALT)[0].split(","), list(snpLine.AF)[0].split(","))}
        altCount = {altAl:int(altValue) for altAl, altValue in zip(list(snpLine.ALT)[0].split(","), list(snpLine.AC)[0].split(","))}
        totalCount = max([round(ac / af) for ac, af in zip(altCount.values(), altFreq.values()) if af != 0])
        altSum = sum(altCount.values())
        refSum = totalCount - altSum
        refAllele = list(snpLine.REF)[0]

        if len(refAllele) == 1:
            # Create an est-sfs dict for the reference allele
            vcfDict = []
            refCountDict = tuple((tuple(x*refSum for x in snpDict.get(refAllel)) for refAllel in  refAllele if refAllele in snpDict.keys()))
            if refCountDict:
                vcfDict.append(refCountDict[0])
            # Create an est-sfs dict for the alternate alleles
            altCountDict = list((tuple(x*altAlleleCount for x in snpDict.get(altAllele))
                                 for altAllele, altAlleleCount in zip(altCount.keys(), altCount.values())
                                 if altAllele in snpDict.keys()))
            vcfDict = vcfDict + altCountDict

            # Zip and sum all the tuples within the vcfDict list
            aMelCount = tuple([sum(x) for x in zip(*vcfDict)])

            # Add est-sfs code for the SNP into the dictionary
            # Create a dict for the three outgroups from the mauve output
            snpCountDict = [snpDict[x.upper()] for x in snpPattern[1:]]
            # Create a dict for the Apis mellifera from the mauve output
            aMelSnpCountDict = snpDict[snpPattern[0].upper()]
            # Sum the mauve Apis mellifera and vcf Apis mellifera alleles
            aMelCount = tuple([sum(x) for x in zip(aMelCount, aMelSnpCountDict)])
            # Add the amel summed counts as the first element of the list
            snpCountDict.insert(0, aMelCount)
            # Assign the list of tuples for each Apis species to the dictionary under the key = SNPName
            specieValues[snpPos] = snpCountDict
            snpSpecDict[snpPos] = snpCountDict
            # Write the dictionary of tupples into a dataframe - for one SNP at a time! mode = "a" --> append
            #pd.DataFrame.from_dict(snpSpecDict, orient='index').to_csv("EstSfs_DictEachSnp.csv", index=None, header=None, sep="\t", mode="a")
            #pd.DataFrame({"ID": [snpPos]}).to_csv("SnpNames.csv", index=None, header=None, mode="a")
			

# Write the dictionary of tupples into a dataframe
pd.DataFrame.from_dict(specieValues, orient='index').to_csv("EstSfs_Dict" + str(noCycle) + ".csv", header=None, sep="\t")

# Remove the spaces in the output file
# os.system("""sed -i "s/ //g" EstSfs_Dict.csv""")
# # Set the correct separators for the file
# os.system("""awk -F "\t" '{print $1"\t"$2" "$3" "$4}' EstSfs_Dict.csv > tmp && mv tmp EstSfs_Dict.csv""")
# # Remove the parenthesis from the file
# os.system("""sed -i "s/(//g" EstSfs_Dict.csv""")
# os.system("""sed -i "s/)//g" EstSfs_Dict.csv""")


