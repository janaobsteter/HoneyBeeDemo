import pandas as pd
import os
import numpy as np
import sys

workingDir="/home/v1jobste/jobsteter/Honeybees/Relate/"
vcfDir="/home/v1jobste/jobsteter/Honeybees/GenomicData/Wragg/Imputed/"
relateDir="/home/v1jobste/jobsteter/Honeybees/Relate/relate_v1.1.8_x86_64_dynamic/"
ancestralFile="/home/v1jobste/jobsteter/Honeybees/MultipleGenomeAlignment/Cactus/EstSfs/EstsfsOutput_3o/AncestralAlleles_HiFiVcf.csv"
ancestralDF = pd.read_csv(ancestralFile, header=None)



os.chdir(workingDir)
chr = sys.argv[1]

os.system("bcftools convert " + vcfDir + "Chr" + str(chr) + "_phased.vcf.gz --hapsample " + workingDir + "Chr" + str(chr))
os.system("gunzip " + workingDir + "Chr" + str(chr) + ".hap.gz")
os.system('sed -i "s/*//g" Chr' + str(chr) + '.hap')

aa = pd.DataFrame()
counter = 1
for chunk in pd.read_table("Chr" + str(chr) + ".hap", sep=" ", header=None, dtype=str, chunksize=1000):
    print(counter)
    counter += 1
    snpPos = ["_".join(x.split("_")[0].split(":")) for x in list(chunk[1])]
    refAl = list(chunk[3])
    chunk['snpPos'] = snpPos
    snpAAInfo = list(set(snpPos).intersection(set(ancestralDF[0])))

    snpAAInfoDF = chunk[chunk['snpPos'].isin(snpAAInfo)]

    # Which match/dont match the ancestral
    refMatchAnc = [snp for snp, ref in zip(snpAAInfoDF.snpPos, snpAAInfoDF[3]) if ref == ancestralDF[1][ancestralDF[0] == snp].item()]
    refMatchAncDF = snpAAInfoDF[snpAAInfoDF['snpPos'].isin(refMatchAnc)]
    refNonMatchAnc = list(set(snpAAInfoDF.snpPos) - set(refMatchAnc))
    refNonMatchAncDF = snpAAInfoDF[snpAAInfoDF['snpPos'].isin(refNonMatchAnc)]

    aa = aa.append(refMatchAncDF)
    aa = aa.append(refNonMatchAncDF.replace(to_replace = {'0': '1', '1': '0'}))


aa.to_csv("AncestralChr" + str(chr) + ".haps", sep=" ", index = None)
    #snpAANonInfo = list(set(snpPos) - set(snpAAInfo))
    #snpAANonInfoDF = chunk[chunk['snpPos'].isin(snpAANonInfo)]


    # for index, row in chunk.iterrows():
    #     print(index)
    #     name = row[1].split("_")[0].split(":")
    #     name = "_".join(name)
    #     if name in list(ancestralDF[0]):
    #         AA = ancestralDF[1][ancestralDF[0] == name].item()
    #         if row[3] == AA:
    #             continue
    #         else:
    #             newRow = row.replace(to_replace = {'0': '1', '1': '0'})
    #             chunk.loc[[index]] = newRow
