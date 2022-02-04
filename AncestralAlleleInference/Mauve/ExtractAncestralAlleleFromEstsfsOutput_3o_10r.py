#!/usr/bin/env python3

#Script to transform est-sfs output into ancestral allele states file

#Usage: ./script.py est-sfs-pvalues.txt ancestral_alleles.txt

import os
import pandas as pd
import sys
from math import isnan

cycle = sys.argv[1] #str(input("Name of est-sfs output file: "))


print("Reading est-sfs output")
est = pd.read_csv("EstsfsOutput_3o_10r/output-" + str(cycle) + "-pvalues.txt", sep = " ", header = None, skiprows = 8)
est.columns = ["LineNumber", "ConfigIndex", "pMajAnc", "pAA", "pAC", "pAG", "pAT", "pCA", "pCC", "pCG", "pCT", "pGA", "pGC", "pGG", "pGT", "pTA", "pTC", "pTG", "pTT", "last"]
namefile = pd.read_csv("EstSfsNames" + str(cycle) + ".csv", header=None)
namefile.columns = ["SnpPos"]

alleles= ["A", "C", "G", "T"]
ancestral_alleles = []
for index, site in est.iterrows():
    if not all([isnan(x) for x in site[3:19]]):
        #print("Site_0", site['LineNumber'])
        #print("Site's probability of MAJ == ANC")
        #print(site['pMajAnc'])
        #print("Site_2", site['pMajAnc'])
        probA = site['pAA'] + site['pAC'] + site['pAG'] + site['pAT']
        probC = site['pCA'] + site['pCC'] + site['pCG'] + site['pCT']
        probG = site['pGA'] + site['pGC'] + site['pGG'] + site['pGT']
        probT = site['pTA'] + site['pTC'] + site['pTG'] + site['pTT']
        #print("Prob A: ", probA)
        #print("Prob C: ", probC)
        #print("Prob G: ", probG)
        #print("Prob T: ", probT)
        probs = [probA, probC, probG, probT]
        larger = max(probs)
        larger_index = probs.index(larger)
        #print(larger)
        #print(larger_index)
        #print(alleles[larger_index])
        #print(larger, " is the ancestra allele state")
        ancestral_alleles.append(alleles[larger_index])
        #print("End of site\n")
    else:
        ancestral_alleles.append("")
#print(ancestral_alleles)


namefile.loc[:, "ancAl"] = ancestral_alleles
namefile.ancAl.value_counts()
namefile.to_csv("EstsfsOutput_3o_10r/AncestralAllele" + str(cycle) + "_3o_10r.csv", index=None)
