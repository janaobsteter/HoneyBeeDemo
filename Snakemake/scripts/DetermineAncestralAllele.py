#!/usr/bin/env python3

#Script to transform est-sfs output into ancestral allele states file
import os
import pandas as pd
import sys
from math import isnan

cycle = snakemake.wildcards[0]


print("Reading est-sfs output")
est = pd.read_csv("EstsfsOutput_3o/output-" + str(cycle) + "-pvalues.txt", sep = " ", header = None, skiprows = 8)
est.columns = ["LineNumber", "ConfigIndex", "pMajAnc", "pAA", "pAC", "pAG", "pAT", "pCA", "pCC", "pCG", "pCT", "pGA", "pGC", "pGG", "pGT", "pTA", "pTC", "pTG", "pTT", "last"]
namefile = pd.read_csv("EstSfsNames" + str(cycle) + ".csv", header=None)
namefile.columns = ["SnpPos"]

alleles= ["A", "C", "G", "T"]
ancestral_alleles = []
for index, site in est.iterrows():
    if not all([isnan(x) for x in site[3:19]]):
        probA = site['pAA'] + site['pAC'] + site['pAG'] + site['pAT']
        probC = site['pCA'] + site['pCC'] + site['pCG'] + site['pCT']
        probG = site['pGA'] + site['pGC'] + site['pGG'] + site['pGT']
        probT = site['pTA'] + site['pTC'] + site['pTG'] + site['pTT']
        probs = [probA, probC, probG, probT]
        larger = max(probs)
        larger_index = probs.index(larger)
        ancestral_alleles.append(alleles[larger_index])
    else:
	    ancestral_alleles.append("")

namefile.loc[:, "ancAl"] = ancestral_alleles
namefile.ancAl.value_counts()
namefile.to_csv("EstsfsOutput/AncestralAllele" + str(cycle) + "_3o.csv", index=None)
