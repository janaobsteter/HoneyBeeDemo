# %%
import allel
import sys
import pandas as pd
import numpy as np

# %%
popType=sys.argv[1]
pop=sys.argv[2]



for chromosome in range(1, 17):
    vcf = allel.read_vcf("/home/jobste/Data/Remapped/LowRecRegionsRemoved/Chr" + str(chromosome) + "_masked.vcf")

    # %%
    meta = pd.read_csv("/home/jobste/Data/TsinferMetaData_Vcf.txt", index_col = "ID")

    # %%

    meta.loc[:, "Index"] = range(len(meta))
    s = list(vcf['samples'])
    meta = meta.reindex(s)
    meta = meta.set_index("Index")

    gt = allel.GenotypeArray(vcf['calldata/GT'])

    #ACS
    # Pop
    if popType == "pop":
        pops = set(meta.Pop)
        popDict = {pop: meta[meta.Pop == pop].index for pop in pops}
        acsPop = gt.count_alleles_subpops(popDict)

    # Subpop
    if popType == "subpop":
        pops = set(meta.SubPop)
        popDict = {subpop: meta[meta.SubPop == subpop].index for subpop in pops}
        acsPop = gt.count_alleles_subpops(popDict)

    if chromosome == 1:
        countArray = np.zeros(len(popDict[pop]) + 1, dtype=int, order='C')



    # %%
    popArray = np.zeros(len(popDict[pop]) + 1, dtype=int, order='C')
    popSFS = allel.sfs(np.array([x[1] for x in acsPop[pop]]))
    popArray[range(len(popSFS))] = popSFS

    countArray = countArray + popArray



pd.DataFrame({"SFS": countArray[1:-1]}).to_csv("Genome_pop" + str(pop) + "_unfold_count.sfs", index = None, header=None)


