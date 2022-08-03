import pandas as pd
import sys

"""This script read in the data computed with PrepareRelateInputFiles and edits them to be used by Relate"""
chr = sys.argv[1]

# Read in the .haps file, leave in only one of the two identical alleles at a locus (make it haploid)
haps = pd.read_table("AncestralChr" + str(chr) + ".haps", sep=" ")
nCol = haps.shape[1]
nRow = haps.shape[0]
selCols = [0,1,2,3,4] + [x for x in range(5, nCol-1) if x % 2 == 1]
selCols = [x for x in map(str, selCols)]
hapsSel = haps[selCols]
hapsSel = hapsSel.sort_values(by="2")
hapsSel.to_csv("AncestralChr" + str(chr) + "_haploid.haps", sep=" ", header=None, index = None)

# Compute missingness by SNP for the .sample table
snpCols = [x for x in map(str, [x for x in range(5, nCol-1) if x % 2 == 1])]
missing = hapsSel[snpCols].apply(lambda col: sum(col == "?") / nRow, axis = 0)

# Read in .samples and add the missingness
samples = pd.read_table("Chr" + str(chr) + ".sample", sep = " ")
samples['missing'] = [0] + list(missing)
samples['ID_2'] = "NA"
samples.to_csv("Chr" + str(chr) + "_missing.sample", index = None, sep= " ")

meta = pd.read_csv("SampleMetaDataFull.csv")
samplePop = pd.DataFrame({ "sample": samples.ID_1,
                           "population":[list(meta.Type[meta.Name == x])[0] if x in list(meta.Name) else "hybrid" for x in samples.ID_1],
                           "group": [list(meta.Lineage[meta.Name == x])[0] if x in list(meta.Name) else "H" for x in samples.ID_1],
                            "country": [list(meta.Country[meta.Name == x])[0] if x in list(meta.Name) else "hybrid" for x in samples.ID_1]})
samplePop = samplePop.drop(0)
samplePop.loc[:, "sex"] = 1
samplePop.to_csv("SamplePop.csv", index = None, sep=" ")
