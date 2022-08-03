import pandas as pd

genMapsDir = "/home/v1jobste/jobsteter/Honeybees/GenomicData/Wragg/Imputed/GenMaps/"
genMapsDir = "/home/x/EddieDir/Honeybees/GenomicData/Wragg/Imputed/GenMaps/"

_recombination_rate = {
    "CM009931.2": 23.9e-8,
    "CM009932.2": 24.6e-8,
    "CM009933.2": 24.1e-8,
    "CM009934.2": 27.6e-8,
    "CM009935.2": 21.4e-8,
    "CM009936.2": 21.2e-8,
    "CM009937.2": 23.4e-8,
    "CM009938.2": 20.9e-8,
    "CM009939.2": 24.6e-8,
    "CM009940.2": 24.8e-8,
    "CM009941.2": 20.3e-8,
    "CM009942.2": 21.2e-8,
    "CM009943.2": 23.4e-8,
    "CM009944.2": 24.6e-8,
    "CM009945.2": 22.1e-8,
    "CM009946.2": 22.8e-8,
    "CM009947.2": 0,
}

rec = [x * 1e8 for x in list(_recombination_rate.values())]

#chr = 16
for chr in range(1, 17):
    genMap = pd.read_csv(genMapsDir + "Chr" + str(chr) + ".gmap", header = None, sep="\t")
    genMap.columns = ["Pos", "Chr", "GenPos"]
    relateGenMap = pd.DataFrame({"pod": genMap.Pos, "COMBINED_rate": rec[chr-1], "Genetic_Map": genMap.GenPos})
    relateGenMap.to_csv("RelateGenMaps/Chr" + str(chr) + ".gmap", sep=" ", index = None)
    # haps = pd.read_table("AncestralChr" + str(chr) + ".haps", sep=" ")
    # nCol = haps.shape[1]
    # nRow = haps.shape[0]
    # selCols = [0,1,2,3,4] + [x for x in range(5, nCol-1) if x % 2 == 1]
    # selCols = [x for x in map(str, selCols)]
    # hapsSel = haps[selCols]
    # hapsSel = hapsSel.sort_values(by="2")
    #
    # snps = haps[['2']]
    # snps = snps.assign(RecRate = [22.8] * len(snps))
    # snps.loc[:, "GenDis"] = ""
    # snps.loc[0, "GenDis"] = 0.110058
    #
    # for row in range(0, len(snps)-1):
    #     snps.loc[row+1, "GenDis"] = (22.8 / 1e6) * (list(snps.iloc[[row+1]]['2'])[0] - list(snps.iloc[[row]]['2'])[0]) + list(snps.iloc[[row]].GenDis)[0]
    #
    # snps.columns = ["pos" ,"COMBINED_rate" ,"Genetic_Map"]
    # snps.to_csv("GeneticMap_chr16.txt", sep=" ", index = None)

