import pandas as pd

fst = pd.DataFrame()
tajimaD = pd.DataFrame()
divergence = pd.DataFrame()
diversity = pd.DataFrame()
afs = pd.DataFrame()
f2 = pd.DataFrame()
f4 = pd.DataFrame()
gnn = pd.DataFrame()
grel = pd.DataFrame()
meanDesc = pd.DataFrame()

for chr in range(1, 17):
    print(chr)
    #fst
    tmp = pd.read_csv("FST" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    fst = fst.append(tmp)


    #TajimaD
    tmp = pd.read_csv("TajimasD" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    tajimaD = tajimaD.append(tmp)


    #Divergence
    tmp = pd.read_csv("Divergence" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    divergence = divergence.append(tmp)


    #Diversity
    tmp = pd.read_csv("Diversity" + str(chr) + ".csv")
    tmp.columns = ['Diversity', "Group"]
    tmp.loc[:, "Chromosome"] = chr
    diversity = diversity.append(tmp)


    #f2
    tmp = pd.read_csv("f2" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    f2 = f2.append(tmp)


    #GNN
    tmp = pd.read_csv("GNN_drones" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    gnn = gnn.append(tmp)


    #grel
    tmp = pd.read_csv("GeneticRelatedness" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    grel = grel.append(tmp)


    #MeanDesc
    tmp = pd.read_csv("MeanDescendants" + str(chr) + ".csv")
    tmp.loc[:, "Chromosome"] = chr
    meanDesc = meanDesc.append(tmp)

    #AFS
    for lineage in ["A", "C", "M", "O"]:
        tmp = pd.read_csv("AFS_Alineage" + str(chr) + ".csv", header=None)
        tmp.loc[:, "Lineage"] = lineage
        tmp.loc[:, "Chromosome"] = chr
        afs = afs.append(tmp)



fst.to_csv("FstCombined.csv", index = False)
fst.groupby(['Group1', 'Group2'])['Fst'].describe().to_csv("FstSummary.csv", index=True)
tajimaD.to_csv("TajimasDCombined.csv", index = False)
tajimaD.groupby(['Group'])['TajimasD'].describe().to_csv("TajimasDSummary.csv", index=True)
divergence.to_csv("DivergenceCombined.csv", index = False)
divergence.groupby(['Group1', 'Group2'])['Divergence'].describe().to_csv("DivergenceSummary.csv", index=True)
diversity.to_csv("DiversityCombined.csv", index = False)
diversity.groupby(['Group'])['Diversity'].describe().to_csv("DiversitySummary.csv", index=True)
f2.to_csv("F2Combined.csv", index = False)
f2.groupby(['Group1', 'Group2'])['f2'].describe().to_csv("F2Summary.csv", index=True)
gnn.to_csv("GnnCombined.csv", index = False)
gnn.groupby(['Names'])[('A', 'C', 'M', 'O')].mean().to_csv("GnnSummary.csv", index=True)
grel.to_csv("GrelCombined.csv", index = False)
grel.groupby(['Group1', 'Group2'])['GeneticRel'].describe().to_csv("GrelSummary.csv", index=True)
meanDesc.to_csv("MeanDescCombined.csv", index = False)
meanDesc.groupby(['NodeName'])[('A', 'C', 'M', 'O')].mean().to_csv("MeanDescSummary.csv", index=True)
afs.to_csv("AFSCombined.csv", index = False)