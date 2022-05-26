setwd("/home/x/EddieDir/Honeybees/tsinfer/")
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(Rmisc)
# A script to plot the output of the trees

#FST
fst <- read.csv("FstCombined.csv")
fst$GroupPair <- paste0(fst$Group1, "_", fst$Group2)
fst$Chromosome <- as.factor(fst$Chromosome)
fstPlot <- ggplot(data = fst, aes(x = GroupPair, y = Fst, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T)
fstSum <- read.csv("FstSummary.csv")
fstSum$GroupPair <- paste0(fstSum$Group1, "_", fstSum$Group2)
ggplot(data = fstSum, aes(x = GroupPair, y = mean)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean Fst")


#DIvergence
divergence <- read.csv("DivergenceCombined.csv")
divergence$GroupPair <- paste0(divergence$Group1, "_", divergence$Group2)
divergence$Chromosome <- as.factor(divergence$Chromosome)
divergencePlot <- ggplot(data = divergence, aes(x = GroupPair, y = Divergence, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T)
divergenceSum <- read.csv("DivergenceSummary.csv")
divergenceSum$GroupPair <- paste0(divergenceSum$Group1, "_", divergenceSum$Group2)
ggplot(data = divergenceSum, aes(x = GroupPair, y = mean)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean divergence")

#Genetic relatendess
grel <- read.csv("GrelCombined.csv")
grel$GroupPair <- paste0(grel$Group1, "_", grel$Group2)
grel$Chromosome <- as.factor(grel$Chromosome)
grelPlot <- ggplot(data = grel, aes(x = GroupPair, y = GeneticRel, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T)
grelSum <- read.csv("GrelSummary.csv")
grelSum$GroupPair <- paste0(grelSum$Group1, "_", grelSum$Group2)
ggplot(data = grelSum, aes(x = GroupPair, y = mean)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean genetic relatedness")


fstL <- fst %>% pivot_longer(Fst, values_to = "Value", names_to = "Metric")
divergenceL <- divergence %>% pivot_longer(Divergence, values_to = "Value", names_to = "Metric")
grelL <- grel %>% pivot_longer(GeneticRel, values_to = "Value", names_to = "Metric")
intraGroup <- rbind(fstL, divergenceL, grelL)

ggplot(data = intraGroup, aes(x = GroupPair, y = Value, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T) + 
  facet_grid(rows = vars(Metric), scales = "free_y")



#TajimasD
tajimaD <- read.csv("TajimasDCombined.csv")
tajimaD$Chromosome <- as.factor(tajimaD$Chromosome)
ggplot(data = tajimaD, aes(x = Group, y = TajimasD, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T)
tajimaDSum <- read.csv("TajimasDSummary.csv")
ggplot(data = tajimaDSum, aes(x = Group, y = mean, fill = Group)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean TajimasD") + 
  scale_fill_viridis(discrete=T)

#Diversity
diversity <- read.csv("DiversityCombined.csv")
diversity$Chromosome <- as.factor(diversity$Chromosome)
ggplot(data = diversity, aes(x = Group, y = Diversity, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T)
diversitySum <- read.csv("DiversitySummary.csv")
ggplot(data = diversitySum, aes(x = Group, y = mean)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean diversity")


tajimaDL <- tajimaD %>% pivot_longer(TajimasD, values_to = "Value", names_to = "Metric")
diversityL <- diversity %>% pivot_longer(Diversity, values_to = "Value", names_to = "Metric")
interGroup <- rbind(tajimaDL, diversityL)
ggplot(data = interGroup, aes(x = Group, y = Value, fill = Chromosome)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 16) + 
  scale_fill_viridis(discrete=T) + 
  facet_grid(rows = vars(Metric), scales = "free_y")

#GNN of H lineage
gnn <- read.csv("GnnCombined.csv")
meta <- read.csv("Meta.csv")
gnn <- merge(gnn, meta, by.x = "Names", by.y = "ID", all.x = TRUE)
# Check whether all the ones with missing Count5ry ared REUd --> these are from France
unique(gnn$Names[is.na(gnn$Country)])
gnn$Country[is.na(gnn$Country)] <- "FRA"
sum(is.na(gnn$Country))

table(gnn$Type)
gnn$Chromosome <- as.factor(gnn$Chromosome)
gnnL <- gnn %>%  pivot_longer(c(A, C, M, O), names_to = "Lineage", values_to = "GNN")

ggplot(data = gnnL, aes(x = Country, y = GNN, fill = Lineage)) + 
  geom_boxplot() + 
  theme_bw(base_size = 16)

#MeanDesc of H lineage
meanDesc <- read.csv("MeanDescCombined.csv")
meanDesc <- merge(meanDesc, meta, by.x = "NodeName", by.y = "ID", all.x = TRUE)



#AFS
afs <- read.csv("AFSCombined.csv")
ggplot(data = afs, aes(x = X0)) + geom_histogram(bins = 100) + 
  facet_grid(cols = vars(Lineage))


################################################3
#Weighted by chromosomes

# Fst plot weighted by chromosome
fst <- read.csv("~/EddieDir/Honeybees/tsinfer/Statistics/FstSubspecieChrPer.csv")
fst$GroupPair1 <- paste0(fst$Subspecie1, "_", fst$Subspecie2)
fst$GroupPair2 <- paste0(fst$Subspecie2, "_", fst$Subspecie1)


subspecies = c("capensis", "scutellata", "unicolor",
               "mellifera", "carnica", "ligustica", "caucasica", "hybrid")
fstHeatMap <- expand.grid(subspecies, subspecies)
fstHeatMap$GroupPair <- paste0(fstHeatMap$Var1, "_", fstHeatMap$Var2)
fstHeatMap <- merge(fstHeatMap, fst[, c("GroupPair1", "Fst")], by.x = "GroupPair", by.y = "GroupPair1", all.x=T)
fstHeatMap <- merge(fstHeatMap, fst[, c("GroupPair2", "Fst")], by.x = "GroupPair", by.y = "GroupPair2", all.x=T)
fstHeatMap$Fst <- ifelse(!is.na(fstHeatMap$Fst.x), fstHeatMap$Fst.x, fstHeatMap$Fst.y)


triAngNames <- c()
for (no in 1:length(subspecies)) {
  for (s2 in subspecies[(no+1):8]) {
    triAngNames <- c(triAngNames, paste0(subspecies[no], "_", s2))
  }
}

fstNonHybrid <- fstHeatMap[(fstHeatMap$Var1 != "hybrid") & (fstHeatMap$Var2 != "hybrid"),]
subspecieFst <- ggplot(fstNonHybrid[fstNonHybrid$GroupPair %in% triAngNames,], aes(x = Var2, y = Var1)) +
                  geom_tile(aes(fill = Fst)) + scale_fill_viridis(na.value="black") + 
                  theme_bw(base_size = 18) + 
                  theme(panel.grid = element_blank(), axis.title = element_blank(),
                        axis.text.x = element_text(angle = 90)) + theme(legend.position = "none")



# Fst plot weighted by chromosome
fst <- read.csv("~/EddieDir/Honeybees/tsinfer/FstLineageChrPer.csv")
fst$GroupPair1 <- paste0(fst$Subspecie1, "_", fst$Subspecie2)
fst$GroupPair2 <- paste0(fst$Subspecie2, "_", fst$Subspecie1)


lineage = c("A", "M", "C", "O", "H")
fstHeatMapL <- expand.grid(lineage, lineage)
fstHeatMapL$GroupPair <- paste0(fstHeatMapL$Var1, "_", fstHeatMapL$Var2)
fstHeatMapL <- merge(fstHeatMapL, fst[, c("GroupPair1", "Fst")], by.x = "GroupPair", by.y = "GroupPair1", all.x=T)
fstHeatMapL <- merge(fstHeatMapL, fst[, c("GroupPair2", "Fst")], by.x = "GroupPair", by.y = "GroupPair2", all.x=T)
fstHeatMapL$Fst <- ifelse(!is.na(fstHeatMapL$Fst.x), fstHeatMapL$Fst.x, fstHeatMapL$Fst.y)


triAngNamesL <- c()
for (no in 1:length(lineage)) {
  for (s2 in lineage[(no+1):5]) {
    triAngNamesL <- c(triAngNamesL, paste0(lineage[no], "_", s2))
  }
}

fstNonHybridL <- fstHeatMapL[(fstHeatMapL$Var1 != "H") & (fstHeatMapL$Var2 != "H"),]
lineageFst <- ggplot(fstNonHybridL[fstNonHybridL$GroupPair %in% triAngNamesL,], aes(x = Var2, y = Var1)) +
                  geom_tile(aes(fill = Fst)) + scale_fill_viridis(na.value="black") + 
                  theme_bw(base_size = 18) + 
                  theme(panel.grid = element_blank(), axis.title = element_blank(),
                        axis.text.x = element_text(angle = 90)) + theme(legend.position = "bottom")


fstNonHybrid$Type <- "subspecie"
fstNonHybridL$Type <- "lineage"

library(Rmisc)
multiplot(subspecieFst, lineageFst, cols=2)




# Fst plot weighted by chromosome
tajimaDL <- read.csv("~/EddieDir/Honeybees/tsinfer/TajimaDLineageChrPer.csv")
tajimaDS <- read.csv("~/EddieDir/Honeybees/tsinfer/TajimaDSubspecieChrPer.csv")

lineage = c("A", "M", "C", "O", "H")

lineageTajimaD <- ggplot(tajimaDL[tajimaDL$Lineage != "H",], aes(x = Lineage, y = TajimaD, fill = Lineage)) +
  geom_col() + scale_fill_viridis("", na.value="black", discrete = T) + 
  theme_bw(base_size = 18) + 
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)) + 
  theme(legend.position = "right")

subspecies = c("capensis", "scutellata", "unicolor",
               "mellifera", "carnica", "ligustica", "caucasica", "hybrid")
tajimaDS$Subspecie <- factor(tajimaDS$Subspecie, levels =  subspecies)
subspecieeTajimaD <- ggplot(tajimaDS[tajimaDS$Subspecie != "hybrid",], aes(x = Subspecie, y = TajimaD, fill = Subspecie)) +
  geom_col() + scale_fill_viridis("", na.value="black", discrete = T) + 
  theme_bw(base_size = 18) + 
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)) + 
  theme(legend.position = "left") + 
  guides(fill = guide_legend(nrow=7))



multiplot(subspecieeTajimaD, lineageTajimaD, cols=2)


#Gnn
samples <- read.csv("SampleMetaData.csv")
meta <- read.csv("~/JANA/HoneybeeGeno/tsinfer/Meta.csv")
length(intersect(samples$Name, meta$ID))
samples$Name[!samples$Name %in% intersect(samples$Name, meta$ID)]
samples <- merge(samples, meta, by.x = "Name", by.y = "ID", all.x = T)
nrow(samples)
write.csv(samples, "~/JANA/HoneybeeGeno/tsinfer/SampleMetaDataFull.csv", row.names=F, quote=F)

gnn <- read.csv("GNN_hybridDronesCombined.csv")
gnn <- merge(gnn, samples, by =  "Name", all.x=T)

gnnA <- gnn %>%  group_by(Subspecie, Lineage, Name, Country) %>%  dplyr::summarise(A = mean(A),
                                                                          C = mean(C),
                                                                          M = mean(M),
                                                                          O = mean(O))
gnnAL <- gnnA %>% pivot_longer(c(A, C, M, O))
gnnAL[is.na(gnnAL$Country), "Country"] <- "REU"

ggplot(gnnAL, aes(x = Country, y = name , fill = value)) + geom_tile() + 
  scale_fill_viridis("GNN") + theme_bw(base_size = 18) + 
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 90)) + 
  theme(legend.position = "left") + 
  guides(fill = guide_legend(nrow=7))
  