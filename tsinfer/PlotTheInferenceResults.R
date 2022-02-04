setwd("/home/x/EddieDir/")
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
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
ggplot(data = tajimaDSum, aes(x = Group, y = mean)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - std, max = mean + std)) + 
  theme_bw(base_size = 16) + 
  ylab("Mean TajimasD")

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




# Fst plot weighted by chromosome
fst <- read.csv("~/EddieDir/Honeybees/tsinfer/FstSubspecieChrPer.csv")
ggplot(fst, aes(Subspecie1, Subspecie2)) +
  geom_raster(aes(fill = Fst))

fst$GroupPair1 <- paste0(fst$Subspecie1, "_", fst$Subspecie2)
fst$GroupPair2 <- paste0(fst$Subspecie2, "_", fst$Subspecie1)
ggplot(fst, aes(x = GroupPair, y = Fst, fill = Subspecie1)) +
  geom_col()

fst <- fst[order(fst$Fst),]
fst$Fst <- as.numeric(fst$Fst)
fst$Subspecie1 <- factor(fst$Subspecie1, levels = c("hybrid", "scutellata", "capensis", "unicolor", "ligustica", 
                                                    "caucasica", "mellifera", "carnica"))
fst$Subspecie2 <- factor(fst$Subspecie2, levels =  c("hybrid", "scutellata", "capensis", "unicolor", "ligustica", 
                                                     "caucasica", "mellifera", "carnica"))
ggplot(fst[fst$Subspecie1 != "hybrid", ], aes(Subspecie1, Subspecie2)) +
  geom_tile(aes(fill = Fst))

subspecies = c("capensis", "scutellata", "unicolor",
               "mellifera", "carnica", "ligustica", "caucasica", "hybrid")
fstHeatMap <- expand.grid(subspecies, subspecies)
fstHeatMap$GroupPair <- paste0(fstHeatMap$Var1, "_", fstHeatMap$Var2)
fstHeatMap <- merge(fstHeatMap, fst[, c("GroupPair1", "Fst")], by.x = "GroupPair", by.y = "GroupPair1", all.x=T)
fstHeatMap <- merge(fstHeatMap, fst[, c("GroupPair2", "Fst")], by.x = "GroupPair", by.y = "GroupPair2", all.x=T)
fstHeatMap$Fst <- ifelse(!is.na(fstHeatMap$Fst.x), fstHeatMap$Fst.x, fstHeatMap$Fst.y)

ggplot(fstHeatMap[(fstHeatMap$Var1 != "hybrid") & (fstHeatMap$Var2 != "hybrid"),], aes(Var1, Var2)) +
  geom_tile(aes(fill = Fst)) + scale_fill_viridis()

ggplot(fstHeatMap[(fstHeatMap$Var1 == "hybrid") | (fstHeatMap$Var2 == "hybrid"),], aes(Var1, Var2)) +
  geom_tile(aes(fill = Fst)) + scale_colour_viridis()


ggplot(data = fst, aes(x = Subspecie1, y = Fst, fill = Subspecie1)) + 
  geom_col() + 
  facet_grid(rows = vars(Subspecie2)) + 
  scale_colour_manual(breaks = subspecies, values = viridis(7))


fst %>%
  mutate(Subspecie1 = fct_relevel(Subspecie1, 
                            subspecies)) %>%
  ggplot(aes(x = Subspecie1, y = Fst, fill = Subspecie1)) + 
  geom_col() + 
  facet_grid(rows = vars(Subspecie2)) + 
  scale_colour_manual(breaks = subspecies, values = viridis(7))

