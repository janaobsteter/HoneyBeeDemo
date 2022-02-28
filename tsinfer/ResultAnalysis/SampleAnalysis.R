##############################################
#-- This is a script to compare the honeybee genome data received from Brock Harpur and David Wragg
##############################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


# Harpur includes 39 diploid worker data
harpur <- read.csv("~/HoneybeeGeno/PopSamples.txt", header=F)
# Wragg include 734 haploid drone data
wragg <- read.csv("~/HoneybeeGeno/DavidWragg_DroneGeno/Runs.csv")

# Read in ancestry percentage from the DWragg samples
anc <- read.csv("~/HoneybeeGeno/DavidWragg_DroneGeno/Wragg_ancestryPercentage.csv", skip=3)
# Read in the sample meta-data
meta <- read.csv("~/HoneybeeGeno/DavidWragg_DroneGeno/Meta.csv")
table(meta$Type)
# There is 558 individuals in the anc file and 744 in the meta file
# Combine the ancestry percentage data and the meta data on the specie
metaAnc <- merge(anc, meta[,c("ID", "Type")], by.x = "Sample", by.y = "ID", all.y = T)
nrow(metaAnc)
head(metaAnc)


# First check how many of the drones, that were not in the ancestry analysis, has a "Type"
sum(is.na(metaAnc$Cluster1))
table(metaAnc$Type[is.na(metaAnc$Cluster1)])
nrow(metaAnc) - 173 # --> 571: this is the number of samples that have either cluster analysis or are not marker as hybrid

# Produce a longer version of the table
metaAncL <- metaAnc %>% pivot_longer(cols = c(paste0("Cluster", c(1, 2, 3, 4, 5))), names_to = "Cluster", values_to = "per") %>% 
                        pivot_longer(cols = c(paste0("Cluster", c(1.1, 2.1, 3.1, 4.1, 5.1))), names_to = "ClusterA", values_to = "perA")

# Look at the cluster percentage estimated with Amel Reference per type
metaAncLS <- metaAncL %>% group_by(Type, Cluster) %>% dplyr::summarise(MeanPer = mean(per, na.rm=T),
                                                                MeanPerA = mean(perA, na.rm=T))
# Vizualize
metaAncLS$Cluster <- as.factor(metaAncLS$Cluster)
# ggplot(data = metaAncLS, aes(x = Type, y = MeanPer, fill = Cluster)) + 
#   geom_col(position = "dodge") + 
#   scale_fill_viridis(discrete = T)
  
# # Look at the cluster percentage estimated with ancestral allele
# Vizualize
tapply(metaAncLS$MeanPer, metaAncLS$Cluster, FUN= function(x) mean(x, na.rm=T))
tapply(metaAncLS$MeanPerA, metaAncLS$Cluster, FUN= function(x) mean(x, na.rm=T))
# ggplot(data = metaAncLS, aes(x = Type, y = MeanPerA, fill = Cluster)) + 
#   geom_col(position = "dodge") + 
#   scale_fill_viridis(discrete = T)

# A function to return the cluster value
returnValueCluster <- function(x) {
  if (is.na(x[1])) {
    return(NA) }
      else {
       return(colnames(metaAnc)[which(x == max(x)) +1]) 
      }
}

# Asign the individual to the cluster with maximum probability
# Compute the maximum cluster value
metaAnc$MAxValue <- apply(metaAnc[,2:6], MARGIN = 1, FUN= function(x) max(x))
summary(metaAnc$MAxValue)
# qplot(metaAnc$MAxValue, geom="histogram")
# Check how many samples has a cluster with probability > 0.9
sum(metaAnc$MAxValue > 0.9, na.rm=T)
# Asign the most likely cluster
metaAnc$ClusterMax <- apply(metaAnc[,2:6], MARGIN = 1, FUN= function(x) returnValueCluster(x))
# Check the cluster
table(metaAnc$ClusterMax)
# Apply the more stringent criteria
#metaAnc$ClusterMax[metaAnc$MAxValue < 0.9] <- ""
table(metaAnc$ClusterMax)
# Check concordance of the cluster and the type
table(metaAnc$ClusterMax, metaAnc$Type)


# Cluster with ancestral allele -------------------------------------------

# Do the same for the Cluster1.1, 2.1 ...
# For the clusters with the internal ancestral allele
metaAnc$MAxValueA <- apply(metaAnc[,7:11], MARGIN = 1, FUN= function(x) max(x))
# qplot(metaAnc$MAxValueA, geom="histogram")
sum(metaAnc$MAxValueA > 0.9)
nrow(metaAnc)
metaAnc$ClusterMaxA <- apply(metaAnc[,7:11], MARGIN = 1, FUN= function(x) returnValueCluster(x))
metaAnc$ClusterMaxA[metaAnc$MAxValueA < 0.9] <- ""
table(metaAnc$ClusterMaxA)
table(metaAnc$ClusterMaxA, metaAnc$Type)


# Asign the lineages ------------------------------------------------------


# Assign the lineages
# Asign the lineages accordind to type
metaAnc$Lineage <- ifelse(metaAnc$Type %in% c("carnica", "ligustica"), "C", 
                          ifelse(metaAnc$Type %in% c("scutellata", "capensis", "unicolor"), "A", 
                                 ifelse(metaAnc$Type %in% c("mellifera"), "M", 
                                        ifelse(metaAnc$Type == "caucasica", "O", NA))))
# Check
table(metaAnc$Lineage)
# Asign the lineages according to cluster analysis
metaAnc$ClusterLineage <- ifelse(metaAnc$ClusterMax == "Cluster1", "C",
                                 ifelse(metaAnc$ClusterMax == "Cluster2", "M", 
                                        ifelse(metaAnc$ClusterMax == "Cluster3", "A",
                                               ifelse(metaAnc$ClusterMax == "Cluster4", "", 
                                                      ifelse(metaAnc$ClusterMax == "Cluster5", "O", NA)))))
# Check
table(metaAnc$ClusterLineage)
metaAnc$ClusterLineage <- as.factor(metaAnc$ClusterLineage)
# Check concordance
table(metaAnc$ClusterLineage, metaAnc$Lineage)
sum(table(metaAnc$ClusterLineage, metaAnc$Lineage))

# Merge the two types of information
table(metaAnc$Lineage, metaAnc$ClusterLineage)
table(metaAnc$Type, metaAnc$ClusterLineage)
metaAnc$ClusterLineage <- as.factor(metaAnc$ClusterLineage)
metaAnc$Lineage <- as.factor(metaAnc$Lineage)
metaAnc$CombinedLineage <- ifelse((metaAnc$Type == "hybrid"), as.character(metaAnc$ClusterLineage), as.character(metaAnc$Lineage))
metaAnc$CombinedLineage[is.na(metaAnc$CombinedLineage)] <- ""
table(metaAnc$CombinedLineage)
sum(table(metaAnc$CombinedLineage))


#Write the table
write.csv(metaAnc[, c("Sample", "CombinedLineage")], "~/HoneybeeGeno/DavidWragg_DroneGeno/LineageInfomation.csv", quote=F, row.names=F)

