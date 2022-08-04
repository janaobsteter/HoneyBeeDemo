setwd("~/EddieDir/Relate/EffectivePopSize/Subspecies/")
library(ggplot2)
library(viridis)
library(dplyr)

popSizeTotal <- data.frame()

for (species in c("carnica", "ligustica", "mellifera", "capensis", "unicolor", "scutellata", "caucasica")) {
  for (chr in c(2,3,4,5,6,11,15)) {
    print(paste0(species, chr))
    argv <- commandArgs(trailingOnly = T)
    filename <- paste0("Chr", chr, "Ne_", species)
    years_per_gen <- 1
    
    #read in population size
    groups   <- as.matrix(read.table(paste(filename, ".coal", sep = ""), nrow = 1))
    t        <- years_per_gen*t(as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = 1, nrow = 1)))
    pop_size <- data.frame(time = numeric(0), pop_size = numeric(0), groups = numeric(0))
    num_pops <- round(sqrt(dim(read.table(paste(filename, ".coal", sep = ""), skip = 2))[1]))
    #num_pops <- (-1 + sqrt(1 + 8 * num_pops))/2
    
    for(p1 in 1:num_pops){
      for(p2 in 1:p1){
        c        <- as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = (p1-1) * num_pops + p2 + 1, nrow = 1))[-c(1:2)]
        str      <- rep(paste(groups[p1]," - ",groups[p2], sep = ""),length(c))
        pop_size <- rbind(pop_size, data.frame(time = t, pop_size = 0.5/c, groups = str))
      }
    }
    pop_size$time[which(pop_size$time > 1e7)] <- 1e7
    pop_size$chr <- chr
    pop_size$species <- species
    popSizeTotal <- rbind(popSizeTotal, pop_size)
  }
}

popSizeTotal$species <- factor(popSizeTotal$species, levels = c("capensis", "scutellata", "unicolor",
                                                                "caucasica", "mellifera", "carnica", "ligustica"))

write.table(popSizeTotal, "Ne_subspecies.csv", quote=F, row.names=F)

#plot
ggplot(popSizeTotal[popSizeTotal$chr == 15,]) +
  geom_step(aes(time, pop_size, color = species), lwd = 1.2) +
  scale_x_continuous(limits = c(1e3,1e7), trans="log10") + annotation_logticks(sides = "bl") +
  scale_y_continuous(trans ="log10") +
  ylab("population size") +
  xlab("years ago")+ 
  scale_colour_manual("", values = c("#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                     "#ad2793", "#310597"))


popSizeTotal <- popSizeTotal[popSizeTotal$pop_size != Inf,]
popSizeTotalA <- popSizeTotal %>% group_by(time, groups, species) %>% summarize(meanPopSize = mean(pop_size))
summary(popSizeTotalA$meanPopSize[popSizeTotalA$species == "unicolor"])
popSizeTotal[(popSizeTotal$species == "unicolor") & (popSizeTotal$time ==  13895),]
popSizeTotalA[(popSizeTotalA$species == "unicolor") & (popSizeTotalA$time ==  13895),]
ggplot(popSizeTotalA) +
  geom_step(aes(time, meanPopSize, color = species), lwd = 1.2) +
  scale_x_continuous(limits = c(1e3,1e7), trans="log10") + annotation_logticks(sides = "bl") +
  scale_y_continuous(trans ="log10", limits = c(1e2, 1e7)) +
  ylab("Population size") +
  xlab("Years ago") + 
  theme_bw(base_size=30) +
  scale_colour_manual("", values = c("#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                     "#ad2793", "#310597"))
