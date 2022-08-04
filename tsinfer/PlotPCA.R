library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)


pca <- read.csv("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/HiFiTrees/WeightedPCA_coordinated.csv")
colnames(pca) <- c("ID", paste0("PCA", 1:6))
meta <- read.csv("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/SampleMetaDataFull.csv")
idNames <- read.csv("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/SampleIdNames.csv")
meta <- merge(meta, idNames, by="Name")

pca <- merge(pca, meta, by="ID")
pca$Subspecie <- factor(pca$Subspecie, levels = c("hybrid", "capensis", "scutellata", "unicolor",
                                                  "caucasica", "mellifera", "carnica", "ligustica"))
ggplot(data = pca[pca$Subspecie != "hybrid",], aes(x = PCA1, y = PCA2, colour = Subspecie)) +
  geom_point(size=2) + 
  theme_bw(base_size=16) + 
  scale_colour_manual("", values = c("#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                    "#ad2793", "#310597"
                                    ))

pca$Lineage <- factor(pca$Lineage, levels = c("A", "O", "M", "C"))
ggplot(data = pca[pca$Subspecie != "hybrid",], aes(x = PCA1, y = PCA2, colour = Country)) +
  geom_point(size=2) + 
  theme_bw(base_size=16) 


pca <- pca[order(match(pca$Subspecie, c("hybrid", "capensis", "scutellata", "unicolor",
                                        "caucasica", "mellifera", "carnica", "ligustica"))),]
ggplot(data = pca, aes(x = PCA1, y = PCA3, colour = Subspecie)) +
  geom_point(size=2) + 
  theme_bw(base_size=16) + 
  scale_colour_manual("", values = c("grey",
                                     "#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                    "#ad2793", "#310597"
                                    ))

pca$Lineage <- factor(pca$Lineage, levels = c("A", "O", "M", "C", "H"))
ggplot(data = pca, aes(x = PCA1, y = PCA3, colour = Lineage)) +
  geom_point(size=2) + 
  theme_bw(base_size=16) 

table(pca$Subspecie)


CAhybrids <- pca[(pca$PCA2 < 150) & (pca$PCA2 > 50),]
CAhybrids$Country


#Plot for the WCGALP
pca12 <- ggplot(data = pca, aes(x = PCA1, y = PCA2, colour = Subspecie)) +
  geom_point(size=2) +
  scale_colour_manual("", values = c("grey",
                                     "#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                     "#ad2793", "#310597")) + 
  guides(colour = guide_legend(keywidth = unit(1, "cm"), nrow=2)) + 
  ggtitle("a")

pca13 <- ggplot(data = pca, aes(x = PCA1, y = PCA3, colour = Subspecie)) +
  geom_point(size=2) +
  scale_colour_manual("", values = c("grey",
                                     "#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                     "#ad2793", "#310597")) + 
  guides(colour = guide_legend(keywidth = unit(1, "cm"), nrow=2)) + 
  ggtitle("b")

pca23 <- ggplot(data = pca, aes(x = PCA2, y = PCA3, colour = Subspecie)) +
  geom_point(size=2) +
  scale_colour_manual("", values = c("grey",
                                     "#f7e225", "#feb72d", "#f07f4f",
                                     "#a8db34",
                                     "#228c8d",
                                     "#ad2793", "#310597")) + 
  guides(colour = guide_legend(keywidth = unit(1, "cm"), nrow=2)) + 
  ggtitle("c")


png("/home/x/JANA/HoneybeeGeno/tsinfer/HiFiTrees/PCA/PCASubspecies.png", res=800, width=170, height=80, units="mm")
combined <- pca12 + pca13 + pca23& theme(legend.position = "bottom")
combined + plot_layout(guides = "collect") 
dev.off()

