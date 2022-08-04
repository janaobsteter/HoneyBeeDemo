library(readr)
library(ggplot2)
edges1 <- read_delim("~/EddieDir/tsinfer/Trees/chr_1_Edges.txt")
edges3 <- read_delim("~/EddieDir/tsinfer/Trees/chr_3_Edges.txt")
edges12 <- read_delim("~/EddieDir/tsinfer/Trees/chr_12_Edges.txt")
edges16 <- read_delim("~/EddieDir/tsinfer/Trees/chr_16_Edges.txt")


summary(edges3$left)
ggplot(data = edges1, aes(x = left)) + geom_histogram(bins = 100)
ggplot(data = edges3, aes(x = left)) + geom_histogram(bins = 100) +
  geom_vline(xintercept =  11771679) + geom_vline(xintercept = 11781139)
ggplot(data = edges12, aes(x = left)) + geom_histogram(bins = 100)
ggplot(data = edges16, aes(x = left)) + geom_histogram(bins = 100) +
  geom_vline(xintercept =  3850299) + geom_vline(xintercept = 4068211)


trees1 <- read_delim("~/EddieDir/tsinfer/Trees/Tables/Chr1TreeInterval.csv")
trees3 <- read_delim("~/EddieDir/tsinfer/Trees/Tables/Chr3TreeInterval.csv")
trees3$y <- 1

ggplot(trees1, aes(x = start)) + geom_histogram(bins = 100)
ggplot(trees3, aes(x = start, y = y)) + geom_linerange(aes(xmin = start, xmax = end))


# Breakpoints
breakpoints <- data.frame()
for (chr in 1:16) {
  chrBreak <- read.csv(paste0("/home/janao/EddieDir/tsinfer/Trees/Tables/Chr", chr, "Breakpoints.csv"))
  chrBreak$Chr <- chr
  breakpoints <- rbind(breakpoints, chrBreak)
}

write.csv(breakpoints, "~/Documents/1Projects/HoneybeeDemography/tsinfer/Breakpoints.csv", quote=F, row.names=F)

ggplot(data = breakpoints[breakpoints$Chr], aes(x = X0)) + geom_histogram(bins = 300) + 
  facet_wrap(. ~ Chr, scales = "free")
  #geom_vline(xintercept =  11771679) + geom_vline(xintercept = 11781139)

ggplot(data = breakpoints[breakpoints$Chr == 10,], aes(x = X0)) + geom_histogram(bins = 300) + 
  facet_wrap(. ~ Chr, scales = "free") + 
  geom_vline(xintercept = 11506149) + geom_vline(xintercept = 11513111)

#Gene yellow-y
#Chromosome 10: 11506149-11513111 
break10 <- breakpoints[breakpoints$Chr == 10,]
yellowy <- break10[break10$X0 > 11506149 & break10$X0 < 11513111,]
yellowLength <- 11513111 - 11506149
yellowLength
nrow(yellowy)
yellowy[1,]
yellowy[nrow(yellowy),]

#Gene csd
#Chromosome 3: 11771679-11781139
break3 <- breakpoints[breakpoints$Chr == 3,]
csd <- break3[break3$X0 > 11771679 & break3$X0 < 11781139,]
csdLength <- 11781139 - 11771679
csdLength
nrow(csd)
csd[1,]
csd[nrow(csd),]
# 
