library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/x/EddieDir/Honeybees/YearCycleSimulation/")

# Apiary size = 20, 20K workers
a20w20 <- read.csv("Apiary20_20Kw/CsdVariability.txt")
a20w20 <- a20w20[!is.na(a20w20$Rep),]
a20w20 <- a20w20[a20w20$Rep != "Rep",]


a20w20$nCSDage0 <- as.numeric(a20w20$nCSDage0)
a20w20$totalCSDage0 <- as.numeric(a20w20$totalCSDage0)
a20w20$year <- as.factor(as.numeric(a20w20$year))
a20w20A <-  a20w20 %>% group_by(Rep, year) %>%  dplyr::summarise(Colony = mean(nCSDage0),
                                                                 Apiary = mean(totalCSDage0)) %>%
  pivot_longer(c(Colony, Apiary))


ggplot(data = a20w20A, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_wrap(. ~ name, scales = "free_y") +
  ggtitle("Apiary with 20 colonies, 20K workers/colony")

# Apiary size = 40, 20K workers
a20w60 <- read.csv("Apiary20_60Kw/CsdVariability.txt")
a20w60 <- a20w60[!is.na(a20w60$Rep),]
a20w60 <- a20w60[a20w60$Rep != "Rep",]


a20w60$nCSDage0 <- as.numeric(a20w60$nCSDage0)
a20w60$totalCSDage0 <- as.numeric(a20w60$totalCSDage0)
a20w60$year <- as.factor(as.numeric(a20w60$year))
a20w60$Rep <- as.factor(as.numeric(a20w60$Rep))
a20w60A <-  a20w60 %>% group_by(Rep, year) %>%  dplyr::summarise(Colony = mean(nCSDage0),
                                                                 Apiary = mean(totalCSDage0)) %>%
  pivot_longer(c(Colony, Apiary))


ggplot(data = a20w60A, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_wrap(. ~ name, scales = "free_y") +
  ggtitle("Apiary with 20 colonies, 60K workers/colony")


#Apiary with 20 colonies, 20K or 60K workers
a20w20A$Scenario <- "20K workers"
a20w60A$Scenario <- "60K workers"
a20 <- rbind(a20w20A, a20w60A)

ggplot(data = a20, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_grid(cols = vars(name), rows = vars(Scenario), scales = "free_y") +
  ggtitle("Apiary with 20 colonies, 60K workers/colony")


#########################3
#Only one chromosome
#############################
a20w60_1chr <- read.csv("Apiary20_60Kw_1chr/CsdVariability.txt")
a20w60_1chr <- a20w60_1chr[!is.na(a20w60_1chr$Rep),]
a20w60_1chr <- a20w60_1chr[a20w60_1chr$Rep != "Rep",]


a20w60_1chr$nCSDage0 <- as.numeric(a20w60_1chr$nCSDage0)
a20w60_1chr$totalCSDage0 <- as.numeric(a20w60_1chr$totalCSDage0)
a20w60_1chr$year <- as.factor(as.numeric(a20w60_1chr$year))
a20w60_1chr$Rep <- as.factor(as.numeric(a20w60_1chr$Rep))
a20w60_1chrA <-  a20w60_1chr %>% group_by(Rep, year) %>%  dplyr::summarise(Colony = mean(nCSDage0),
                                                                           Apiary = mean(totalCSDage0)) %>%
  pivot_longer(c(Colony, Apiary))


ggplot(data = a20w60_1chrA, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_wrap(. ~ name, scales = "free_y") +
  ggtitle("Apiary with 20 colonies, 60K workers/colony, 1 chromosome")

# Apiary wih 40 colonies, 60K workers per colony, 1 chromosome
a40w60_1chr <- read.csv("Apiary40_60Kw_1chr/CsdVariability.txt")
a40w60_1chr <- a40w60_1chr[!is.na(a40w60_1chr$Rep),]
a40w60_1chr <- a40w60_1chr[a40w60_1chr$Rep != "Rep",]


a40w60_1chr$nCSDage0 <- as.numeric(a40w60_1chr$nCSDage0)
a40w60_1chr$totalCSDage0 <- as.numeric(a40w60_1chr$totalCSDage0)
a40w60_1chr$year <- as.factor(as.numeric(a40w60_1chr$year))
a40w60_1chr$Rep <- as.factor(as.numeric(a40w60_1chr$Rep))
a40w60_1chrA <-  a40w60_1chr %>% group_by(Rep, year) %>%  dplyr::summarise(Colony = mean(nCSDage0),
                                                                           Apiary = mean(totalCSDage0)) %>%
  pivot_longer(c(Colony, Apiary))


ggplot(data = a40w60_1chrA, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_wrap(. ~ name, scales = "free_y") +
  ggtitle("Apiary with 20 colonies, 20K workers/colony")


a20w20A$Scenario <- "20 colonies, 20K workers"
a20w60A$Scenario <- "20 colonies, 60K workers"
a20w60_1chrA$Scenario <- "20 colonies, 60K workers, 1 chr"
a40w60_1chrA$Scenario <- "40 colonies, 60K workers, 1 chr"
aA <- rbind(a20w20A, a20w60A, a20w60_1chrA, a40w60_1chrA)

ggplot(data = aA, aes(x = year, y = value, colour = Rep, group = Rep)) +
  geom_line() +
  scale_colour_viridis(discrete = T) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_wrap(. ~ name + Scenario, ncol=2, dir="v", scales = "free_y")




a1chr <- rbind(a20w60_1chrA, a40w60_1chrA)
a1chrA <- a1chr %>%  group_by(year, name, Scenario) %>%  summarise(value = mean(value))
a1chrA$Scenario <- plyr::revalue(a1chrA$Scenario, c("40 colonies, 60K workers, 1 chr"="40 colonies", "20 colonies, 60K workers, 1 chr"="20 colonies"))
a1chr$Scenario <- plyr::revalue(a1chr$Scenario, c("40 colonies, 60K workers, 1 chr"="40 colonies", "20 colonies, 60K workers, 1 chr"="20 colonies"))
ggplot() +
  geom_line(data = a1chr, aes(x = year, y = value, group = Rep, colour = Rep), alpha = 0.2) +
  geom_line(data = a1chrA, aes(x = year, y = value, group = Scenario)) +
  theme_bw() +
  ylab("Number of distinct csd alleles") + xlab("Year") +
  facet_grid(rows = vars(name), scales = "free")
