#!/usr/bin/env Rscript

args <- commandArgs(T)

#args <- NULL
#args[1] <- "imputed.vcf_1.map"
#args[2] <- "geneticMap_LG1.txt"

vcf <- read.table(args[1], header=F, stringsAsFactors=F)
names(vcf) <- c("LG", "POS")
map <- read.table(args[2], header=F, stringsAsFactors=F)
names(map) <- c("ACC", "POS", "GEN", "RATE", "LG")
vcf$GEN <- map[match(vcf$POS, map$POS),]$GEN
write.table(vcf, gsub(".map", ".genmap", args[1]), row.names=F, col.names=F, append=F, quote=F)
