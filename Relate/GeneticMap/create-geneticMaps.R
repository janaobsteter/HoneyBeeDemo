#!/usr/bin/env Rscript

args <- commandArgs(T)

#args <- NULL
#args[1] <- "remapped_Liu_crossovers_NCBI.txt" # as output by NCBI remap
#args[2] <- "output path"

outpath <- args[2]

Liu <- read.table(args[1], header=F, strings=F)
Liu$chr <- as.numeric(as.factor(Liu[,1]))
colnames(Liu) <- c("GenBank", "start", "end", "chr")

# HAv3.1 (Assembled molecule size per chromosome) as reported by NCBI genome
chr.len <- c(27754200, 16089512, 13619445, 13404451, 13896941, 17789102, 14198698, 12717210,
             12354651, 12360052, 16352600, 11514234, 11279722, 10670842, 9534514, 7238532)

Mb <- 1000000

for(chr in seq(1,16)) {
  dat <- Liu[which(Liu$chr==chr),]
  breaks <- seq(1, chr.len[chr], by=Mb)
  if(breaks[length(breaks)] < chr.len[chr]) breaks[length(breaks)] <- chr.len[chr]
  recomb <- NULL
  for(i in 1:(length(breaks)-1)) recomb[i] <- nrow( dat[which(dat$start >= breaks[i] & dat$start <= breaks[i+1]),] )
  recomb[which(recomb==0)] <- 1
  # Adjust value because Liu et al sequenced 43 drones (100 meioses are observed for centiMorgan)
  recomb <- recomb/43*100
  cm <- NULL
  rr <- NULL
  for(i in 1:length(recomb))   cm <- c(cm, rep( recomb[i] / Mb, breaks[i+1]-breaks[i]))
  for(i in 1:length(recomb))   rr <- c(rr, rep( recomb[i], breaks[i+1]-breaks[i]))
  gp <- NULL
  gp <- c(0,cumsum(cm))
  map <- NULL
  map <- data.frame(CHR=as.factor(Liu[which(Liu$chr==chr),]$GenBank[1]), Physical=seq(1,chr.len[chr]), Genetic=gp, Rate=c(rr[1],rr), LG=chr, stringsAsFactors=F)
  # Check if any duplication in genetic map position
  if(length(names(table(duplicated(map$Genetic))))>1) cat("duplicate position(s) on chr ", chr, "\n")
  # Check approximate length in Mb
  write.table(map, paste(outpath, "/geneticMap_LG", chr, ".txt", sep=""), sep="\t", append=F, quote=F, row.names=F, col.names=F)
}

