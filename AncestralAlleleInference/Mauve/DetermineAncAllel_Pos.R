# This is a script to determinethe ancestral allele from the est-sfs output, combine it with the name and save it back into a file
#setwd("~/EddieDir/Honeybees/MultipleGenomeAlignment/EstSfs/")
#Est-sfs uSFS
#sfs1 <- as.data.frame(t(read.csv("outputEstsfs1.txt", header=F)))
#ggplot(sfs1, aes(x=V1)) + geom_histogram(bins=100)


args = commandArgs(trailingOnly=TRUE)
cycle = args[1]

# Pvalues
p <- read.csv(paste0("EstsfsOutput_2o/output-", cycle, "-pvalues.txt"), sep=" ", skip=8, header=F)
colnames(p) <- c("Site", "Code", "PmajAnc", "A", "C", "G", "T")       
p <- p[,-8]



# Define the functio to determine the ancestral allele
maxProbAllele <- function(x) {
  if (!is.nan(sum(x[4:7]))) {
    ancAl <- c("A", "C", "G", "T")[which.max(x[4:7])]
  } else {
    ancAl <- NA
  }
  return(ancAl)
}

# namesProbAllele <- function(x) {
#   if (!is.nan(sum(x[4:7]))) {
#     site <- x[1]
#   } else {
#     site = 
#   }
#   
# }


# Write the ancestral allele to the matrix
p$ancAl <- unlist(apply(p, MARGIN=1, FUN=function(x) maxProbAllele(x)))
#p$Site <- unlist(apply(p, MARGIN=1, FUN=function(x) namesProbAllele(x)))


# Read in the names
names <- read.table(paste0("EstSfsNames", cycle, ".csv"))
# Make sure the length of the SNPs and the names is the same --> if so, add ancestral alleles to the names
if (nrow(names) == nrow(p)) {
  names$ancAl <- p$ancAl
  colnames(names) <- c("SnpPos", "ancAl")
  write.csv(names, paste0("EstsfsOutput_2o/AncestralAllele", cycle, "_2o.csv"), row.names=F)
} else {
  return(paste0("Length of names and SNP for cycle", cycle, " does not match!!!"))
}

