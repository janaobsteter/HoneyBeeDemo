anc2o <- read.csv("EstsfsOutput_2o/AncestralAllele_COMBINED.csv", header=F)
colnames(anc2o) <- c("SnpPos", "ancAl")
anc3o <- read.csv("EstsfsOutput_3o/AncestralAllele_COMBINED.csv", header=F)
colnames(anc3o) <- c("SnpPos", "ancAl")
anc3o_10r <- read.csv("EstsfsOutput_3o_10r/AncestralAllele_COMBINED.csv", header=F)
colnames(anc3o_10r) <- c("SnpPos", "ancAl")
anc3o_rate6 <- read.csv("EstsfsOutput_rate6/AncestralAllele_COMBINED.csv", header=F)
colnames(anc3o_rate6) <- c("SnpPos", "ancAl")
# 
# anc2o <- read.csv("EstsfsOutput_2o/AncestralAllele1.csv")
# anc3o <- read.csv("EstsfsOutput_3o/AncestralAllele1_3o.csv")

anc3o$ancAl[anc3o$ancAl == ""] <- NA
anc3o_10r$ancAl[anc3o_10r$ancAl == ""] <- NA
anc3o_rate6$ancAl[anc3o_rate6$ancAl == ""] <- NA



print("Comparing 2outgroups vs 3outgroups, random=5.")
if (nrow(anc2o) == nrow(anc3o)) {
  print("The files have the same number of lines")
} else {
  print("The files DO NOT have the same number of lines!!!")
}

print(paste0("The number of rows is ", nrow(anc2o)))
print(paste0("The number of common positions is ", sum(anc2o$SnpPos == anc3o$SnpPos)))

print(paste0("Number of not NA lines in 2o file is ", nrow(anc2o[!is.na(anc2o$ancAl),])))
print(paste0("Number of not NA lines in 3o file is ", nrow(anc3o[!is.na(anc3o$ancAl),])))
per = sum(anc2o$ancAl == anc3o$ancAl, na.rm=T) / nrow(anc2o[!is.na(anc2o$ancAl),])
print(paste0("Number of common ancestral alleles is ", sum(anc2o$ancAl == anc3o$ancAl, na.rm=T), " or ", round(per,4)*100, "%"))


print("Comparing 3outgroups,random=5 vs 3outgroups, random=10.")
if (nrow(anc3o) == nrow(anc3o_10r)) {
  print("The files have the same number of lines")
} else {
  print("The files DO NOT have the same number of lines!!!")
}

print(paste0("The number of rows is ", nrow(anc3o)))
print(paste0("The number of common positions is ", sum(anc3o$SnpPos == anc3o_10r$SnpPos)))

print(paste0("Number of not NA lines in 3o file is ", nrow(anc3o[!is.na(anc3o$ancAl),])))
print(paste0("Number of not NA lines in 3o_10r file is ", nrow(anc3o_10r[!is.na(anc3o_10r$ancAl),])))
per = sum(anc3o$ancAl == anc3o_10r$ancAl, na.rm=T) / nrow(anc3o[!is.na(anc3o$ancAl),])
print(paste0("Number of common ancestral alleles is ", sum(anc3o$ancAl == anc3o_10r$ancAl, na.rm=T), " or ", round(per,4)*100, "%"))


print(paste0("The number of rows is ", nrow(anc3o_rate6)))
print(paste0("The number of common positions is ", sum(anc3o_rate6$SnpPos == anc3o$SnpPos)))

print(paste0("Number of not NA lines in 3orate6 file is ", nrow(anc3o_rate6[!is.na(anc3o_rate6$ancAl),])))
print(paste0("Number of not NA lines in 3o file is ", nrow(anc3o[!is.na(anc3o$ancAl),])))
per = sum(anc3o$ancAl == anc3o_rate6$ancAl, na.rm=T) / nrow(anc3o[!is.na(anc3o$ancAl),])
print(paste0("Number of common ancestral alleles is ", sum(anc3o$ancAl == anc3o_rate6$ancAl, na.rm=T), " or ", round(per,4)*100, "%"))



#which difer
# colnames(anc3o) <- c("SnpPos", "ancAl3o")
# combined <- merge(anc2o, anc3o, by="SnpPos")
# head(combined[combined$ancAl != combined$ancAl3o,])
# 
# dict <- read.table("EstSfs_Dict1.csv")
# dict[dict$V1 %in% c("carnica.LG1_1003128","carnica.LG1_1003855", "carnica.LG1_1005067"),]



#Check the NA ones
head(anc2o[is.na(anc2o$ancAl),])
dict <- read.table("EstSfs_Dict0.csv")
dict[dict$V1 %in% c("carnica.LG1_7841", "carnica.LG1_8071"),]


# Do the PCA of the missing sites
dict <- read.csv("EstSfs_DictTotal_PCA.csv", header=F)
head(dict)
pca <- prcomp(as.matrix(dict), scale. = TRUE)
head(pca$x)
nrow(pca$x)
missingSnps <- data.frame(missing = is.na(anc3o$ancAl))
colour <- c("red", "black")[missingSnps$missing + 1]
length(colour)

ggplot(data = as.data.frame(pca$x), aes(x = PC1, y = PC2, colour = colour)) + geom_point()
