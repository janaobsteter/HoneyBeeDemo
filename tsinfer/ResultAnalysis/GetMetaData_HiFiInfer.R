library(dplyr)

hifi <- read.table("~/JANA/HoneybeeGeno/SampleNames_HiFi.txt")
#nrow(hifi)
#hifi$name <- gsub("0_", "", hifi$V1)
#head(hifi)
meta <- read.csv("~/JANA/HoneybeeGeno/DavidWragg_DroneGeno/Meta.csv")
samplesAll <- read.table("~/EddieDir/SamplesAll.txt")



samplesAllRuns <- merge(samplesAll, runs, by.x="V1", by.y="GVCF")
nrow(samplesAll)
nrow(samplesAllRuns)

projectWork <- project$study_accession[project$Ploidy == "Workers"]
remove <- samplesAllRuns$V1[samplesAllRuns$study_accession %in% projectWork]

samplesAll <- samplesAll[!samplesAll$V1 %in% remove,]
write.table(samplesAll, "~/EddieDir/Drones.txt", quote=F, row.names=F)


nrow(samplesAll)
nrow(hifi)
length(intersect(samplesAll$V1, hifi$V1))
length(intersect(samplesAll$V1, meta$IID))


length(intersect(hifi$V1, meta$ID))
length(intersect(hifi$V1, meta$IID))
length(intersect(hifi$V1, meta$FID))

meta$VCFId <- ifelse(meta$IID %in% hifi$V1, meta$IID, meta$FID)
length(intersect(hifi$V1, meta$VCFId))
present <- intersect(hifi$V1, meta$VCFId)
missing <- hifi$V1[!hifi$V1 %in% intersect(hifi$V1, meta$VCFId)]
meta$IID[!meta$IID %in% intersect(hifi$V1, meta$VCFId)]


runs <- read.csv("~/JANA/HoneybeeGeno/DavidWragg_DroneGeno/Runs.csv")
project <- read.csv("~/JANA/HoneybeeGeno/DavidWragg_DroneGeno/Project.csv", sep=";")
runsGvcf <- merge()
head(project)
runsMissing <- runs[runs$GVCF %in% missing,]
nrow(runsMissing)
projectMissing <- project[project$sam %in% runsMissing$study_accession,]
nrow(projectMissing)
projectMissing

runsMissing <- merge(runsMissing, project, by="study_accession")
nrow(runsMissing)
head(runsMissing)
table(runsMissing$study_title)
write.table(runsMissing$GVCF[runsMissing$Ploidy == "Workers"], "~/JANA/HoneybeeGeno/DavidWragg_DroneGeno/WorkersGZVF.txt", quote=F, row.names=F)

# runsPresent <- runs[runs$GVCF %in% present,]
runsPresent <- merge(runsPresent, project, by="study_accession")
table(runsPresent$Ploidy)
runsPresent[runsPresent$Ploidy == "",]



# metaVCF <- meta[meta$VCFId %in% hifi$name,]
nrow(metaVCF)

head(metaVCF)
table(metaVCF$Country)
table(metaVCF$Type)


length(intersect(meta$VCFId, hifi$name))


#CHECK THE ORDER
sum(metaVCF$VCFId == asdSample$V1)
length(intersect(metaVCF$VCFId, asdSample$V1))

#ORDER THE COUTNRIES AND SUBSPECIES FILE
metaVCF <- metaVCF[match(asdSample$V1, metaVCF$VCFId),]
length(unique(metaVCF$VCFId))
sum(metaVCF$VCFId == asdSample$V1)

write.table(metaVCF$Country, "Countries.csv", row.names=F, col.names=F, quote=F)
write.table(metaVCF$Type, "Subspecies.csv", row.names=F, col.names=F, quote=F)

write.table(t(unique(metaVCF$Country)), "CountryList.csv", row.names=F, col.names=F, quote=F, sep=",")
write.table(t(unique(metaVCF$Type)), "SubspeciesList.csv", row.names=F, col.names=F, quote=F, sep=",")


# Write the coordinate file
head(metaVCF)
coordinate <- metaVCF[,c("Type", "VCFId", "Country"  )]
table(coordinate$Country)
coordWorld <- read.csv("~/JANA/Downloads/world_country_and_usa_states_latitude_and_longitude_values.csv")
sampleContries <- coordWorld[coordWorld$country %in% c("Canada", "Germany", "Finland",
                                                       "France", "Italy", "Madagascar",
                                                       "Mauritius", "Scotland", "Slovenia",
                                                       "Switzerland", "Sweden", "United States",
                                                       "ZAF"),]
write.table(sampleContries[,c("country_code", "latitude", "longitude", "country")], "SampleCountries.csv", quote=F, row.names=F)

sampleCountries <- read.csv("SampleCountries.csv")

#Merge with samples
head(metaVCF)
metaVCF <- merge(metaVCF, sampleCountries, by.x="Country", by.y = "country_code")
nrow(metaVCF)

cities = read.csv("~/JANA/Downloads/worldcities.csv")
table(cities$country)
cities$country[cities$country == "United States"] = "UnitedStates"
cities$country[cities$country == "South Africa"] = "SouthAfrica"
cities$country[cities$country == "United Kingdom"] = "Scotland"
citiesSamples <- cities[cities$country %in% metaVCF$country,]
table(citiesSamples$country) > table(metaVCF$country[metaVCF$country %in% cities$country])

unique(metaVCF$country[!metaVCF$country %in% cities$country])

coordinatesDF <- data.frame()
for (country in unique(metaVCF$country)) {
  print(country)
  countryDF <- metaVCF[metaVCF$country == country,]
  citiesDF <- cities[cities$country == country,]
  if (nrow(citiesDF) > nrow(countryDF)) {
    sampleCities <- citiesDF %>% sample_n(nrow(countryDF)) %>% select(c(lat, lng))
  } else if (nrow(citiesDF) < nrow(countryDF) & nrow(citiesDF != 0)) {
    sampleCities <- citiesDF %>% sample_n(nrow(countryDF), replace = T) %>% select(c(lat, lng))
  } else {
    sampleCities <- data.frame("lat" = countryDF$latitude, "lng" = countryDF$longitude)
  }
  coordinatesDF <- rbind(coordinatesDF, cbind(countryDF, sampleCities))
}

table(coordinatesDF$country)

#ORDER THE COUTNRIES AND SUBSPECIES FILE
nrow(coordinatesDF)
coordinatesDF <- coordinatesDF[match(asdSample$V1, coordinatesDF$VCFId),]
length(unique(coordinatesDF$VCFId))
sum(coordinatesDF$VCFId == asdSample$V1)

# Write coordinates
write.csv(metaVCF[,c("country", "VCFId", "latitude", "longitude")], "Coordinates.csv", quote=F, row.names=F)
write.csv(coordinatesDF[,c("country", "VCFId", "lat", "lng")], "Coordinates1.csv", quote=F, row.names=F)
