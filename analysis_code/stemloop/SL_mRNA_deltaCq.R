##Analysis of qPCR data for stemloop-citrine constructs

library(dplyr)
library(tidyverse)
library(readxl)

setwd("../../data/stemloop/SL_mRNA/")

#Assign the sample name to each well so that I can get a Cq value for them. 
#Function to be able to match each plate to its Cq value sheet, 
#then combining them into a single dataframe:

LabelData <- function(smap, Cqvalues){
names(smap)[1] <- "letter"
#pivot plate map to tidy format
smap <- smap %>%
pivot_longer(cols = 2:ncol(smap), names_to = "number") %>%
transmute(plate_well = paste(letter, number, sep = ""), name = value)
#remove empty wells
smap <- subset(smap, is.na(name) == FALSE)
#parse conditions from names
for (i in (1:nrow(smap))) {
temp_name <- str_split_fixed(smap[i,2],"-", n = 4)
split_name <- as_tibble(temp_name, .name_repair = "minimal")
smap[i,3] <- split_name[1,1]
smap[i,4] <- split_name[1,2]
smap[i,5] <- split_name[1,3]
smap[i,6] <- split_name[1,4]
}
#import Cq values from qPCR Data Sheet
names(Cqvalues)[1] <- "plate_well"
names(Cqvalues)[4] <- "Cq"
names(Cqvalues)[2] <- "PrimerTarget"
#remove unnecessary columns and empty wells
Cqvalues <- subset(Cqvalues, is.na("Cq")==FALSE, c(1,2,4))
#assign sample names to Cq values
qPCRdata <- left_join(Cqvalues, smap, by = "plate_well")
names(qPCRdata)[4] <- "Name"
names(qPCRdata)[5] <- "Strain"
names(qPCRdata)[6] <- "Cit"
names(qPCRdata)[7] <- "Clone"
names(qPCRdata)[8] <- "Replicate"
#convert Cq values to numeric class instead of character class
qPCRdata[, c("Cq")] <- sapply(qPCRdata[, c("Cq")], as.numeric)
return(qPCRdata)
}

#Organizing directories:
mainfolderpath <- ".../data/stemloop/SL_mRNA/"
minzerofolderpath <- ".../data/stemloop/SL_mRNA/Citmin_Cit0"
threesixfolderpath <- ".../data/stemloop/SL_mRNA/Cit3_Cit6"
ninemaxfolderpath <- ".../data/stemloop/SL_mRNA/Cit9_CitMax"
  
mapname <- "PlateMap.xlsx"
samplenames <- file.path(mainfolderpath, mapname)
  
minzerodatasheet <- file.path(minzerofolderpath, "CitminandCit0 WT vs HP -  Quantification Cq Results.xlsx")
threesixdatasheet <- file.path(threesixfolderpath, "Cit3andCit6 WT vs HP -  Quantification Cq Results.xlsx")
ninemaxdatasheet <- file.path(ninemaxfolderpath, "Cit9andCitX WT vs HP -  Quantification Cq Results.xlsx")

#labelling Cq values using above function
smap1 <- read_excel(samplenames, sheet = "citMinandcit0", col_names = TRUE)
Cqvalues1 <- read_excel(minzerodatasheet, col_names = TRUE)
  
smap2 <- read_excel(samplenames, sheet = "cit3andcit6", col_names = TRUE)
Cqvalues2 <- read_excel(threesixdatasheet, col_names = TRUE)
  
smap3 <- read_excel(samplenames, sheet = "cit9andcitMax", col_names = TRUE)
Cqvalues3 <- read_excel(ninemaxdatasheet, col_names = TRUE)
  
data1 <- LabelData(smap = smap1, Cqvalues = Cqvalues1)
data2 <- LabelData(smap = smap2, Cqvalues = Cqvalues2)
data3 <- LabelData(smap = smap3, Cqvalues = Cqvalues3)
  
qPCRdatamerged <- rbind(data1, data2, data3)
  
#reorder the columns to be more like the final table format
qPCRdataCleaned <- qPCRdatamerged[, c("Strain", "Cit", "Clone", "Replicate", "PrimerTarget", "Cq")]

#Average the technical replicates of each biological sample and take the standard deviation as well.
TechRepMeans <- aggregate(Cq ~ Strain + Cit + Clone + PrimerTarget, qPCRdataCleaned, mean)
TechRepMeansCqStd <- aggregate(Cq ~ Strain + Cit + Clone + PrimerTarget, qPCRdataCleaned, sd)
names(TechRepMeansCqStd)[5] <- "CqStd"

TechRepMeans <- merge(TechRepMeans, TechRepMeansCqStd, by = c("Strain", "Cit", "Clone", "PrimerTarget"))
TechRepMeans <- subset(TechRepMeans, Strain != "NTC", drop = TRUE)

#Using 2 as primer efficiency to do the deltaCq analysis
TechRepMeans$PrimerEfficiency <- 2
TechRepMeans$LogTransform <- TechRepMeans$PrimerEfficiency^(-TechRepMeans$Cq)

#Normalize the citrine to the mCherry for each biological sample: first separating the Cqs that target Ci and C9 
#from those that target mCherry by putting them into two separate dataframes:
mCherryCqsdf <- TechRepMeans[TechRepMeans$PrimerTarget == "mCh",]
CitrineCqsdf <- TechRepMeans[TechRepMeans$PrimerTarget != "mCh",]
SeparateCqs <- merge(CitrineCqsdf, mCherryCqsdf, by = c("Strain", "Cit", "Clone"))

SeparateCqs = select(SeparateCqs, -PrimerTarget.x, -PrimerTarget.y, -PrimerEfficiency.x, -PrimerEfficiency.y)

SeparateCqs <- SeparateCqs %>% rename("LogTransform.Cit" = "LogTransform.x")
SeparateCqs <- SeparateCqs %>% rename("Cq.Cit" = "Cq.x")
SeparateCqs <- SeparateCqs %>% rename("CqStd.Cit" = "CqStd.x")

SeparateCqs <- SeparateCqs %>% rename("LogTransform.mCh" = "LogTransform.y")
SeparateCqs <- SeparateCqs %>% rename("Cq.mCh" = "Cq.y")
SeparateCqs <- SeparateCqs %>% rename("CqStd.mCh" = "CqStd.y")

BioReps <- subset(SeparateCqs, select = c(Strain, Cit, Clone))
BioReps$NormalizedLogTransform <- SeparateCqs$LogTransform.Cit/SeparateCqs$LogTransform.mCh

write.csv(BioReps, "J.WT.mRNA.csv", row.names=FALSE) 

