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
    transmute(Well = paste(letter, number, sep = ""), Name = value)
  
  #remove empty wells
  smap <- subset(smap, !is.na(Name))
  
  smap_split <- do.call(rbind, sapply( smap$Name, function(x)(strsplit(x, "-"))))
  smap <- smap %>% add_column(  Strain = smap_split[,1], Cit = smap_split[,2], Clone = smap_split[,3], Replicate = smap_split[,4]  )
  
  qPCRdata <- left_join( Cqvalues, smap, by = "Well" )
  return(qPCRdata)
}


#Organizing directories:

mapname <- "PlateMap.xlsx"

minzero <- "Citmin_Cit0/CitminandCit0 WT vs HP -  Quantification Cq Results_0.csv"
threesix <- "Cit3_Cit6/Cit3andCit6 WT vs HP -  Quantification Cq Results_0.csv"
ninemax <- "Cit9_CitMax/Cit9andCitX WT vs HP -  Quantification Cq Results_0.csv"

#labelling Cq values using above function
smap1 <- read_excel(mapname, sheet = "citMinandcit0", col_names = TRUE)
Cqvalues1 <- read_csv(minzero, col_names = T, col_select = c("Well", "Target", "Cq")) # tibble
  
smap2 <- read_excel(mapname, sheet = "cit3andcit6", col_names = TRUE)
Cqvalues2 <- read_csv(threesix, col_names = T, col_select = c("Well", "Target", "Cq")) # tibble

smap3 <- read_excel(mapname, sheet = "cit9andcitMax", col_names = TRUE)
Cqvalues3 <- read_csv(ninemax, col_names = T, col_select = c("Well", "Target", "Cq")) # tibble

data1 <- LabelData(smap = smap1, Cqvalues = Cqvalues1)
data2 <- LabelData(smap = smap2, Cqvalues = Cqvalues2)
data3 <- LabelData(smap = smap3, Cqvalues = Cqvalues3)
  
qPCRdatamerged <- rbind(data1, data2, data3)
  
#reorder the columns to be more like the final table format
qPCRdataCleaned <- qPCRdatamerged[, c("Strain", "Cit", "Clone", "Replicate", "Target", "Cq")]

#Average the technical replicates of each biological sample and take the standard deviation as well.
TechRepMeans <- aggregate(Cq ~ Strain + Cit + Clone + Target, qPCRdataCleaned, mean)
TechRepMeansCqStd <- aggregate(Cq ~ Strain + Cit + Clone + Target, qPCRdataCleaned, sd)
names(TechRepMeansCqStd)[5] <- "CqStd"

TechRepMeans <- merge(TechRepMeans, TechRepMeansCqStd, by = c("Strain", "Cit", "Clone", "Target"))
TechRepMeans <- subset(TechRepMeans, Strain != "NTC", drop = TRUE)

#Using 2 as primer efficiency to do the deltaCq analysis
TechRepMeans$PrimerEfficiency <- 2
TechRepMeans$Exp <- TechRepMeans$PrimerEfficiency^(-TechRepMeans$Cq)

#Normalize the citrine to the mCherry for each biological sample: first separating the Cqs that target Ci and C9 
#from those that target mCherry by putting them into two separate dataframes:
mCherryCqsdf <- TechRepMeans[TechRepMeans$Target == "mCh",]
CitrineCqsdf <- TechRepMeans[TechRepMeans$Target != "mCh",]
SeparateCqs <- merge(CitrineCqsdf, mCherryCqsdf, by = c("Strain", "Cit", "Clone"))

SeparateCqs = select(SeparateCqs, -Target.x, -Target.y, -PrimerEfficiency.x, -PrimerEfficiency.y)


SeparateCqs <- SeparateCqs %>% rename("Exp.Cit" = "Exp.x")
SeparateCqs <- SeparateCqs %>% rename("Cq.Cit" = "Cq.x")
SeparateCqs <- SeparateCqs %>% rename("CqStd.Cit" = "CqStd.x")

SeparateCqs <- SeparateCqs %>% rename("Exp.mCh" = "Exp.y")
SeparateCqs <- SeparateCqs %>% rename("Cq.mCh" = "Cq.y")
SeparateCqs <- SeparateCqs %>% rename("CqStd.mCh" = "CqStd.y")

BioReps <- subset(SeparateCqs, select = c(Strain, Cit, Clone))
BioReps$Cit_mCh <- SeparateCqs$Exp.Cit/SeparateCqs$Exp.mCh

write.csv(BioReps, "J.WT.mRNA.csv", row.names=FALSE) 

