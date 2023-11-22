# TE calculations for the individual crispr experiments

library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/"

ord <- c(FUN12 = 3, RPG1 = 2, HO = 1)
se <- function(x) sd(x)/sqrt(length(x))

## citrine fluorescence

data <- read.csv( file.path( datadir, "CRISPRi_protein/crispri_normgated_data_bg_corrected.csv"), header=T)
data <- data[ with( data, (strain %in% c(102, 105, 111) & media == "tet") ), ] # crispri targets of interest: HO control, RPG1, FUN12

fluor_ratios <- aggregate( ratio ~ clone + cit + gene, data, median ) # median of the green/red ratio of all points per sample
fluor_avgs <- aggregate( ratio ~  cit + gene, fluor_ratios, mean ) # average of the median green/red ratios of three clones 
fluor_avgs$se <- aggregate( ratio ~  cit + gene, fluor_ratios, se)$ratio

citmins <- fluor_avgs[ with( fluor_avgs, cit == "ci"), ]
cit9s <- fluor_avgs[ with( fluor_avgs, cit == "c9"), ]

fluor_c9_v_ci <- data.frame( gene = citmins$gene )
fluor_c9_v_ci$fc <- cit9s$ratio / citmins$ratio
fluor_c9_v_ci$se <- fluor_c9_v_ci$fc * sqrt( (with(cit9s, se/ratio)^2 + with(citmins, se/ratio)^2) )
# error propagation  

  

## mRNA

fun12 <- read.csv( file.path( datadir, "CRISPRi_mRNA/HOandFUN12_20230622/SCD_CRISPRi(HOandFUN12)_20230622 -  Quantification Summary_0.csv" ), header = T, row.names = 2)
rpg1 <- read.csv( file.path( datadir, "CRISPRi_mRNA/HOandRPG1_20230622/SCD_CRISPRi_20230622 -  Quantification Summary_0.csv" ), header = T, row.names = 2)

fun12map <- read.csv( file.path( datadir, "CRISPRi_mRNA/qPCR_map_20230622/HOandFUN12-Table 1.csv" ), header=T, row.names=1)
fun12labels <- c(t( fun12map[1:7,] ))[1:81]
fun12labels <- str_split(fun12labels, "-", simplify=T)
colnames(fun12labels) <- c("expt", "cit", "biorep", "techrep")

fun12 <- cbind(fun12, fun12labels)

rpg1map <- read.csv( file.path( datadir, "CRISPRi_mRNA/qPCR_map_20230622/HOandRPG1-Table 1.csv" ), header=T, row.names=1)
rpg1labels <- c(t( rpg1map[1:7,] ))[1:81]
rpg1labels <- str_split(rpg1labels, "-", simplify=T)
colnames(rpg1labels) <- c("expt", "cit", "biorep", "techrep")

rpg1 <- cbind(rpg1, rpg1labels)

# the average Cqs of the technical replicates for each sample (eg, HO cit9 bio rep 1)
fun12_techavg <- aggregate( Cq ~ Target + biorep + cit + expt, fun12[1:72,], mean)
rpg1_techavg <- aggregate( Cq ~ Target + biorep + cit + expt, rpg1[1:72,], mean)

mrna_techavg <- rbind( fun12_techavg, rpg1_techavg)

citrines <- grep("cit", mrna_techavg$Target)
mcherries <- grep("mCherry", mrna_techavg$Target)

dcq <- mrna_techavg[ citrines,c("biorep", "cit", "expt") ]
dcq$dCq <- mrna_techavg$Cq[ citrines ] - mrna_techavg$Cq[ mcherries ]
dcq$ratio <- 2^-dcq$dCq

mrna_avg <- aggregate( ratio ~ cit + expt, dcq, mean ) 
mrna_avg$se <- aggregate( ratio ~ cit + expt, dcq, se)$ratio

citmins <- mrna_avg[ with( mrna_avg, cit == "Citmin"), ]
cit9s <- mrna_avg[ with( mrna_avg, cit == "Cit9"), ]

mrna_c9_v_ci <- data.frame( gene = citmins$expt )
mrna_c9_v_ci$fc <- cit9s$ratio / citmins$ratio
mrna_c9_v_ci$se <- mrna_c9_v_ci$fc * sqrt( (with(cit9s, se/ratio)^2 + with(citmins, se/ratio)^2) )
# error propagation  

## assemble everything to plot
fluor_c9_v_ci <- fluor_c9_v_ci[ order( ord[fluor_c9_v_ci$gene]), ]
mrna_c9_v_ci <- mrna_c9_v_ci[ order( ord[mrna_c9_v_ci$gene]), ]

plot_fc <- data.frame( mRNA = mrna_c9_v_ci$fc, protein = fluor_c9_v_ci$fc, row.names = fluor_c9_v_ci$gene)
plot_se <- data.frame( mRNA = mrna_c9_v_ci$se, protein = fluor_c9_v_ci$se, row.names = fluor_c9_v_ci$gene)

write.csv(plot_fc, file.path( datadir, "ratios_for_CRISPRi_te.csv" )) 
write.csv(plot_se, file.path( datadir, "standard_error_for_CRISPRi_te.csv" ))
