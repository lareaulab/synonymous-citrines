##Plotting the mRNA of the SL-citrine constructs in three different ways
library(dplyr)
library(tidyverse)
library(readxl)

setwd(".../data/stemloop/SL_mRNA/")

#load in data
datafile <- (".../J.WT.mRNA.csv.csv")
BioReps = read.csv( datafile, header=T)

#prepping data for plotting
citrine_order <- c("citMin", "cit0", "cit3", "cit6", "cit9", "citMax")
rates <- c( citMin = 164.6625671,
            cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919,
            cit9 = 302.6383892, 
            citMax = 388.8861452)
BioReps$elongation_rate <- rates[BioReps$Cit] 
ratios <- BioReps

wt.avg = aggregate( NormalizedLogTransform ~ Cit + Strain + elongation_rate, ratios[ratios$Strain == "WT",], mean)
j.avg = aggregate( NormalizedLogTransform ~ Cit + Strain + elongation_rate, ratios[ratios$Strain == "HC1J",], mean)

#main figure mRNA plot
pdf("HPvsWT_mrna_withaverages.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( ratios$elongation_rate, ratios$NormalizedLogTransform, 
      col = ifelse(ratios$Strain == "WT", "#3584c6", "#e95c64"),
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0,1.2),
      xlab = "",
      ylab = "citrine / mCherry\nrelative mRNA ratio"
)
axis( 1 )
axis( 2 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( wt.avg$elongation_rate, wt.avg$NormalizedLogTransform, lwd = 1.5, col = "#3584c6")
lines( j.avg$elongation_rate, j.avg$NormalizedLogTransform, lwd = 1.5, col = "#e95c64")

legend( "topright", pch = 20, 
        legend = c("no hairpin", "strong hairpin"), 
        col=c("#3584c6", "#e95c64"),
        cex=0.55,
        bty = "n")

dev.off()


#Plotting ratio of SL to WT mRNA for each citrine to look at the mRNA stabilization by codon optimality:

#average all of the biological replicates (A-C) of each type of sample (strain-cit)
BioRepAvg <- aggregate(NormalizedLogTransform ~ Strain + Cit, BioReps, mean)
#dividing hairpin by wt
divided <- BioRepAvg %>%
  pivot_wider(
    names_from = Strain,
    values_from = c(NormalizedLogTransform)
  )
divided$ratio <- divided$HC1J/divided$WT

rates <- c( cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919,
            cit9 = 302.6383892, 
            citMax = 388.8861452,
            citMin = 164.6625671)

divided$elongation_rate <- rates[divided$Cit] 

HP.WTratios <- divided$ratio
trend <- lm(HP.WTratios~rates)

#plot:
pdf("HPvsWT_mrna_ratios.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( divided$elongation_rate, divided$ratio, 
      col = "#e95c64",
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 2),
      xlab = "",
      ylab = "SL / WT\n mRNA ratio"
)
#abline(trend,col='#b2b2b2', lty = 2)

axis( 1 )
axis( 2, seq(0,2,0.33), labels = c(0.0, NA, 0.6, NA, 1.3, NA, 2))
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

dev.off()


#mRNA datapoints normalized to the average of citmax for each respective condition (WT vs SL):
#taking the avg citmax value for each strain category, then dividing each bio rep with the respective citmax:
JcitMax.avg <- BioRepAvg[BioRepAvg$Cit == "citMax"&BioRepAvg$Strain == "HC1J",]
JcitMax.avg <- as.numeric(JcitMax.avg[1,3])
WTcitMax.avg <- BioRepAvg[BioRepAvg$Cit == "citMax"&BioRepAvg$Strain == "WT",]
WTcitMax.avg <- as.numeric(WTcitMax.avg[1,3])

normedBioReps <- BioReps[, c("Strain", "Cit", "Clone")]

normedBioReps$normedmRNA <- ifelse(BioReps$Strain <= "HC1J",  BioReps$NormalizedLogTransform/JcitMax.avg, 
                                   ifelse(BioReps$Strain <= "WT", BioReps$NormalizedLogTransform/WTcitMax.avg, NA)
)
rates <- c( citMin = 164.6625671,
            cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919,
            cit9 = 302.6383892, 
            citMax = 388.8861452)
normedBioReps$elongation_time <- rates[normedBioReps$Cit] 

#With average lines:
normed.wt.avg = aggregate( normedmRNA ~ Cit + Strain + elongation_time, normedBioReps[normedBioReps$Strain == "WT",], mean)
normed.j.avg = aggregate( normedmRNA ~ Cit + Strain + elongation_time, normedBioReps[normedBioReps$Strain == "HC1J",], mean)

#plot
pdf("HPvsWT.mRNA.normedtocitMax.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( normedBioReps$elongation_time, normedBioReps$normedmRNA, 
      col = ifelse(normedBioReps$Strain == "WT", "#3584c6", "#e95c64"),
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0,5),
      xlab = "",
      ylab = "mRNA (cit/mCh)\nnormalized to citMax avg"
)
axis( 1 )
axis( 2 )
#axis( 2, seq(0,10,1.666666666), labels = c(0.0, NA, NA, 5, NA, NA, 10))
title( main = "TE normed to citMax", xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( normed.wt.avg$elongation_time, normed.wt.avg$normedmRNA, lwd = 1.5, col = "#3584c6")
lines( normed.j.avg$elongation_time, normed.j.avg$normedmRNA, lwd = 1.5, col = "#e95c64")

legend( "topright", pch = 20, 
        legend = c("no hairpin", "strong hairpin"), 
        col=c("#3584c6", "#e95c64"),
        cex=0.55,
        bty = "n")

dev.off()