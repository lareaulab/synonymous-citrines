# # Plotting the mRNA of the SL-citrine constructs
# library(dplyr)
# library(tidyverse)

datadir <- "data/stemloop/SL_mRNA/"
figdir <- "figures"

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

#load in data
datafile <- file.path( datadir, "J.WT.mRNA.csv" )
ratios <- read.csv( datafile, header=T)

ratios$elongation_time <- cit[ tolower(ratios$Cit), 1 ] # look up elongation times by citrine name

wt.avg <- aggregate( Cit_mCh ~ Cit + Strain + elongation_time, ratios[ratios$Strain == "WT",], mean)
j.avg <- aggregate( Cit_mCh ~ Cit + Strain + elongation_time, ratios[ratios$Strain == "HC1J",], mean)

cols <- c( WT = "#3584c6", HC1J = "#e95c64" )

#main figure mRNA plot
pdf( file.path( figdir, "stemloop_mrna.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( ratios$elongation_time, ratios$Cit_mCh,
      col = cols[ ratios$Strain ],
      pch = 20,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 1),
      xlab = NA,
      ylab = "citrine / mCherry\nrelative mRNA ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
lines( wt.avg$elongation_time, wt.avg$Cit_mCh, lwd = 1.5, col = cols[ "WT" ])
lines( j.avg$elongation_time, j.avg$Cit_mCh, lwd = 1.5, col = cols[ "HC1J" ])
legend( "topright", pch = 20, inset = c( -0.2, -0.6 ), 
        legend = c("no stem loop", "strong stem loop"), 
        col = cols, 
        cex = 0.8,
        bty = "n")
dev.off()



#Plotting ratio of SL to WT mRNA for each citrine to look at the mRNA stabilization by codon optimality:

#average all of the biological replicates (A-C) of each type of sample (strain-cit)
BioRepAvg <- aggregate( Cit_mCh ~ Strain + Cit + elongation_time, ratios, mean)

#dividing hairpin by wt
divided <- BioRepAvg %>%
  pivot_wider(
    names_from = Strain,
    values_from = c(Cit_mCh)
  )
divided$ratio <- divided$HC1J/divided$WT
divided$elongation_time <- cit[ tolower( divided$Cit), 1 ]

#plot:
pdf( file.path( figdir, "stemloop_mrna_ratios.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( divided$elongation_time, divided$ratio, 
      col = "darkgrey",
      pch = 20,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 2),
      xlab = "",
      ylab = "SL / WT\n mRNA ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()
