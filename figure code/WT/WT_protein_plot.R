## plot the three biological isolates that match the qPCR mRNA measurements
library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

setwd(".../data/WT/WT_protein/")

WT_protein_ratios <- read.csv("WT_protein_ratios.csv")
ratios <- as.data.frame(WT_protein_ratios)

speeds <- cit[ratios$Strain,1]
#cols <- c( MIN = "magenta3", Y000 = "royalblue2", Y333 = "green3", Y666 = "gold1", Y999 = "darkorange2", MAX = "red2")
cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

# protein output
pdf("../tunney_citrine_fluor_normgate.newspeeds2.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( 
  speeds, 
  ratios$ratio,
  col = cols[ratios$Strain],
  pch = 20,
  #cex = 0.6,
  #lwd = 1.5,
  axes = F,
  #      xlim = c(150,400),
  xlim = c(175,400),
  ylim = c(0, 0.3),
  xlab = "",
  ylab = "citrine / mCherry\nfluorescence ratio"
)
#axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 1, lwd = 0.75, at = seq(175, 400, by = 25), labels = c(NA, "200", NA, NA, NA, "300", NA, NA, NA, "400") )
axis( 2, lwd = 0.75, at = c(0,0.1,0.2,0.3) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()