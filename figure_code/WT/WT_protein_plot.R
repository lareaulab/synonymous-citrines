## plot the three biological isolates that match the qPCR mRNA measurements
library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

setwd("../../data/WT/WT_protein/")

WT_protein_ratios <- read.csv("WT_protein_ratios.csv")
ratios <- as.data.frame(WT_protein_ratios)

speeds <- ratios$citscore
cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

# protein output
pdf("../../../figures/tunney_citrine_fluor.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( 
  speeds, 
  ratios$ratio,
  col = cols[ratios$Strain],
  pch = 20,
  axes = F,
  xlim = c(150,400),
  ylim = c(0, 0.3),
  xlab = "",
  ylab = "citrine / mCherry\nfluorescence ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = c(0,0.1,0.2,0.3) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()
