#load useful libraries
library(dplyr)
library(tidyverse)

#set up the directories (this file will go in the plotting code)
datadir <- "data/stemloop/SL_protein/"
figdir <- "figures/"

#load in the data from the csv's in the SL folder:
SL_protein_ratios <- read.csv(file.path( datadir, "normgated_data_bg_corrected.csv" ))
SLratios <- as.data.frame(SL_protein_ratios)

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

medians <- aggregate( ratio ~ clone + cit + strain, SLratios, median)
medians <- medians[ medians$strain != "BY4743", ]

medians$elongation_time <- cit[ tolower(medians$cit), 1 ] # look up elongation times by citrine name

#select only the WT datapoints from those dfs
WTmedians <- medians[medians$strain == "WT", ]

#protein output plot
speeds <- WTmedians$elongation_time
cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

pdf( file.path( figdir, "efl_citrine_fluor.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( 
  WTmedians$elongation_time, 
  WTmedians$ratio,
  col = cols[WTmedians$cit],
  pch = 20,
  axes = F,
  xlim = c(150,400),
  ylim = c(0, 0.4),
  xlab = "",
  ylab = "citrine / mCherry\nfluorescence ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = c(0,0.1,0.2,0.3,0.4) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()