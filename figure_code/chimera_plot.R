# plot the citrine : mCherry fluorescence ratio of the two 'chimera' transcripts: 
# first half fast second half slow, and vice versa
# flow data is analyzed by ko_flow_processing_normgate.R

library(plyr)
library(tidyverse)

datadir <- "data/RQC_knockouts_and_chimeras/"
figdir <- "figures/"

namemap <- c( ci = 'citmin', c9 = 'cit9', Sphinx = 'slowfast', Lion = 'fastslow')

gated <- read.csv( paste0( datadir, "knockouts_chimeras_normgated_bg_corrected.csv" ), header=T)

ratios <- aggregate( ratio ~ sample + clone + cit + strain, gated, median)

ratios <- ratios[ with( ratios, strain %in% c('Lion', 'Sphinx', 'WT')), ]
ratios$cit <- namemap[ ratios$cit ]

averages <- aggregate( ratio ~ cit, ratios, mean) # average of the three points

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")


xpoints <- cit[ ratios$cit, ]

#desc <- c( ci = 'fast', c9 = 'slow', Sphinx = 'slow:fast', Lion = 'fast:slow')
desc <- c( citmin = 'fast', cit9 = 'slow', slowfast = 'slow:fast', fastslow = 'fast:slow')

cols <- c( citmin = 'magenta3', cit9 = 'darkorange2', slowfast = 'darkgray', fastslow = 'darkgray' )
firsthalf <- c( citmin = 'magenta3', cit9 = 'darkorange2', slowfast = 'darkorange2', fastslow = 'magenta3' )
secondhalf <- c( citmin = 'magenta3', cit9 = 'darkorange2', slowfast = 'magenta3', fastslow = 'darkorange2' )


cairo_pdf( paste0( figdir, "chimeras.pdf"), width = 2.01, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,7) ) # added 4 lines to right margin, each line is 0.065 inches
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( xpoints, ratios$ratio, 
      ylim = c(0,0.3),
      xlim = c(150,400),
      pch = 3, 
      cex = 0.6, lwd = 1.5,
      col = cols[ ratios$cit ],
      xlab = NA, 
      ylab = "citrine / mCherry\nfluorescence ratio",
      axes = F
      )
#axis( 1, lwd = 0.75, at = seq(150, 400, by = 50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 1, lwd = 0.75, at = seq(150, 350, by = 50), labels = c(NA, "200", NA, "300", NA) )
axis( 2, lwd = 0.75, at = c(0,0.1,0.2,0.3) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
rect( xleft = 370, ybottom = averages$ratio - 0.005, xright = 380, ytop = averages$ratio + 0.005, col = firsthalf[averages$cit], border = NA )
rect( xleft = 380, ybottom = averages$ratio - 0.005, xright = 390, ytop = averages$ratio + 0.005, col = secondhalf[averages$cit], border = NA )
text( x = 400, y = averages$ratio + 0.005, labels = paste(desc[averages$cit]), adj = 0, col = cols[averages$cit] )
dev.off()
