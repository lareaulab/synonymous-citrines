library(tidyverse)
library(dplyr)

datadir <- "data/stemloop/SL_protein/"
figdir <- "figures"

gated = read.csv( file.path( datadir, "normgated_data_bg_corrected.csv" ), header=T, row.names=1)

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

medians <- aggregate( ratio ~ clone + cit + strain, gated, median)
medians <- medians[ medians$strain != "BY4743", ]

medians$elongation_time <- cit[ tolower(medians$cit), 1 ] # look up elongation times by citrine name

wt.avg = aggregate( ratio ~ cit + strain + elongation_time, medians[medians$strain == "WT",], mean)
j.avg = aggregate( ratio ~ cit + strain + elongation_time, medians[medians$strain == "HC1j",], mean)
g.avg = aggregate( ratio ~ cit + strain + elongation_time, medians[medians$strain == "HC1g",], mean)

cols <- c( WT = "#3584c6", HC1g = "#fbb615", HC1j = "#e95c64" )

pdf( file.path( figdir, "stemloop_protein.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

plot( wt.avg$elongation_time, wt.avg$ratio, 
      type = "l", lwd = 1.5, col = cols[ "WT" ],
      bty = "n",
      ylim = c(0,0.4),
      xlim = c(150, 400),
      axes = F,
      xlab = NA,
      ylab = "citrine / mCherry\nfluorescence ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by=50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( g.avg$elongation_time, g.avg$ratio, lwd = 1.5, col = cols[ "HC1g" ])
lines( j.avg$elongation_time, j.avg$ratio, lwd = 1.5, col = cols[ "HC1j" ])

points( medians$elongation_time, medians$ratio, pch = 20, col = cols[ medians$strain ] )

legend( "topright", pch = 20, inset = c( -0.2, -0.6 ), 
        legend = c("no stem loop", "weak stem loop", "strong stem loop"), 
        col = cols, 
        cex = 0.8,
        bty = "n")

dev.off()



