library(tidyverse)
library(dplyr)

setwd("/Users/liana/Dropbox/Berkeley/CodonCRISPR/Citrine Codon Paper and Final Data/figures data and code/figure 2 - hairpins, chimeras, start context/")

gated = read.csv("hairpin_flow_20221116_v005/normgated_data_bg_corrected.csv", header=T, row.names=1)

medians <- aggregate( ratio ~ clone + cit + strain, gated, median)

medians <- medians[ medians$strain != "BY4743", ]

# Add elongation rate for each citrine
rates <- c( citmin = 164.6625671,
            cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919, 
            cit9 = 302.6383892, 
            citmax = 388.8861452 )

medians$elongation_rate = 400 / rates[medians$cit] 
medians$elongation_time = rates[medians$cit] 

wt.avg = aggregate( ratio ~ cit + strain + elongation_rate + elongation_time, medians[medians$strain == "WT",], mean)
j.avg = aggregate( ratio ~ cit + strain + elongation_rate + elongation_time, medians[medians$strain == "HC1j",], mean)
g.avg = aggregate( ratio ~ cit + strain + elongation_rate + elongation_time, medians[medians$strain == "HC1g",], mean)

palette(c("#F8766D","#00BA38","#619CFF"))

pdf("hairpins_rate.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

plot( wt.avg$elongation_rate, wt.avg$ratio, 
      type = "l", lwd = 1.5, col = palette()[3],
      bty = "n",
      ylim = c(0,0.4),
      #xlim = c(150, 400),
      xlim = c(1, 2.5),
      xlab = "predicted elongation rate (a.u.)",
      ylab = "citrine / mCherry\nfluorescence ratio"
)

lines( g.avg$elongation_rate, g.avg$ratio, lwd = 1.5, col = palette()[1])
lines( j.avg$elongation_rate, j.avg$ratio, lwd = 1.5, col = palette()[2])
#lines( wt.avg$elongation_rate, wt.avg$ratio, lwd = 2, col = palette()[3])

points( medians$elongation_rate, medians$ratio, pch = 20, col = as.factor(medians$strain))

# legend( "topright", pch = 20, 
#         legend = c("no hairpin", "weak hairpin", "strong hairpin"), 
#         col=palette()[c(3,1,2)] )

dev.off()

pdf("hairpins_time.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

plot( wt.avg$elongation_time, wt.avg$ratio, 
      type = "l", lwd = 1.5, col = palette()[3],
      bty = "n",
      ylim = c(0,0.4),
      xlim = c(150, 400),
      xlab = "",
      ylab = "citrine / mCherry\nfluorescence ratio"
)
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( g.avg$elongation_time, g.avg$ratio, lwd = 1.5, col = palette()[1])
lines( j.avg$elongation_time, j.avg$ratio, lwd = 1.5, col = palette()[2])
#lines( wt.avg$elongation_time, wt.avg$ratio, lwd = 2, col = palette()[3])

points( medians$elongation_time, medians$ratio, pch = 20, col = as.factor(medians$strain))

# legend( "topright", pch = 20, bty = "n",
#        legend = c("no hairpin", "weak hairpin", "strong hairpin"), 
#        col=palette()[c(3,1,2)] )

text(300, c(0.475, 0.4, 0.325 ), labels = c("no hairpin", "weak hairpin", "strong hairpin"), col = palette()[c(3,1,2)], adj=0 )

dev.off()

