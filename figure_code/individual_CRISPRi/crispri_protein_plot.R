# plot the flow results for the individual crispr validations

library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/CRISPRi_protein/"
figdir <- "figures"

data <- read.csv( file.path( datadir, "crispri_normgated_data_bg_corrected.csv"), header=T)

data <- data[ with( data, (strain %in% c(100:105, 111) & media == "tet") ), ] # crispri targets of interest from the CiBER-seq experiment

ratios <- aggregate( ratio ~ clone + cit + alias + gene, data, median ) # median of the green/red ratio of all points per sample
avgs <- aggregate( ratio ~  cit + alias + gene, ratios, mean ) # average of the median green/red ratios of three clones 
se <- aggregate( ratio ~  cit + alias + gene, ratios, function(x) sd(x)/sqrt(length(x)) )
avgs$se <- se$ratio

citmins <- avgs[ with( avgs, cit == "ci"), ]
cit9s <- avgs[ with( avgs, cit == "c9"), ]

c9ci_ratios <- citmins[ , c("gene", "alias") ]
c9ci_ratios$c9ci <- cit9s$ratio / citmins$ratio
c9ci_ratios$se <- c9ci_ratios$c9ci * sqrt( (with(cit9s, se/ratio)^2 + with(citmins, se/ratio)^2) )
  
c9ci_ratios$alias[ c9ci_ratios$gene == "HO" ] <- "control"

row.names(c9ci_ratios) <- c9ci_ratios$gene

ord <- c( "HO", c9ci_ratios$gene[ order(c9ci_ratios$c9ci, decreasing = T) ][1:6] )
c9ci_ratios <- c9ci_ratios[ ord, ]

cairo_pdf( file.path( figdir, "crispri_flow.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
#cairo_pdf( file.path( figdir, "crispri_flow.pdf"), width = 2, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,3,3) )
par( mar = c(8,6.5,3,1) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
mp <- barplot( c9ci_ratios$c9ci,
               ylim = c(0,0.5),
               space = 0.6,
               axes = F,
               #ylab = "slow citrine as\nfraction of fast citrine",
               ylab = "slow / fast citrine",
               xlab = NA,
               border = NA,
               col = "darkgrey"
               )
abline( h = c9ci_ratios["HO","c9ci"], lty = "dotted", col = "grey30")
axis(1, labels = c9ci_ratios$gene, at = mp, font = 3, tick = F, las = 2, lwd = 0.75 )
axis(2, lwd = 0.75 )
title( xlab = "CRISPRi target", line = 6 )
arrows( mp, with( c9ci_ratios, c9ci - se ), mp, with( c9ci_ratios, c9ci + se ), angle=90, code=3, length = 0.02)
dev.off()


