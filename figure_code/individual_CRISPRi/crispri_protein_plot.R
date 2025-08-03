# plot the flow results for the individual crispr validations

library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/CRISPRi_protein/"
figdir <- "figures"

data <- read.csv( file.path( datadir, "crispri_normgated_data_bg_corrected.csv"), header=T)

data <- data[ with( data, (strain %in% c(101:105, 111) & media == "tet") ), ] # crispri targets of interest from the CiBER-seq experiment

ratios <- aggregate( ratio ~ clone + cit + alias + gene, data, median ) # median of the green/red ratio of all points per sample

citmins <- ratios[ with( ratios, cit == "ci"), ]
cit9s <- ratios[ with( ratios, cit == "c9"), ]

c9ci_ratios <- merge( citmins, cit9s, by = c("gene", "alias", "clone"), suffixes = c(".ci", ".c9") )
c9ci_ratios$c9ci <- with( c9ci_ratios, ratio.c9 / ratio.ci)

avgs <- aggregate( c9ci ~ gene + alias, c9ci_ratios, mean ) # average of the cit9/citmin ratios of three pairs of clones

genes <- unique( c9ci_ratios$gene )

# do welch's t-test for the difference between the three citmin/cit9 values for each target gene vs HO control
pvals <- sapply( genes, function(x){ t.test( c9ci_ratios$c9ci[c9ci_ratios == x], c9ci_ratios$c9ci[c9ci_ratios == "HO"])$p.value })
sig <- rep( "", length(genes))
names(sig) <- names(pvals)
sig[ pvals < 0.05 ] <- "*"
sig[ "HO" ] <-  ""

# hack to put HO control first
sort_order = 1:6
names(sort_order) <- avgs$gene[order(avgs$c9ci, decreasing = T)]
sort_order["HO"] <- 0
avgs$order <- sort_order[ avgs$gene ]
avgs <- avgs[ order(avgs$order), ]

jitter <- rep( c(-0.1, 0, 0.1), length(genes) )
               
cairo_pdf( file.path( figdir, "crispri_flow.pdf"), width = 1.75, height = 1.3, pointsize = 6.5)
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,2,2) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
mp <- barplot( avgs$c9ci,
               ylim = c(0,0.5),
               space = 0.6,
               axes = F,
               ylab = "slow citrine as\nfraction of fast citrine",
               xlab = NA,
               border = NA,
               col = "grey80"
               )
abline( h = avgs$c9ci[avgs$gene == "HO"], lty = "dotted", col = "grey50")
axis(1, labels = avgs$gene, at = mp, font = 3, tick = F, las = 2, lwd = 0.75 )
axis(2, lwd = 0.75 )
title( xlab = "CRISPRi target", line = 6 )

# add individual points
par( xpd = NA)
mps <- mp[,1]
names(mps) <- avgs$gene
points( mps[ c9ci_ratios$gene ] + jitter, 
        c9ci_ratios$c9ci, 
        pch = 19, cex = 0.5, 
        col = "darkorchid4" )
        #col = "grey30" )

# add significance
text( mps[ avgs$gene ], 0.6, labels = sig[ avgs$gene ], cex = 1.5 )

xlim <- par('usr')[1:2] # save the x coordinates to make the next plot match
dev.off()


cols <- c( ci = "#9FDA3AFF", c9 = "#365C8DFF")

cairo_pdf( file.path( figdir, "crispri_flow_raw.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,1,2) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
par( xaxs = 'i' )
plot( mps[ratios$gene],
      ratios$ratio,
      xlim = xlim,
      ylim = c(0,0.1),
      axes = F,
      ylab = "citrine / mCherry\nfluorescence ratio",
      xlab = NA,
      col = cols[ ratios$cit]
)
axis(1, labels = names(mps), at = mps, font = 3, tick = F, las = 2, lwd = 0.75 )
axis(2, lwd = 0.75 )
title( xlab = "CRISPRi target", line = 6 )
dev.off()

