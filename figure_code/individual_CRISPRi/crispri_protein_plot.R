# plot the flow results for the individual crispr validations

library(plyr)
library(tidyverse)

setwd(".../data/individual_CRISPRi/CRISPRi_protein/FCS_files/")
data <- read.csv("crispri_normgated_data_bg_corrected.csv", header=T)

data.scd <- data[ with( data, (strain %in% c(100:105, 111) & media == "SCD") ), ] # uninduced - don't analyze these
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
#c9ci_ratios$alias <-  paste0( "(", c9ci_ratios$alias, ")" )

row.names(c9ci_ratios) <- c9ci_ratios$gene

ord <- c( "HO", c9ci_ratios$gene[ order(c9ci_ratios$c9ci, decreasing = T) ][1:6] )
c9ci_ratios <- c9ci_ratios[ ord, ]

pdf("crispri_flow.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
mp <- barplot( c9ci_ratios$c9ci,
               ylim = c(0,0.5),
               space = 0.6,
               axes = F,
               ylab = "slow citrine / fast citrine",
               xlab = NA,
               border = NA,
               col = "#ACAAFF"
#               col = "#CAA5E2"
               )
abline( h = c9ci_ratios["HO","c9ci"], lty = "dotted", col = "darkgrey")
axis(1, labels = c9ci_ratios$gene, at = mp, font = 3, tick = F, las = 2)
axis(2)
title( xlab = "CRISPRi target", line = 6 )
arrows( mp, with( c9ci_ratios, c9ci - se ), mp, with( c9ci_ratios, c9ci + se ), angle=90, code=3, length = 0.02)
dev.off()


pdf("crispri_flow_aliases.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,5,3,1) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
mp <- barplot( c9ci_ratios$c9ci,
               ylim = c(0,0.5),
               space = 0.6,
               axes = F,
               ylab = "slow citrine / fast citrine",
               xlab = NA,
               border = NA,
               col = "#ACAAFF"
#               col = "#CAA5E2"
#               col = "#8583E2"
)
abline( h = c9ci_ratios["HO","c9ci"], lty = "dotted", col = "darkgrey")
# axis(1, labels = c9ci_ratios$gene, at = mp - 0.3, font = 3, tick = F, cex.axis = 1, las = 2)
# axis(1, labels = c9ci_ratios$alias, at = mp + 0.3, tick = F, cex.axis = 0.6, las = 2)
axis(1, labels = c9ci_ratios$gene, at = mp, font = 3, tick = F, cex.axis = 1)
axis(1, labels = c9ci_ratios$alias, at = mp, tick = F, cex.axis = 0.8, line = 2)
axis(2)
title( xlab = "CRISPRi target", line = 6 )
arrows( mp, with( c9ci_ratios, c9ci - se ), mp, with( c9ci_ratios, c9ci + se ), angle=90, code=3, length = 0.02)
dev.off()


pdf("crispri_flow_all.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
# ugly hack to match x axis
mp <- barplot( c9ci_ratios$c9ci,
               ylim = c(0,0.08),
               space = 0.6,
               axes = F,
               ylab = "citrine / mCherry\nfluorescence",
               xlab = NA,
               border = NA,
               col = NA)
axis(2)
axis(1, labels = ord, at = mp, font = 3, tick = F, las = 2)
title( xlab = "CRISPRi target", line = 6 )
points( xpts, ratios$ratio,
      col = cols,
      pch = 20 )
dev.off()



row.names(mp) = ord
ratios.scd <- aggregate( ratio ~ clone + cit + alias + gene, data.scd, median ) # median of the green/red ratio of all points per sample

cols <- rep( c("darkorange2", "magenta3"), each = 3)
cols2 <- rep( c("#FED976", "plum"), each = 3)
xpts <- rep(mp[unique(ratios$gene),1], each=6)

pdf("crispri_flow_all_withscd.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp <- barplot( c9ci_ratios$c9ci,
               ylim = c(0,0.08),
               space = 0.6,
               axes = F,
               ylab = "citrine / mCherry\nfluorescence",
               xlab = NA,
               border = NA,
               col = NA)
axis(2)
axis(1, labels = ord, at = mp, font = 3, tick = F, las = 2)
title( xlab = "CRISPRi target", line = 6 )
points( xpts - 0.05, ratios.scd$ratio,
        col = cols2,
        pch = 20 )
points( xpts + 0.05, ratios$ratio,
      col = cols,
      pch = 20 )
dev.off()




# data for the three strains we followed up on
fhr <- ratios[ with( ratios, gene %in% c("HO", "FUN12", "RPG1")), ]
fhr_avgs <- aggregate( ratio ~  cit + alias + gene, fhr, mean ) # average of the median green/red ratios of three clones 
fhr_se <- aggregate( ratio ~  cit + alias + gene, fhr, function(x) sd(x)/sqrt(length(x)) )
fhr_avgs$se <- fhr_se$ratio
fhr <- rbind( fhr[ with( fhr, gene == "HO"), ],
              fhr[ with( fhr, gene == "RPG1"), ],
              fhr[ with( fhr, gene == "FUN12"), ])
fhr_avgs <- rbind( fhr_avgs[ with( fhr_avgs, gene == "HO"), ],
                   fhr_avgs[ with( fhr_avgs, gene == "RPG1"), ],
                   fhr_avgs[ with( fhr_avgs, gene == "FUN12"), ])

fhr_fluor_ratios <- data.frame( HO = fhr_avgs$ratio[with( fhr_avgs, gene == "HO")],
                                RPG1 = fhr_avgs$ratio[with( fhr_avgs, gene == "RPG1")],
                                FUN12 = fhr_avgs$ratio[with( fhr_avgs, gene == "FUN12")])
row.names(fhr_fluor_ratios) <- fhr_avgs$cit[with( fhr_avgs, gene == "HO")]

fhr_mrna_ratios <- data.frame( HO = bio_avgs$dcq[bio_avgs$expt == "HO"], RPG1 = bio_avgs$dcq[bio_avgs$expt == "RPG1"], FUN12 = bio_avgs$dcq[bio_avgs$expt == "FUN12"])
row.names(fhr_mrna_ratios) <- bio_avgs$Target[with( bio_avgs, expt == "HO")]
fhr_mrna_ratios <- 2^-fhr_mrna_ratios

fhr_te <- fhr_fluor_ratios / fhr_mrna_ratios
