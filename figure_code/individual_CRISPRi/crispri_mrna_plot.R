# mRNA and protein for the individual crispr experiments

library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/"
figdir <- "figures"

## citrine fluorescence

data <- read.csv( file.path( datadir, "CRISPRi_protein/crispri_normgated_data_bg_corrected.csv"), header=T)
data <- data[ with( data, (strain %in% c(102, 111) & media == "tet") ), ] # crispri targets of interest: HO control, RPG1

ratios <- aggregate( ratio ~ clone + cit + gene, data, median ) # median of the green/red ratio of all points per sample

# rearrange to match up cit9 (slow) and citmin (fast) pairs
prot_c9ci <- merge( ratios[ with( ratios, cit == "ci"), ], 
                    ratios[ with( ratios, cit == "c9"), ], 
                    by = c("gene", "clone"), suffixes = c(".ci", ".c9") )
prot_c9ci$c9ci <- with( prot_c9ci, ratio.c9 / ratio.ci)
prot_c9ci <- prot_c9ci[ , c("gene","clone","c9ci")]

avgs <- aggregate( c9ci ~ gene, prot_c9ci, mean ) # average of the cit9/citmin ratios of three pairs of clones


## citrine mRNA abundance
mrna <- read.csv( file.path( datadir, "CRISPRi_mRNA/HOandRPG1_20230622/SCD_CRISPRi_20230622 -  Quantification Summary_0.csv" ), header = T, row.names = 2)
qpcr_map <- read.csv( file.path( datadir, "CRISPRi_mRNA/HOandRPG1_20230622/HOandRPG1-Table 1.csv" ), header=T, row.names=1)
qpcr_labels <- c( t(qpcr_map[1:7,]) )[1:81]
qpcr_labels <- str_split( qpcr_labels, "-", simplify=T )
colnames(qpcr_labels) <- c( "expt", "cit", "biorep", "techrep" )

mrna <- cbind( mrna, qpcr_labels )

mrna$Abundance <- 2^-mrna$Cq

techmeans <- aggregate( Abundance ~ biorep + Target + expt + cit, mrna, mean)
mch_target <- techmeans[ techmeans$Target == "mCherry", ]
cit_target <- techmeans[ grep( "cit", techmeans$Target ), ]

mrna_merged <- merge( cit_target, mch_target, by = c("expt", "cit", "biorep"), suffixes = c(".cit", ".mch"))
mrna_merged$Relative <- mrna_merged$Abundance.cit / mrna_merged$Abundance.mch

mrna_citmins <- mrna_merged[ with( mrna_merged, cit == "Cit9"), ]
mrna_cit9s <- mrna_merged[ with( mrna_merged, cit == "Citmin"), ]

mrna_c9ci <- merge( mrna_citmins, mrna_cit9s, by = c("expt", "biorep"), suffixes = c(".c9", ".ci"))
mrna_c9ci$c9ci <- with( mrna_c9ci, Relative.c9 / Relative.ci )
mrna_c9ci <- mrna_c9ci[, c("expt", "biorep", "c9ci") ]

mrna_avgs <- aggregate( c9ci ~ expt, mrna_c9ci, mean ) # average of the cit9/citmin ratios of three pairs of clones

## to plot:
# protein datapoints are in prot_c9ci$c9ci
# protein averages are avgs$c9ci

# mRNA datapoints are mrna_c9ci$c9ci
# mRNA averages are mrna_avgs$c9ci

mrna_avgs_plot <- merge( mrna_avgs, avgs, by.x = "expt", by.y = "gene", suffixes = c(".mrna", ".protein"))
row.names(mrna_avgs_plot) <- mrna_avgs_plot$expt
mrna_avgs_plot <- mrna_avgs_plot[,2:3]

names(prot_c9ci) <- names(mrna_c9ci) # make the names match

mrna_plot <- merge( mrna_c9ci, prot_c9ci, by = c("expt", "biorep"), suffixes = c(".mrna", ".protein"))

# statistical significance
pvals <- c( mRNA = t.test( mrna_plot$c9ci.mrna[ mrna_plot$expt == "HO"], mrna_plot$c9ci.mrna[ mrna_plot$expt == "RPG1"] )$p.value,
            protein = t.test( mrna_plot$c9ci.protein[ mrna_plot$expt == "HO"], mrna_plot$c9ci.protein[ mrna_plot$expt == "RPG1"] )$p.value )
sig <- c( mRNA = "n.s.", protein = "n.s")
sig[ pvals < 0.5 ] <- "*"


cols <- c(mRNA = "lightgrey", protein = "darkgrey")
jitter <- c( -0.1, 0, 0.1 )

cairo_pdf( file.path( figdir, "rpg1_mrna_protein.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,4,0.5) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
max <- max(mrna_plot[, 3:4])
mp <- barplot( as.matrix(mrna_avgs_plot), beside = T,
               space = c(0.2, 0.8),
               col = cols[ rep( c("mRNA", "protein"), each = 2) ],
               names.arg = c( NA, NA),
               ylim = c(0,max), 
               border = NA, axes = F, ylab = NA, font = 3)
points( rep(mp[,1], each = 3) + jitter, mrna_plot$c9ci.mrna )
points( rep(mp[,2], each = 3) + jitter, mrna_plot$c9ci.protein )

title( ylab = "slow citrine as\nfraction of fast citrine")
axis( 2, lwd = 0.75 )
axis( 1, at = mp, labels = rep( c("HO", "RPG1"), 2), font = 3, tick = F, line = 0 )
axis( 1, at = apply(mp, 2, mean), labels = c("mRNA", "protein"), tick = F, line = 2 )

segments( mp[1,], max * 1.1, mp[2,], lwd = 0.75  )
text( apply(mp, 2, mean), max * 1.1, labels = sig, pos = 3)
dev.off()
