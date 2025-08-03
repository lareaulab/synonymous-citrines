library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/growth_defect_controls/"
figdir <- "figures"

datafile <- file.path( datadir, "crispri_growth_defect_controls_normgated_data.csv" )
data = read.csv( datafile, header=T)

ratios <- aggregate( ratio ~ clone + cit + gene + media, data, median) #medians of cit/mcherry ratios for each strain in each condition for each clone

citmins <- ratios[ with( ratios, cit == "ci"), ]
cit9s <- ratios[ with( ratios, cit == "c9"), ]

c9ci_ratios <- merge( citmins, cit9s, by = c("gene", "media", "clone"), suffixes = c(".ci", ".c9") )
c9ci_ratios$c9ci <- with( c9ci_ratios, ratio.c9 / ratio.ci)


avgs <- aggregate( c9ci ~ gene + media, c9ci_ratios, mean ) # average of the cit9/citmin ratios of three pairs of clones

c9ci_matrix <- cbind( CDC6 = avgs$c9ci[ with( avgs, gene == "CDC6") ], 
                      IRR1 = avgs$c9ci[ with( avgs, gene == "IRRI") ] )
row.names( c9ci_matrix ) <- avgs$media[ with( avgs, gene == "CDC6") ]

genes <- c("CDC6", "IRR1")

pvals <- sapply( unique(c9ci_ratios$gene), function(gene){ t.test( c9ci_ratios$c9ci[c9ci_ratios == gene & c9ci_ratios$media == "SCD"], 
                                                                   c9ci_ratios$c9ci[c9ci_ratios == gene & c9ci_ratios$media == "tet"] )$p.value })

sig <- pvals
sig[ pvals >= 0.05 ] <- "n.s."
names(sig)[names(sig) == "IRRI"] <- "IRR1" # typo in gated data

mediacols <- c( tet = "#39c0c4", SCD = "#99cccc")
labels <- c( tet = "induced", SCD = "uninduced" )

cairo_pdf( file.path( figdir, "growthcontrols.pdf" ), width = 1.5, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(4,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( c9ci_matrix, 
              beside = TRUE,
              col = mediacols[ row.names(c9ci_matrix) ],
              ylab = "slow citrine as\nfraction of fast citrine", 
              xlab = NA,
              names.arg = c(NA,NA),
              border = NA, axes = F,
              ylim = c(0, 0.5) )
axis( 1, lwd = 0.75, labels = genes, at = apply( mp, 2, mean), font = 3, tick = F )
axis( 2, lwd = 0.75 )

jitter <- rep( c(-0.1, 0, 0.1), 4 )
points( rep( as.vector(mp), each = 3 ) + jitter, 
        c9ci_ratios$c9ci,
        pch = 19, cex = 0.5, col = "grey30" )

mp <- as.data.frame( mp )
row.names(mp) <- rownames(c9ci_matrix)
names(mp) <- colnames(c9ci_matrix) 
max <- max(c9ci_ratios$c9ci)

sapply( colnames(mp), function(gene){
  segments(mp[1,gene], max*1.1, mp[2,gene], lwd = 0.75)
  text( mean(mp[,gene]), max*1.1, labels = sig[gene], pos = 3)
})

legend( 'topright', legend = labels[ row.names(c9ci_matrix) ],
        fill = mediacols[ row.names(c9ci_matrix) ],
        border = NA, bty = "n", inset = c( -0.2, -0.6 ))
dev.off()