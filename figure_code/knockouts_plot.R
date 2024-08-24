# plot the results of knocking out several genes known to be involved in ribosome quality control 
# or related processes. flow data is analyzed by ko_flow_processing_normgate.R

library(plyr)
library(tidyverse)

datadir <- "data/RQC_knockouts_and_chimeras/"
figdir <- "figures/"

gated <- read.csv(paste0( datadir, "knockouts_chimeras_normgated_bg_corrected.csv"), header = T)

ratios <- aggregate( ratio ~ sample + clone + cit + strain, gated, median)

# isolates that turned out to be heterozygous for the knockout
hets <- c("Gcn1-C-ci", "Gcn1-C-c9", "Hel2Syh1-B-ci", "Hel2Syh1-B-c9" )
# chimera strains - will analyze in another figure
exclude <- c('4743', 'Lion', 'Sphinx')

ratios <- ratios[ !(ratios$strain %in% exclude ), ]
ratios[ ratios$sample %in% hets, "ratio" ] <- NA

xpoints <- c(jitter( rep(1:8, each = 6), amount = 0.15 ), rep(9, 2))

cols <- c("darkorange2", "magenta3")
collist <- c( rep( rep( cols, each = 3), 8), cols) 

labels <- c(expression(paste(italic('gcn1'), Delta)),
            expression(paste(italic('hel2'), Delta)),
            expression(paste(italic('hel2'), Delta, italic(' syh1'), Delta)),
            expression(paste(italic('mbf1'), Delta)),
            expression(paste(italic('rkr1'), Delta)),
            expression(paste(italic('smy2'), Delta)),
            expression(paste(italic('smy2'), Delta, italic(' syh1'), Delta)),
            expression(paste(italic('syh1'), Delta)),
            "WT")

cairo_pdf( paste0(figdir, "rqc_knockouts.pdf"), width = 3.5, height = 1.4, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(9,6.5,2,10) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( xpoints, ratios$ratio, 
      ylim = c(0,0.4),
      pch = 20,
      col = collist,
      xlab = NA, 
      ylab = "citrine / mCherry\nfluorescence ratio",
      axes = F
)
axis(2, lwd = 0.75)
axis(1, at = 1:9, labels = labels,
     col = NA,
     lwd = 0.75,
     las = 2)
text(10.5, ratios$ratio[50], label = "fast citrine", col = cols[2])
text(10.5, ratios$ratio[49], label = "slow citrine", col = cols[1])
dev.off()
