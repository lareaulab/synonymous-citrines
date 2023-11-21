## generate figure re-analyzing the citrine mRNA qPCR data from Tunney, McGlincey, et al 2018

setwd("../../data/WT/")

WT_protein_ratios <- read.csv("WT_protein/WT_protein_ratios.csv")
protein <- as.data.frame(WT_protein_ratios)
names(protein) <- tolower( names( protein ))


WT_mRNA_ratios <- read.csv("WT_mRNA/WT_mRNA_ratios.csv")
mrna <- as.data.frame(WT_mRNA_ratios)

te <- merge( protein, mrna, by = c("isolate", "strain", "citscore"), suffixes = c(".prot", ".mrna") )
te$te <- with( te, ratio.prot / ratio.mrna)

cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

pdf("../../figures/tunney_citrine_te.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( te$citscore, te$te,
      col = cols[te$strain],
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 0.3),
      pch = 20,
      xlab = "",
      ylab = "citrine / mCherry\nrelative mRNA ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = c(0,0.1,0.2,0.3) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()
