## generate figure re-analyzing the citrine mRNA qPCR data from Tunney, McGlincey, et al 2018

WT_mRNA_ratios <- read.csv("../../data/WT/WT_mRNA/WT_mRNA_ratios.csv")
ratios <- as.data.frame(WT_mRNA_ratios)

cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

pdf("../../figures/tunney_citrine_mrna.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
plot( ratios$citscore, ratios$ratio, 
      col = cols[ratios$strain],
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0,1.5),
      pch = 20,
      xlab = "",
      ylab = "citrine / mCherry\nrelative mRNA ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()
