##plotting TE with error bars for HO, RPG1, and FUN12 CRISPRi
setwd("../../data/individual_CRISPRi/")

plot_fc <- read.csv("ratios_for_CRISPRi_te.csv")
plot_se <- read.csv("standard_error_for_CRISPRi_te.csv")

cols <- c(mRNA = "lightgrey", protein = "darkgrey")

pdf("../../figures/crispri_te.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp <- barplot( t(plot_fc), beside = T,
               space = c(0.2, 1.5),
               col = cols,
               ylim = c(0,0.6), border = NA, axes = F, ylab = NA, font = 3)
title( ylab = "slow citrine as\nfraction of fast citrine",
       xlab = "CRISPRi target" )
axis(2)
## plot error bars
segments( x0 = mp[1,],
          y0 = plot_fc$mRNA - plot_se$mRNA,
          y1 = plot_fc$mRNA + plot_se$mRNA )
segments( x0 = mp[2,],
          y0 = plot_fc$protein - plot_se$protein,
          y1 = plot_fc$protein + plot_se$protein )
legend( "topright", legend = c("mRNA", "protein"), 
        fill = cols, border = NA, bty = "n", inset = c( 0, -0.5 ))
dev.off()
