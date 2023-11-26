##plotting TE with error bars for HO, RPG1, and FUN12 CRISPRi

datadir <- "data/individual_CRISPRi/"
figdir <- "figures"

plot_fc <- read.csv( file.path( datadir, "ratios_for_CRISPRi_te.csv" ), row.names = 1)
plot_se <- read.csv( file.path( datadir, "standard_error_for_CRISPRi_te.csv" ), row.names = 1 )

cols <- c(mRNA = "lightgrey", protein = "darkgrey")

pdf( file.path( figdir, "crispri_te.pdf"), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
#par( mar = c(6,6.5,5,3) )
par( mar = c(8,6.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp <- barplot( t(plot_fc), beside = T,
               space = c(0.2, 1.5),
               col = cols[ c("mRNA", "protein") ],
               ylim = c(0,0.6), border = NA, axes = F, ylab = NA, font = 3)
title( ylab = "slow citrine as\nfraction of fast citrine",
       xlab = "CRISPRi target" )
axis( 2, lwd = 0.75 )
## plot error bars
arrows( mp, t(plot_fc - plot_se), mp, t(plot_fc + plot_se ), angle=90, code=3, length = 0.02)

legend( "bottomright", legend = c("mRNA", "protein"), 
        fill = cols[ c("mRNA", "protein") ], border = NA, bty = "n", inset = c( -0.2, -1))
dev.off()
