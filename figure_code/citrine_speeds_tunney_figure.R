# plot the decoding time of each citrine

library("fields")

datadir <- "data/codon_scores/"
figdir <- "figures/"

profile <- read.delim( file.path( datadir, "citrine_profiles_full_model.tsv" ), sep = "\t", row.names = 1, header = F) # from ixnos model: /mnt/lareaulab/rtunney/iXnos/expts/weinberg/lasagne_nn/full_cod_n3p2_nt_n9p8_rep0

citrines <- c("citmin", "cit0", "cit3", "cit6", "cit9", "citmax")
profile <- profile[rev(citrines),] # just the original 6

cols <- colorRampPalette(c("white", "blue"))(128)[1:96] # just use the lighter blue part of this palette
#citrinecols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")
citrinecols <- rev(viridis(8)[2:7])
names(citrinecols) <- c( "citmin", "cit0", "cit3", "cit6", "cit9", "citmax" )


cairo_pdf( file.path( figdir, "citrines_tunney.pdf"), width = 3.5, height = 1, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6,2,2) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
image( x = 1:238, y = 1:6, 
       z = t(profile),
       zlim = range(profile),
       col = cols,
       xlab = "codon position",
       ylab = NA,
       axes = F
)
axis(1, lwd = 0, lwd.ticks = 0.75)
axis(2, lwd = 0, line = 1, lwd.ticks = 0, labels = rev(c("fastest", "slowest")), at = c(1,6), las = 1)
par( xpd = NA )
imagePlot( legend.only=TRUE, 
           zlim = range(profile), 
           col=cols, horizontal=T, 
           smallplot = c(0.875, 0.95, 0.3, 0.35),
           legend.lab = "dwell time (au)",
           axis.args = list ( lwd = 0, lwd.ticks = 0.75, #at = c(1, 2), 
                              tcl = -0.3, mgp = c(3, 0.5, 0))
)
par( xpd = NA )
points( rep(-5, 6), 1:6, pch = 20, col = rev(citrinecols))
dev.off()