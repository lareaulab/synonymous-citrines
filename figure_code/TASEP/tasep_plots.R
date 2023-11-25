options(stringsAsFactors=FALSE)
library("RColorBrewer")
pal <- brewer.pal(11, "RdGy")

datadir <- "data/TASEP/"
figdir <- "figures/"

init <- read.delim( file.path( datadir, "scan_init_withlife.txt" ), header=TRUE)
init$te_ratio <- with( init, (nprot_slow / mrna_slow) / (nprot_fast / mrna_fast) )
    
cols <- c( citmin = 'magenta3', cit9 = 'darkorange2')

pdf( file.path( figdir, "tasep_protein.pdf" ), width=1.75, height=1.3, pointsize = 6.5, useDingbats=FALSE, bg = "white")
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )

plot( init$init_rate * 60, init$nprot_ratio,
      xlim = c(0, 6),
      ylim = c(0,1),
      pch = 20,
      col = "darkgrey",
      xlab = "initiation [ribosomes / min]",
      ylab = "slow citrine as\nfraction of fast citrine",
      axes=F )
axis( 2, at = c(0, 0.5, 1.0), lwd = 0.75)
axis( 2, at = seq(0, 1.0, 0.25), labels = F, lwd = 0.75)
axis( 1, at = seq(0, 6, 2), lwd = 0.75)
axis( 1, at = seq(0, 6), labels = F, lwd = 0.75 )

abline(h=0.65, lwd=0.8, col=pal[[8]], lty = 2)
abline(h=0.32, lwd=0.8, col=pal[[10]], lty = 2)

text(x = 6.25, y = 0.65, adj = c(1, -0.2), labels="ratio expected from\nmRNA abundance alone", col=pal[[8]], cex = 0.7)
text(x = 6.25, y = 0.32, adj = c(1, 1.3), labels="observed ratio", col=pal[[10]], cex = 0.7)
dev.off()


################
pdf( file.path( figdir, "tasep_te.pdf" ), width=1.75, height=1.3, pointsize = 6.5, useDingbats=FALSE, bg = "white")
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )

plot( init$init_rate * 60, init$te_ratio,
      xlim = c(0, 6),
      ylim = c(0,1),
      pch = 20,
      col = "darkgrey",
      xlab = "initiation [ribosomes / min]",
      ylab = "slow citrine TE as\nfraction of fast citrine TE",
      axes=F )
axis( 2, at = c(0, 0.5, 1.0), lwd = 0.75)
axis( 2, at = seq(0, 1.0, 0.25), labels = F, lwd = 0.75)
axis( 1, at = seq(0, 6, 2), lwd = 0.75)
axis( 1, at = seq(0, 6), labels = F, lwd = 0.75 )

# abline(h=0.65, lwd=0.8, col=pal[[8]], lty = 2)
# abline(h=0.32, lwd=0.8, col=pal[[10]], lty = 2)
# 
# text(x = 6.25, y = 0.65, adj = c(1, -0.2), labels="ratio expected from\nmRNA abundance alone", col=pal[[8]], cex = 0.7)
# text(x = 6.25, y = 0.32, adj = c(1, 1.3), labels="observed ratio", col=pal[[10]], cex = 0.7)
dev.off()


################
pdf( file.path( figdir, "tasep_init_occlusion.pdf" ), width=1.75, height=1.3, pointsize = 6.5, useDingbats=FALSE, bg = "white")
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

plot( init$init_rate * 60, 
      1 - init$collide_fast,
      xlim = c(0,6),
      ylim= c(0,1),
      pch=20, 
      col=cols[["citmin"]],
      xlab="initiation [ribosomes / min]", 
      ylab="fraction successful\ninitiation attempts",
      axes=FALSE )

points(init$init_rate * 60, 
       1 - init$collide_slow,
       pch=20, 
       col=cols[["cit9"]] )

axis( 2, at=c(0, 0.5, 1.0), lwd = 0.75)
axis( 2, at=seq(0,1.0,0.25), labels=FALSE, lwd = 0.75)
axis( 1, at=seq(0,6,2), lwd = 0.75)
axis( 1, at=seq(0,6), labels=FALSE, lwd = 0.75 )

legend(x="bottomright", 
       col=cols, 
       bty="n",
       legend=c("fast", "slow"), 
       text.col=cols,
       pch=20, 
       cex = 0.8)
dev.off()
