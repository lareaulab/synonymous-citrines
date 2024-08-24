# plot the data that was tabulated by mrpa_analysis.R

datadir <- "data/ciBERseq/ciBERseq/"
figdir <- "figures"

fast <- read.table( file.path( datadir, "fast_mpra_results_byguide.txt" ), header = T, row.names = 1)
slow <- read.table( file.path( datadir, "slow_mpra_results_byguide.txt" ), header = T, row.names = 1)

eifs <- c("RPG1", "PRT1", "NIP1", "HCR1", "SUI1", "SUI2", "SUI3", "FUN12", "TIF1", "TIF11", "TIF2", "TIF3", "TIF34", "TIF35", "TIF4631", "TIF4632", "TIF5", "TIF6")

eif_rows <- unlist(sapply( eifs, function(x){grep(paste0(x, "$"), slow$gene, value = F)}))

eif_data <- data.frame( slow = slow$logFC[eif_rows], 
                        fast = fast$logFC[eif_rows], 
                        delta = slow$logFC[eif_rows] - fast$logFC[eif_rows], 
                        gene = fast$gene[eif_rows], 
                        row.names = row.names(slow[eif_rows,]))
eif_hits <- eif_data[ with( eif_data, delta > 0.5 & slow > 1 ), ]
eif_hits <- eif_hits[ with( eif_hits, order(slow) ), ]


cairo_pdf( file.path( figdir, "ciber-seq.pdf" ), width = 1.75, height = 1.75, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
#par( mar = c(7,6.5,1,1) )
par( mar = c(6,6.5,2,1) )
par( oma = c(0,0.5,1,0) )

plot(fast$logFC, slow$logFC, 
     #xlab = NA,
     #ylab = "log fold change of reporter\nupon induction, slow citrine",
     ylab = "log fold change, slow citrine",
     xlab = "log fold change, fast citrine",
     xlim = c(-4.25, 2.5),
     ylim = c(-4.25, 2.5),
     bty = "n", 
     axes = F,
     asp = 1, pch = 20,
     col = "darkgrey")
axis( 1, at = c(-4, -2, 0, 2), lwd = 0.75 )
axis( 2, at = c(-4, -2, 0, 2), lwd = 0.75 )
#title( xlab = "log fold change of reporter\nupon induction, fast citrine", line = 4.5 )

abline( 0, 1, lty = 2, col = "grey" )
# bounds for including significant hits
segments( -2, 1, 0.5, 1, lty = 2, col = "red" )
segments( 0.5, 1, 2, 2.5, lty = 2, col = "red" )

# eIF hits in that region 
points(eif_hits$fast, eif_hits$slow, col = "black", pch = 20)
text(eif_hits$fast, eif_hits$slow, labels = eif_hits$gene,
     pos = c(4,2,4,4,4,2),  col = "black")
dev.off()
