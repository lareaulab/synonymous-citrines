# plot the data that was tabulated by mrpa_analysis.R

setwd(".../data/ciBERseq/ciBERseq/")

fast <- read.table( "fast_mpra_results_byguide.txt", header = T, row.names = 1)
slow <- read.table( "slow_mpra_results_byguide.txt", header = T, row.names = 1)

eifs <- c("RPG1", "PRT1", "NIP1", "HCR1", "SUI1", "SUI2", "SUI3", "FUN12", "TIF1", "TIF11", "TIF2", "TIF3", "TIF34", "TIF35", "TIF4631", "TIF4632", "TIF5", "TIF6")

eif_rows <- unlist(sapply( eifs, function(x){grep(paste0(x, "$"), slow$gene, value = F)}))

eif_data <- data.frame( slow = slow$logFC[eif_rows], 
                        fast = fast$logFC[eif_rows], 
                        delta = slow$logFC[eif_rows] - fast$logFC[eif_rows], 
                        gene = fast$gene[eif_rows], 
                        row.names = row.names(slow[eif_rows,]))
eif_hits <- eif_data[ with( eif_data, delta > 0.5 & slow > 1 ), ]
eif_hits <- eif_hits[ with( eif_hits, order(slow) ), ]

# other <- c("YIH1", "TMA19")
# other_rows <- unlist(sapply( other, function(x){grep(paste0(x, "$"), slow$gene, value = F)}))
# other_data <- data.frame( slow = slow$logFC[other_rows],
#                           fast = fast$logFC[other_rows],
#                           delta = slow$logFC[other_rows] - fast$logFC[other_rows],
#                           gene = fast$gene[other_rows],
#                           row.names = row.names(slow[other_rows,]))
# other_hits <- other_data[ with( other_data, delta > 0.5 & slow > 1 ), ]



pdf("ciber-seq2.pdf", width = 3.33, height = 3.33, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )

plot(fast$logFC, slow$logFC, 
     xlab = "log fold change of reporter upon induction, fast citrine",
     ylab = "log fold change of reporter upon induction, slow citrine",
     xlim = c(-4.25, 2.5),
     ylim = c(-4.25, 2.5),
     bty = "n", 
     axes = F,
     asp = 1, pch = 20,
     col = "darkgrey")
axis( 1, at = c(-4, -2, 0, 2) )
axis( 2, at = c(-4, -2, 0, 2) )

abline( 0, 1, lty = 2, col = "grey" )
# bounds for including significant hits
segments( -2, 1, 0.5, 1, lty = 2, col = "red" )
segments( 0.5, 1, 2, 2.5, lty = 2, col = "red" )

# eIF hits in that region 
points(eif_hits$fast, eif_hits$slow, col = "black", pch = 20)
text(eif_hits$fast, eif_hits$slow, labels = eif_hits$gene,
     pos = c(4,2,4,4,4,2), col = "black")

# points(other_hits$fast, other_hits$slow, col = "darkblue", pch = 20)
# text(other_hits$fast, other_hits$slow, labels = other_hits$gene, pos = 2, col = "darkblue")
dev.off()


pdf("ciber-seq3.pdf", width = 2.5, height = 2.5, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )

plot(fast$logFC, slow$logFC, 
     xlab = "log fold change of reporter\nupon induction, fast citrine",
     ylab = "log fold change of reporter\nupon induction, slow citrine",
     xlim = c(-4.25, 2.5),
     ylim = c(-4.25, 2.5),
     bty = "n", 
     axes = F,
     asp = 1, pch = 20,
     col = "darkgrey")
axis( 1, at = c(-4, -2, 0, 2) )
axis( 2, at = c(-4, -2, 0, 2) )

abline( 0, 1, lty = 2, col = "grey" )
# bounds for including significant hits
segments( -2, 1, 0.5, 1, lty = 2, col = "red" )
segments( 0.5, 1, 2, 2.5, lty = 2, col = "red" )

# eIF hits in that region 
points(eif_hits$fast, eif_hits$slow, col = "black", pch = 20)
text(eif_hits$fast, eif_hits$slow, labels = eif_hits$gene,
     pos = c(4,2,4,4,4,2), col = "black")
dev.off()
