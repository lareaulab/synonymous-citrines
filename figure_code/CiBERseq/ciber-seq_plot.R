# plot the data that was tabulated by mrpa_analysis.R

datadir <- "data/ciBERseq/ciBERseq/"
figdir <- "figures"

fast <- read.table( file.path( datadir, "fast_mpra_results_byguide.txt" ), header = T, row.names = 1)
slow <- read.table( file.path( datadir, "slow_mpra_results_byguide.txt" ), header = T, row.names = 1)

slowsig <- slow$adj.P.Val < 0.05
fastsig <-  fast$adj.P.Val < 0.05

up_cutoff <- log2(2)
diff_cutoff <- log2(2)
#diff_cutoff <- log2(1.5)

slowup <-  slow$logFC > up_cutoff & (slow$logFC - fast$logFC) > diff_cutoff

hits <- slowsig & slowup

eif_up_cutoff <- log2(2)
eif_diff_cutoff <- log2(1.5)
eif_up <-  slow$logFC > eif_up_cutoff & (slow$logFC - fast$logFC) > eif_diff_cutoff

# 9 genes with term "formation of cytoplasmic translation initiation complex"
eifs <- c( "FUN12", "TIF11", "TIF5", "RPG1", "TIF35", "HCR1", "TIF34", "NIP1", "PRT1" )
#eifs <- c("RPG1", "PRT1", "NIP1", "HCR1", "SUI1", "SUI2", "SUI3", "FUN12", "TIF1", "TIF11", "TIF2", "TIF3", "TIF34", "TIF35", "TIF4631", "TIF4632", "TIF5", "TIF6")
eif_rows <- slow$gene %in% eifs 
eif_hits <- eif_rows & eif_up
fun12 <- slow$gene == "FUN12" & slowup

cairo_pdf( file.path( figdir, "ciber-seq.pdf" ), width = 1.75, height = 1.75, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,2,1) )
par( oma = c(0,0.5,1,0) )

plot( fast$logFC, slow$logFC, 
      ylab = "log fold change, slow citrine",
      xlab = "log fold change, fast citrine",
      xlim = c(-4.25, 2.5),
      ylim = c(-4.25, 2.5),
      bty = "n", 
      axes = F,
      asp = 1, pch = 16, cex = 0.5, lwd = 0.5, 
      col = "grey80" )
axis( 1, at = c(-4, -2, 0, 2), lwd = 0.75 )
axis( 2, at = c(-4, -2, 0, 2), lwd = 0.75 )

abline( 0, 1, lty = 2, lwd = 0.5, col = "darkgrey" )

points( fast$logFC[slowsig], slow$logFC[slowsig], pch = 19, cex = 0.5, lwd = 0.5, col = "palegreen3" )
points( fast$logFC[fastsig], slow$logFC[fastsig], pch = 1, cex = 0.5, lwd = 0.5, col = "cornflowerblue" )

text( fast$logFC[hits], slow$logFC[hits], labels = slow$gene[hits], pos = 4, col = "palegreen3" )

dev.off()


cairo_pdf( file.path( figdir, "ciber-seq_eifs.pdf" ), width = 1.75, height = 1.75, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,2,1) )
par( oma = c(0,0.5,1,0) )

plot( fast$logFC, slow$logFC, 
      ylab = "log fold change, slow citrine",
      xlab = "log fold change, fast citrine",
      xlim = c(-4.25, 2.5),
      ylim = c(-4.25, 2.5),
      bty = "n", 
      axes = F,
      asp = 1, pch = 16, cex = 0.5, lwd = 0.5,
      col = "grey80" )
axis( 1, at = c(-4, -2, 0, 2), lwd = 0.75 )
axis( 2, at = c(-4, -2, 0, 2), lwd = 0.75 )

abline( 0, 1, lty = 2, lwd = 0.5, col = "darkgrey" )

# bounds for including significant hits
segments( -2, eif_up_cutoff, eif_up_cutoff - eif_diff_cutoff, eif_up_cutoff, lty = 2, col = "grey50", lwd = 0.5 ) # horizontal line at up_cutoff
segments( eif_up_cutoff - eif_diff_cutoff, eif_up_cutoff, 2, 2.5, lty = 2, col = "grey50", lwd = 0.5 ) # diagonal line at diff_cutoff

points( fast$logFC[eif_hits & !fun12], slow$logFC[eif_hits & !fun12], pch = 19, cex = 0.5, lwd = 0.5, col = "grey50" )
points( fast$logFC[fun12], slow$logFC[fun12], pch = 19, cex = 0.5, lwd = 0.5, col = "palegreen3" ) 

text( fast$logFC[eif_hits & !fun12], slow$logFC[eif_hits & !fun12], labels = slow$gene[eif_hits & !fun12], pos = 4, col = "grey50")
text( fast$logFC[fun12], slow$logFC[fun12], labels = slow$gene[fun12], pos = 4, col = "palegreen3" ) 

dev.off()


## write out table of all significant hits for supplemental data
## save the log2fc, adjusted p-val, gene id, and gene symbol
allsig <-  merge( slow[,c(1,5,7,9)], fast[,c(1,5,7,9)],
                  by = c("row.names","yorf", "gene"), suffixes=c(".slow", ".fast"))
names(allsig)[1] <- "guide"
allsig <- allsig[ allsig$adj.P.Val.slow < 0.05 | allsig$adj.P.Val.fast < 0.05, ]
write.table( x = allsig, 
             file = file.path( figdir, "ciber_sig_guides.tsv" ),
             row.names = F, col.names = TRUE, quote = F, sep = "\t")


## write out table of things used for the GO enrichment analysis:

go_out <-  merge( slow[,c(1,5,7,9)], fast[,c(1,5,7,9)],
                  by = c("row.names","yorf", "gene"), suffixes=c(".slow", ".fast"))
names(go_out)[1] <- "guide"

go_cutoffs <- go_out$logFC.slow > eif_up_cutoff & (go_out$logFC.slow - go_out$logFC.fast) > eif_diff_cutoff
go_out <- go_out[ go_cutoffs, ]
write.table( x = go_out, 
             file = file.path( figdir, "ciber_guides_go.tsv" ),
             row.names = F, col.names = TRUE, quote = F, sep = "\t")

