# plot the distribution of endogenous gene scores and citrine scores
# from the ixnos full_cod_n3p2_nt_n9p8_rep0 model

figdir <- "figures"

endog_scores <- "data/endogenous/endogenous_scores_full_model.tsv" # ixnos full_cod_n3p2_nt_n9p8_rep0 model
uncharacterized <- "data/endogenous/uncharacterized_orfs.txt"
citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"

cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

endogenous <- read.delim(endog_scores, header=F)
names(endogenous) <- c("gene", "length", "score")
endogenous$avg <- with( endogenous, score/length)

# clean up genes a bit
# remove "Uncharacterized ORF" genes (these are not dubious orfs, but they include many very short genes of unknown function)
unchar <- read.delim(uncharacterized, header=F)
unchar <- unchar$V1
endogenous <- endogenous[!(endogenous$gene %in% unchar),]
# remove mitochondrially encoded genes
endogenous <- endogenous[ grep("^Q", endogenous$gene, invert = T), ]
# remove "blocked_reading_frame" genes with internal stop codons (mutations relative to other strains)
blocked <- c("YDR134C", "YER109C", "YIL167W", "YIR043C", "YOL153C", "YOR031W")
endogenous <- endogenous[!(endogenous$gene %in% blocked),]
# remove shorter than 50 because of edge effects (scoring can't deal with the ends of genes well)
endogenous <- endogenous[ endogenous$length > 49, ]


names <- c("citmin", "cit0", "cit3", "cit6", "cit9", "citmax")
nn.scores <- cit[names,]

cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")


cairo_pdf( file.path( figdir, "endogenous_histogram.pdf"), width = 1.75, height = 1, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

h = hist( endogenous$avg * 238,
          breaks = 30, 
          col = "gray80",
          border = NA,
          xlim = c(150,400),
          #freq = F,
          axes = F, 
          xlab = NA,
          ylab = "frequency",
          main = NA )
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = c(0,1000) )
points(nn.scores, rep(10, length(names)), col = cols[names], pch = 20)
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
#text(nn.scores[c(1,7)],c(max(h$density)/3, max(h$density)/3),labels = c("fastest","slowest"), col = c("magenta3", "red2"), cex = 0.7)
dev.off()