library(fields)

datadir <- "data/stemloop/SL_protein/"
figdir <- "figures"

stemloop_data <- file.path( datadir, "normgated_data_bg_corrected.csv" )
citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"

utr_data <- "data/endogenous/Akirtava-TL-preprint-data-scer.csv" # each endog. utr expressing yfp
endog_data <- "data/endogenous/endogenous_scores_full_model.tsv" # ixnos full_cod_n3p2_nt_n9p8_rep0 model
uncharacterized <- "data/endogenous/uncharacterized_orfs.txt"

te_data <- "data/endogenous/lahtvee-2017-cellsystems-table-s6.csv" # protein per mRNA per time calculation for metabolic enzymes
id_convert_file <- "data/endogenous/uniprot_geneid.txt"

# fluorescence output from each citrine construct with each stemloop UTR variant
stemloop <- read.csv( stemloop_data, header=T, row.names=1)

# ixnos elongation scores for citrine constructs
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

medians <- aggregate( ratio ~ clone + cit + strain, stemloop, median) # median yellow/red ratio for each clone
medians <- medians[ medians$strain != "BY4743", ]

medians$elongation_time <- cit[ medians$cit, 1 ] # look up elongation times by citrine name

avg = aggregate( ratio ~ cit + strain + elongation_time, medians, mean)

# SL score: assume that citmin is the maximum output, take the citrine/mCherry ratio for citmin with each hairpin, 
# scale by citmin+no-hairpin to get an initiation score for each stem loop
sl_scores <- c( HC1j = avg$ratio[with( avg, which(cit == "citmin" & strain == "HC1j"))],
                HC1g = avg$ratio[with( avg, which(cit == "citmin" & strain == "HC1g"))],
                WT = avg$ratio[with( avg, which(cit == "citmin" & strain == "WT"))]
)
sl_scores <- sl_scores/max(sl_scores)

avg$sl_score <- sl_scores[ avg$strain ]
avg <- avg[ order(avg$sl_score), ]

#####
# ixnos scores of endogenous yeast genes (with ixnos full_cod_n3p2_nt_n9p8_rep0 model)
endog <- read.delim( endog_data, header=F)
names(endog) <- c("gene", "length", "tot_score")
endog$avg_score <- with( endog, tot_score / length )
# remove various problems from the endogenous genes file: uncharacterized ORFs, etc
unchar <- read.delim(uncharacterized, header=F)
unchar <- unchar$V1
endog <- endog[!(endog$gene %in% unchar),]
# remove mitochondrially encoded genes
endog <- endog[ grep("^Q", endog$gene, invert = T), ]
# remove "blocked_reading_frame" genes with internal stop codons (mutations relative to other strains)
blocked <- c("YDR134C", "YER109C", "YIL167W", "YIR043C", "YOL153C", "YOR031W")
endog <- endog[!(endog$gene %in% blocked),]
# remove shorter than 50 because of edge effects (scoring can't deal with the ends of genes well)
endog <- endog[ endog$length > 49, ]
#######

######
# TE measurements for ~1000 metabolic enzymes from Lahtvee et al. 
# Measured protein, protein decay rate, and mRNA level and estimated TE.
te <- read.delim( te_data, header=T, sep=",")
id_convert <- read.delim( id_convert_file, header=F)
names(id_convert) <- c("ID", "gene")
te <- merge( te, id_convert, by = "ID")
######

######
# Transcript leader (UTR) scores (RNA and YFP output) from Akirtava ... McManus preprint, endogenous UTR on YFP reporter
utr <- read.delim( utr_data, header=T, sep=",", row.names=1)
utr$gene <- do.call( rbind, sapply(rownames(utr), strsplit, "\\|"))[,1]

# this dataset has multiple UTRs for some genes. save the highest-output version for each gene
utr <- utr[order(utr$MeanYFP, decreasing = T),]
genes <- unique(utr$gene)

# grabbing the first match for each gene name in the sorted table should find the highest YFP for each
# save just the mean RNA and YFP columns
utr <- utr[ match( genes, utr$gene ), c( "gene", "MeanRNA", "MeanYFP")]
# also calculate protein per mRNA
utr$per_rna <- utr$MeanYFP / utr$MeanRNA
######

######
# combine datasets
# add endogenous CDS scores to reporter UTR scores
utr <- merge( utr, endog, by = "gene")

# just the subset with TE calculation from Lahtvee et al (metabolic enzymes)
endog_with_te <- merge( utr, te, by = "gene")
endog_with_te$per_rna[is.infinite(endog_with_te$per_rna )] <- NA
######

# plotting colors
slcols <- c( WT = "#000000", HC1g = "#BB90C2", HC1j = "#8D5B47" ) # colors for stemloops
citcols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")
pal <- colorRamp(c("blue", "red")) # heatmap
tepal <- colorRamp(c("blue", "red")) # color palette for coloring endogenous genes by Lahtvee TE

####################################################################
## panel 6A: output vs elongation time for each citrine construct ##
####################################################################

cairo_pdf( file.path( figdir, "stemloop_protein.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,1,2) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )

plot( avg$elongation_time[ avg$strain == "WT" ], 
      avg$ratio[ avg$strain == "WT" ], 
      type = "l", lwd = 1.5, col = slcols[ "WT" ],
      bty = "n",
      ylim = c(0,0.4),
      xlim = c(150, 400),
      axes = F,
      xlab = NA,
      ylab = "citrine / mCherry\nfluorescence ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by=50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( avg$elongation_time[ avg$strain == "HC1g" ], avg$ratio[ avg$strain == "HC1g" ], lwd = 1.5, col = slcols[ "HC1g" ])
lines( avg$elongation_time[ avg$strain == "HC1j" ], avg$ratio[ avg$strain == "HC1j" ], lwd = 1.5, col = slcols[ "HC1j" ])

points( medians$elongation_time, medians$ratio, pch = 20, col = slcols[ medians$strain ] )

# indicate which citrine is which
points( cit[ names(citcols), "time"], rep(-0.04, 6),
        col = citcols,
        pch = 17, cex = 1)

dev.off()


#####################################################################
## panel 6B: output vs initiation score for each citrine construct ##
#####################################################################

cairo_pdf( file.path( figdir, "stemloop_init_v_protein.pdf"), width = 1.5, height = 1.3, pointsize = 6.5)
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,1,3) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )

plot( sl_scores[ medians$strain ], medians$ratio, 
      col = citcols[medians$cit ],
      bty = "n",
      ylim = c(0,0.4),
      xlim = c(0.5, 1),
      pch = 20,
      axes = F,
      #xlab = NA,
      xlab = "initiation score",
      ylab = "citrine / mCherry\nfluorescence ratio"
)
axis( 1, lwd = 0.75 )
axis( 2, lwd = 0.75 )

sapply( names(citcols), function(cit) {
  lines( avg$sl_score[avg$cit == cit], avg$ratio[avg$cit == cit], lwd = 1.5, col = citcols[cit])
} )

# label which stem loop
points( sl_scores,  rep(-0.04,3),
        col = slcols[ names(sl_scores) ],
        pch = 18, cex = 1.5)

dev.off()

############################################################
## panel 6C: output vs initiation score + elongation time ##
############################################################

cairo_pdf( file.path( figdir, "init_v_elong_background.pdf"), width = 2, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,1,6) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )
# background: endogenous genes (just the ones with TE measurements, to match panel D)
plot( endog_with_te$avg_score * 238, # scaled to match the length of citrine for easy comparison with other plots
      endog_with_te$MeanYFP, 
      pch = 20, cex = 0.4,
      col = "gray90",
      bty = "n",
      ylim = c(0,1),
      xlim = c(150, 400),
      axes = F,
      xlab = NA,
      ylab = "initiation score"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by=50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = seq(0, 1, by=0.25), labels = c(0, NA, "0.5", NA, "1"))
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
# now plot the actual stem-loop data
points( avg$elongation_time, avg$sl_score,
     pch = 20, 
     col = rgb(pal( avg$ratio/max(avg$ratio)), max=255) )
# indicate which citrine is which
points( cit[ names(citcols), "time"], rep(-0.1, 6),
        col = citcols,
        pch = 17, cex = 1)
# and which stem loop
points( rep(130,3), sl_scores,
        col = slcols[ names(sl_scores) ],
        pch = 18, cex = 1.5)
# colorbar
imagePlot( legend.only=TRUE,
           zlim = c(0, max(avg$ratio)),
           col=colorRampPalette(c("blue","red"))(64),
           smallplot= c(.85, .875, 0.4, 0.75),
           legend.lab = "fluorescence",
           axis.args = list ( lwd = 0, lwd.ticks = 0.75, las = 3, at = c(0, 0.1, 0.2, 0.3),
                              tcl = -0.3, mgp = c(3, 0.5, 0)))
dev.off()



###################################################################
## panel 6D: endogenous TE vs initiation score + elongation time ##
###################################################################

scale <- max( log(endog_with_te$TE) )

cairo_pdf( file.path( figdir, "cds_utr_te.pdf"), width = 2, height = 1.3, pointsize = 6.5)
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,1,6) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )
plot( endog_with_te$avg_score * 238, # scaled to match the length of citrine for easy comparison with other plots
      endog_with_te$MeanYFP, 
      xlim = c(150,400),
#      xlim = c(175,325),
      ylim = c(0,1),
      pch = 20, cex = 0.4,
      col = rgb( tepal( log(endog_with_te$TE)/scale ), max=255 ),
      bty = "n",
      axes = F,
      xlab = NA,
      ylab = "initiation score" #"YFP output, McManus"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by=50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75, at = seq(0, 1, by=0.25), labels = c(0, NA, "0.5", NA, "1"))
title( xlab = "average elongation time\n(arbitrary units)", line = 4.5 )
# show where our citrines are
points( cit[ names(citcols), "time"], rep(-0.1, 6),
        col = citcols,
        pch = 17, cex = 1)
# and our stem loops
points( rep(130,3), sl_scores,
        col = slcols[ names(sl_scores) ],
        pch = 18, cex = 1.5)
imagePlot( legend.only=TRUE,
           zlim = c( 0, scale ),
           col=colorRampPalette(c("blue","red"))(64),
           smallplot= c(.85, .875, 0.4, 0.75),
           legend.lab = "log TE",
           axis.args = list ( lwd = 0, lwd.ticks = 0.75, las = 3,
                              tcl = -0.3, mgp = c(3, 0.5, 0)))
dev.off()



