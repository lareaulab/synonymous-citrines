# new citrine constructs for this paper - TE measurements combining flow and qPCR data

datadir <- "data/stemloop/"
figdir <- "figures/"

##protein
protein_ratios <- read.csv(file.path(datadir, "SL_protein/normgated_data_bg_corrected.csv"))
protein <- as.data.frame(protein_ratios)

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

medians <- aggregate( ratio ~ clone + cit + strain, protein, median)
medians <- medians[ medians$strain != "BY4743", ]

medians$elongation_time <- cit[ tolower(medians$cit), 1 ] # look up elongation times by citrine name

WTprotein <- medians[medians$strain == "WT", ]
names(WTprotein) <- tolower( names( WTprotein ))


##mRNA
SL_mRNA_ratios <- read.csv( file.path( datadir, "SL_mRNA/J.WT.mRNA.csv" ))
SLmRNA <- as.data.frame(SL_mRNA_ratios)

WTmRNA <- SLmRNA[SLmRNA$Strain == "WT", ]

WTmRNA$elongation_time <- cit[ tolower(WTmRNA$Cit), 1 ] # look up elongation times by citrine name
names(WTmRNA) <- tolower( names( WTmRNA ))
WTmRNA$cit <- tolower(WTmRNA$cit)


##translation efficiency
te <- merge( WTprotein, WTmRNA, by = c("strain","clone", "cit", "elongation_time"), suffixes = c(".prot", ".mrna") )
te$te <- with( te, ratio / cit_mch)

cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

cair_pdf( file.path(figdir, "efl_citrine_te.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( te$elongation_time, te$te,
      col = cols[te$cit],
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 0.8),
      pch = 20,
      xlab = "",
      ylab = "translation efficiency"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()