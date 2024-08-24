# new citrine constructs for this paper - mRNA measurements with qPCR

#set up the directories
datadir <- "data/stemloop/SL_mRNA/" 
figdir <- "figures/"

#read in data from stemloop experiments
SL_mRNA_ratios <- read.csv( file.path( datadir, "J.WT.mRNA.csv" ))
SLratios <- as.data.frame(SL_mRNA_ratios)
#select only WT data (EFL yeast)
WTratios <- SLratios[SLratios$Strain == "WT", ]

#add elongation times:
citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")
WTratios$elongation_time <- cit[ tolower(WTratios$Cit), 1 ] # look up elongation times by citrine name


#mRNA levels plot
cols <- c( citMin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citMax = "red2")

cairo_pdf( file.path( figdir, "efl_citrine_mrna.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( WTratios$elongation_time, WTratios$Cit_mCh, 
      col = cols[WTratios$Cit],
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 1),
      pch = 20,
      xlab = "",
      ylab = "citrine / mCherry\nrelative mRNA ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()