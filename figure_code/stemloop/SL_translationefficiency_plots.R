##Calculating and plotting the translation efficiency (protein/mRNA) for 
##each SL-citrine construct.

library(dplyr)
library(plyr)
library(tidyr)

datadir <- "data/stemloop/"
figdir <- "figures"

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")

#qPCR (mRNA) data:
mRNA <- read.csv( file.path( datadir, "SL_mRNA/J.WT.mRNA.csv" ), header = T)
colnames(mRNA)[which(colnames(mRNA) == "Cit_mCh")] <- "normalized.mrna"
mRNA$Cit <- tolower( mRNA$Cit )
colnames(mRNA) <- tolower( colnames(mRNA) )

wt.mrna.avg <- aggregate( normalized.mrna ~ cit + strain, mRNA[mRNA$strain == "WT",], mean)
j.mrna.avg <- aggregate( normalized.mrna ~ cit + strain, mRNA[mRNA$strain == "HC1J",], mean)

#fluorescence (protein) data:
gated = read.csv( file.path( datadir, "SL_protein/normgated_data_bg_corrected.csv" ), header=T, row.names=1)

fluorescence <- aggregate( ratio ~ clone + cit + strain, gated, median)
fluorescence <- fluorescence[ fluorescence$strain != "BY4743", ]
colnames(fluorescence)[which(colnames(fluorescence) == "ratio")] <- "normalized.protein"
fluorescence$strain = toupper(fluorescence$strain)

fluorescence$elongation_time <- cit[ tolower(fluorescence$cit), 1 ] # look up elongation times by citrine name

wt.avg = aggregate( normalized.protein ~ cit + strain + elongation_time, fluorescence[fluorescence$strain == "WT",], mean)
j.avg = aggregate( normalized.protein ~ cit + strain + elongation_time, fluorescence[fluorescence$strain == "HC1J",], mean)
g.avg = aggregate( normalized.protein ~ cit + strain + elongation_time, fluorescence[fluorescence$strain == "HC1G",], mean)

translation <- merge(fluorescence, mRNA, by = c("strain", "cit", "clone"))
translation$te <- translation$normalized.protein / translation$normalized.mrna

wt.te <- merge( wt.avg, wt.mrna.avg, by = c("strain", "cit") )
wt.te$te <- with( wt.te, normalized.protein / normalized.mrna )
wt.te <- wt.te[ order( wt.te$elongation_time ), ]

j.te <- merge( j.avg, j.mrna.avg, by = c("strain", "cit") )
j.te$te <- with( j.te, normalized.protein / normalized.mrna )
j.te <- j.te[ order( j.te$elongation_time ),]

cols <- c( WT = "#3584c6", HC1G = "#fbb615", HC1J = "#e95c64" )

#main figure TE plot
pdf( file.path( figdir, "stemloop_te.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( translation$elongation_time, translation$te,
      col = cols[ translation$strain ],
      pch = 20,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 0.8),
      xlab = NA,
      ylab = "translation efficiency"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( wt.te$elongation_time, wt.te$te, lwd = 1.5, col = cols[ "WT" ])
lines( j.te$elongation_time, j.te$te, lwd = 1.5, col = cols[ "HC1J" ])

legend( "topright", pch = 20, inset = c( -0.2, -0.6 ), 
        legend = c("no stem loop", "strong stem loop"), 
        col = cols[ c("WT", "HC1J") ], 
        cex = 0.8,
        bty = "n")
dev.off()


## slightly different results from these two approaches? check this

#Plotting ratio of SL to WT TE for each citrine 

# #average all of the biological replicates (A-C) of each type of sample (strain-cit)
# BioRepAvg <- aggregate( Cit_mCh ~ Strain + Cit + elongation_time, ratios, mean)
# 
# #dividing hairpin by wt
# divided <- BioRepAvg %>%
#   pivot_wider(
#     names_from = Strain,
#     values_from = c(Cit_mCh)
#   )
# divided$ratio <- divided$HC1J/divided$WT
# divided$elongation_time <- cit[ tolower( divided$Cit), 1 ]

te.ratio <- merge( wt.te, j.te, by = c("cit"), suffixes = c(".wt", ".j") )
te.ratio$ratio <- with( te.ratio, te.j / te.wt )

#plot:
pdf( file.path( figdir, "stemloop_te_ratios.pdf" ), width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( te.ratio$elongation_time.wt, te.ratio$ratio,
      col = "darkgrey",
      pch = 20,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 1),
      xlab = "",
      ylab = "SL / no SL\nTE ratio"
)
axis( 1, lwd = 0.75, at = seq(150, 400, by =50), labels = c(NA, "200", NA, "300", NA, "400") )
axis( 2, lwd = 0.75 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
dev.off()






