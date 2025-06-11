library(tidyverse)
library(dplyr)

datadir <- "data/E2A/"
figdir <- "figures/"

datafile <- file.path( datadir, "E2A_normgated_data_corrected.csv")
gated = read.csv( datafile, header=T)

#verifying that the cit/mch is actually divided that way:
gated$mCh.cit.ratio.A <- gated$Citrine.cor/gated$mCherry.cor
colnames(gated)[colnames(gated) == "mCh.cit.ratio.A"] ="cit.mCh"
colnames(gated)[colnames(gated) == "cit.iRFP.ratio.A"] ="cit.iRFP"
colnames(gated)[colnames(gated) == "mCh.iRFP.ratio.A"] ="mCh.iRFP"


# for each clone, get the median cit/iRFP, median mch/iRFP, and median cit/mch
citmedians <- aggregate( cit.iRFP ~ clone + cit + strain, gated, median)
mChmedians <- aggregate( mCh.iRFP ~ clone + cit + strain, gated, median)
ratiomedians <- aggregate( cit.mCh ~ clone + cit + strain, gated, median)

medians <- merge(citmedians, mChmedians, sort = F)
medians <- merge(medians, ratiomedians, sort = F)

# remove the no-fluorescence control
medians <- medians[ medians$strain != "BY4743", ]
# strain column is now extraneous
medians$strain <- NULL

#making a table with the averages of all clones for each result for each citrine for the bar graph y values:
citmedians.avg <- aggregate( cit.iRFP ~ cit, medians, mean)
mChmedians.avg <-  aggregate( mCh.iRFP ~ cit, medians, mean)
ratiomedians.avg <-  aggregate( cit.mCh ~ cit, medians, mean)

medians.avg <- merge(citmedians.avg, mChmedians.avg)
medians.avg <- merge(medians.avg, ratiomedians.avg)

medians.avg.scaled <- medians.avg
medians.avg.scaled$cit.iRFP = medians.avg$cit.iRFP / medians.avg$cit.iRFP[2]
medians.avg.scaled$mCh.iRFP = medians.avg$mCh.iRFP / medians.avg$mCh.iRFP[2]
medians.avg.scaled$cit.mCh = medians.avg$cit.mCh/ medians.avg$cit.mCh[2]

medians.scaled <- medians
medians.scaled$cit.iRFP = medians$cit.iRFP / medians.avg$cit.iRFP[2]
medians.scaled$mCh.iRFP = medians$mCh.iRFP / medians.avg$mCh.iRFP[2]
medians.scaled$cit.mCh = medians$cit.mCh/ medians.avg$cit.mCh[2]


#colsdark <- c("darkorange2", "magenta3")
#colslight <- c("#ffcc66", "#cb9ac6")

colsdark <- c( cit9 = "#365C8DFF", citmin = "#9FDA3AFF" )
colslight <- c( cit9 = "#6699CCFF", citmin = "#CCFF99FF" )

all.labels <- c( "cit/iRFP", "mCh/iRFP", "cit/mCh" )

scaled.matrix <- as.matrix( medians.avg.scaled[, 2:4])

cairo_pdf( file.path( figdir, "e2a.pdf"), width = 2, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(3,6.5,5,5) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( scaled.matrix, 
              beside = TRUE,
              col = colslight,
              ylab = "fluorescence ratio\n(scaled)",
              ylim = c(0, 1),
              border = NA, axes = F,
              names.arg = all.labels, 
              mgp = c(3, 0, 0) )#, padj = 1)
axis(2, at = c(0, 0.5, 1), lwd = 0.75)
points( x = rep(as.numeric(mp), each = 3), 
        y = c( medians.scaled$cit.iRFP, medians.scaled$mCh.iRFP, medians.scaled$cit.mCh ),
        pch = 20, col = rep( colsdark, each = 3 ))
legend( "topright", legend = c("slow citrine", "fast citrine"), 
         fill = colslight, border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()
