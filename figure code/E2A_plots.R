library(tidyverse)
library(dplyr)

setwd("..../data/E2A/")

datafile <- (".../E2A_normgated_data_corrected.csv")
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


#colsdark <- c("#f26739", "#983794")
colsdark <- c("darkorange2", "magenta3")
colslight <- c("#ffcc66", "#cb9ac6")

all.matrix <- as.matrix( medians.avg[,2:4] )
all.labels <- c( "cit/iRFP", "mCh/iRFP", "cit/mCh" )

ind.matrix <- as.matrix( medians.avg[,2:3] )
ind.labels <- c( "cit/iRFP", "mCh/iRFP")

scaled.matrix <- as.matrix( medians.avg.scaled[, 2:4])

pdf("e2a_scaled.pdf", width = 2, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(4,6.5,4,5) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( scaled.matrix, 
              beside = TRUE,
              col = colslight,
              ylab = "fluorescence ratio\n(scaled)",
              border = NA, axes = F,
              #ylim = c(0, 2.5),
              names.arg = all.labels)
#axis(2, seq(0,2.5,0.5))
axis(2, lwd = 0.75)
points( x = rep(as.numeric(mp), each = 3), 
        y = c( medians.scaled$cit.iRFP, medians.scaled$mCh.iRFP, medians.scaled$cit.mCh ),
        pch = 20, col = rep( colsdark, each = 3 ))
legend( "topright", legend = c("slow citrine", "fast citrine"), 
         fill = colslight, border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()

pdf("e2a_all.pdf", width = 2, height = 1.2, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(4,6.5,4,5) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( all.matrix, 
              beside = TRUE,
              col = colslight,
              ylab = "fluorescence ratio",
              border = NA, axes = F,
              #ylim = c(0, 2.5),
              names.arg = all.labels)
#axis(2, seq(0,2.5,0.5))
axis(2)
points( x = rep(as.numeric(mp), each = 3), 
        y = c( medians$cit.iRFP, medians$mCh.iRFP, medians$cit.mCh ),
        pch = 20, col = rep( colsdark, each = 3 ))
legend( "topright", legend = c("slow citrine", "fast citrine"), 
        fill = colslight, border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()


pdf("e2a_cit_mch.pdf", width = 1.5, height = 1, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(3,6.5,2,4) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( ind.matrix, 
              beside = TRUE,
              col = colslight,
              ylab = "fluorescence ratio",
              border = NA, axes = F,
              ylim = c(0, 2.5),
              names.arg = ind.labels)
axis(2, seq(0,2.5,0.5))
points( x = rep(as.numeric(mp), each = 3), 
        y = c( medians$cit.iRFP, medians$mCh.iRFP ),
        pch = 20, col = rep( colsdark, each = 3 ))
# legend( x=6, y=2, legend = c("slow citrine", "fast citrine"), 
#         fill = colslight, border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()

  
  