library(tidyverse)
library(dplyr)

setwd(".../data/5'standardized/")

gated = read.csv(".../5standardized_normgated_data_bg_corrected.csv", header=T)

citrine_construct_scores_fname <- "figures data and code/citrine_scores_full_model.tsv"

cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")
nm_namemap <- c( cit0 = "cit0", cit3 = "cit3", cit6 = "cit6", cit9 = "cit9", citMax = "citmax", citMin = "citmin")
std_namemap <- c( cit0 = "std0", cit3 = "std3", cit6 = "std6", cit9 = "std9", citMax = "stdmax", citMin = "stdmin")

medians <- aggregate( ratio ~ clone + cit + strain, gated, median)
medians <- medians[ medians$strain != "BY4743", ]

NMmedians <- medians[medians$strain == "NM",]
cit20medians <- medians[medians$strain == "eCit20",]

NMmedians$elongation_time <- cit[ nm_namemap[ NMmedians$cit ], 1 ]
cit20medians$elongation_time <- cit[ std_namemap[ cit20medians$cit], 1 ]

NMmedians <- subset(NMmedians, NMmedians$clone != "A" | NMmedians$cit != "cit9")
NMmedians <- subset(NMmedians, NMmedians$clone != "B" | NMmedians$cit != "cit0")
NMmedians <- subset(NMmedians, NMmedians$clone != "A" | NMmedians$cit != "citMax") # no mCherry signal
cit20medians <- subset(cit20medians, cit20medians$clone != "A" | cit20medians$cit != "cit9")
cit20medians <- subset(cit20medians, cit20medians$clone != "C" | cit20medians$cit != "cit3")
cit20medians <- subset(cit20medians, cit20medians$clone != "A" | cit20medians$cit != "cit3")


#cols <- c( citMin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citMax = "red2")

pdf(".../standardized_citrines.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,2,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( NMmedians$elongation_time, NMmedians$ratio, 
      col = alpha("#3584c6", 1),
      #      col = cols[ NMmedians$cit ],
      pch = 20,
      lwd = 0.75,
      bty = "n",
      xlim = c(150,400),
      ylim = c(0,0.5),
      xlab = "",
      ylab = "citrine / mCherry\nfluorescence ratio"
)
points(cit20medians$elongation_time, cit20medians$ratio, pch = 20, col = alpha("#f58b76", 1))
#points(cit20medians$elongation_time, cit20medians$ratio, pch = 2, col = cols[ cit20medians$cit ] )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

legend( "topright",
        legend = c("original", "5' standardized"),
        col = c("#3584c6", "#f58b76"),
        pch = 20,
        cex=0.7,
        bty = "n",
        inset = c(-0.1, -0.1))

# legend( "topright",
#         legend = c("original", "5' standardized"),
#         #col = black, 
#         pch = c(20, 2),
#         cex=0.7,
#         bty = "n",
#         inset = c(-0.1, -0.1))

dev.off()