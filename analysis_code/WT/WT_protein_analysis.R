## re-gate and re-analyze the raw flow data from Tunney, McGlincy et al 2018.
## plot the three biological isolates that match the qPCR mRNA measurements
library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

datadir <- "data/WT/WT_protein/"

sample_info <- paste0( datadir, "sample-info.csv" )

citrine_construct_scores_fname <- "data/codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")


var.names <- c("FSC.A", "SSC.A", "Citrine.A", "mCherry.A", "time") # FITC, PE-Texas Red
plots_on = F

# Import and gate flow data - run on all .fcs files in directory
# load the fcs files into a set (a list with many flowFrames, each from a different .fcs experiment file)
flowData <- read.flowSet( files = list.files( datadir, pattern = "_[CEF]0[1345678].*fcs", full.names = T ), transformation = F)
# limit this to the isolates with qPCR data -- isolates 3, 5, 6: rows C, E F and columns 1,3:8

flowCore::colnames(flowData) <- var.names

sampleNames(flowData) <- sapply( sampleNames(flowData),
                                 function(x){ substr( x, nchar(x) - 13, nchar(x) - 12 ) } )

# gating: fits 2d normal distribution, selects all events within the Mahalanobis distance
# the sampling process to estimate covariance may explain why we get a very slightly different number of filtered events
# each time it's run on the same data.
normfilter <- norm2Filter( "FSC.A", "SSC.A") 
normresults <- flowCore::filter( flowData, normfilter )
flowDataGated <- flowCore::Subset( flowData, normresults )

# diagnostic plots of flow data
if(plots_on == TRUE) {
   # plot all FSC vs SSC data and superimpose gate boundary
   png( paste0( datadir, "fsc-ssc-normgates.png" ), width = 20, height = 10, units = "in", res = 300)
   print( xyplot( `SSC.A` ~ `FSC.A`, data = flowData, smooth = FALSE, filter = normresults), layout = c(7,3))
   dev.off()
   # plot gated FSC vs SSC data only
   png( paste0( datadir, "fsc-ssc-normgated-only.png" ), width = 20, height = 10, units = "in", res = 300 )
   print( xyplot( `SSC.A` ~ `FSC.A`, data = flowDataGated, smooth = FALSE, layout = c(7,3)))
   dev.off()
   # plot gated fluorescence data only
   png( paste0( datadir, "fluor-normgated-only.png" ), width = 20, height = 10, units = "in", res = 300 )
   print( xyplot( `Citrine.A` ~ `mCherry.A`, data = flowDataGated, smooth = FALSE, 
                  ylim = c(0,3000), xlim = c(0,10000) ), layout = c(7,3))
   dev.off()
   # plot all fluorescence data
   png( paste0( datadir, "fluor-ungated.png" ), width = 20, height = 10, units = "in", res = 300 )
   print( xyplot( `Citrine.A` ~ `mCherry.A`, data = flowData, smooth = FALSE,
          ylim = c(0,3000), xlim = c(0,10000), layout = c(7,3)))
   dev.off()
}

# convert the flowSet data to matrices then dataframes and combine them into one big dataframe
gated_data <- fsApply( flowDataGated, exprs, simplify = F, use.exprs = F )
gated_data <- lapply( gated_data, data.frame )
gated_data <- ldply( gated_data )

# map well to sample
flow_map <- read.csv( sample_info, header = T )
row.names(flow_map) <- flow_map$Well

gated_data <- cbind( gated_data,  flow_map[ gated_data$.id, ])

background <- gated_data %>%
   dplyr::filter(Strain == 'NIY110') %>%
   select(Citrine.A) %>%
   summarise(median = median(Citrine.A))

mchbackground <- gated_data %>%
   dplyr::filter(Strain == 'NIY110') %>%
   select(mCherry.A) %>%
   summarise(median = median(mCherry.A))


gated_data$Citrine.cor <- gated_data$Citrine.A - background$median
gated_data$mCherry.cor <- gated_data$mCherry.A - mchbackground$median

gated_data$ratio <- with( gated_data, Citrine.cor / mCherry.cor)

ratios <- aggregate( ratio ~ Isolate + Strain, gated_data, median)

# we are only using 3,5,6 so this doesn't matter here, but see original paper scripts for notes on the other isolates - a few had FACS failures 
# rep 3 of cit999 was later discarded after testing genomic copy number 
ratios$ratio[ with( ratios, Isolate == 3 & Strain == "999" ) ] <- NA

#means <- aggregate( ratio ~ Strain, ratios, mean)

# remove control strain
ratios <- ratios[ ratios$Strain != "NIY110", ]

# remap names 
rename <- c( "000" = "cit0", "333" = "cit3", "666" = "cit6", "999" = "cit9", "MAX" = "citmax", "MIN" = "citmin")
ratios$Strain <- rename[ratios$Strain]
ratios <- within( ratios, citscore <- cit[Strain,] )

write.csv(ratios, paste0( datadir, "WT_protein_ratios.csv" ), row.names=FALSE)
