library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

datadir <- "data/stemloop/SL_protein/"
fcsdir <- file.path( datadir, "gate_example" )
figdir <- "figures"

flowData <- read.flowSet( files = list.files(path = fcsdir, pattern = "*.fcs", full.names = T), transformation = F )

#Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC-W", "Citrine.A", "mCherry.A", "time")
flowCore::colnames(flowData) <- var.names

sampleNames(flowData) <- gsub(".fcs", "", sampleNames(flowData) )

normfilter <- norm2Filter( "FSC.A", "SSC.A") 
normresults <- flowCore::filter( flowData, normfilter )
flowDataGated <- flowCore::Subset( flowData, normresults )

# 50k events measured per sample; approx 20k events remain per sample after gating

#plot forward and side scatter and draw gate
png( file.path( figdir, "gate_example_fsc_ssc.png" ), width = 6, height = 6, units = "in", res = 300 )
print( xyplot( SSC.A ~ FSC.A, data = flowData, smooth = FALSE, filter = normresults, ylim = c(0,100000), xlim = c(0,200000) ))
dev.off()

# plot gated fluorescence data 
png( file.path( figdir, "gated_fluor_data.png" ), width = 6, height = 6, units = "in", res = 300 )
print( xyplot( Citrine.A ~ mCherry.A, data = flowDataGated, smooth = FALSE, ylim = c(0,2000), xlim = c(0,6000) ))
dev.off()
