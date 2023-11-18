if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("flowCore")
#BiocManager::install("flowStats")
#BiocManager::install("flowViz")

library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

setwd("/Users/loudevanneaux/Desktop/Citrine Codons Paper/Github3/Citrine-Codons-Data/data/5p_standardized/")
flow_map_file = "//Users/loudevanneaux/Desktop/Citrine Codons Paper/Github3/Citrine-Codons-Data/data/5p_standardized/flow_map.csv"

PrimaryDirectory <- getwd()
plots_on <- FALSE

#Import and gate flow data - run on all .fcs files in the working directory folder 
#Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("FSC.A","FSC.H", "SSC.A", "SSC.H", "Citrine.A", "Citrine.H", "mCherry.A", "mCherry.H", "time")


#load the fcs files into a set
flowData <- read.flowSet( files = list.files(path = PrimaryDirectory, pattern = "*.fcs"), transformation = F )
flowCore::colnames(flowData) <- var.names
sampleNames(flowData) <- sapply( list.files( path = PrimaryDirectory, pattern = ".fcs" ),
                                 function(x){ substr( x, nchar(x) - 6, nchar(x) - 4 ) })

# gating: fits 2d normal distribution, selects all events within the Mahalanobis distance
normfilter <- norm2Filter( "FSC.A", "SSC.A") 
normresults <- flowCore::filter( flowData, normfilter )
flowDataGated <- flowCore::Subset( flowData, normresults )

# diagnostic plots of flow data
if(plots_on == TRUE) {
  #plot data with gates
  png( "../fsc-ssc-normgates.png", width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `SSC.A` ~ `FSC.A`, data = flowData, smooth = FALSE, filter = normresults, ylim = c(0,20000), xlim = c(0,100000) ))
  dev.off()
  # plot gated data only
  png( "../fsc-ssc-normgated-only.png", width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `SSC.A` ~ `FSC.A`, data = flowDataGated, smooth = FALSE, ylim = c(0,20000), xlim = c(0,100000) ))
  dev.off()
  # plot gated fluorescence data only
  png("../fluor-normgated-only.png", width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `Citrine.H` ~ `mCherry.H`, data = flowDataGated, smooth = FALSE, ylim = c(0,2000), xlim = c(0,10000) ))
  dev.off()
  png("../fluor-all.png", width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `Citrine.H` ~ `mCherry.H`, data = flowData, smooth = FALSE, ylim = c(0,2000), xlim = c(0,10000) ))
  dev.off()
}

# convert the flowSet data to matrices then dataframes and combine them into one big dataframe
gated_data <- fsApply( flowDataGated, exprs, simplify = F, use.exprs = F)
gated_data <- lapply( gated_data, data.frame )
gated_data <- ldply( gated_data )

###
# map sample names to wells
flow_map <- read.csv( flow_map_file, header = T, row.names = 1)
colnames(flow_map) <- sprintf("%02d", 1:12) 
wells <- apply( expand.grid(row.names(flow_map), colnames(flow_map)), 1, function(x){ paste0(x[1], x[2]) })
flow_map <- data.frame( sample = c(as.matrix(flow_map)), row.names = wells )
#flow_map <- data.frame( .id = wells, sample = c(as.matrix(flow_map)))
#flow_map <- flow_map[flow_map$sample != "",]

split_line <- function( line ) {
  if ( is.na(line) | line == "" ) { NA }
  else { strsplit( line, "-") }
}
splitcols <- do.call(rbind, sapply( flow_map$sample, split_line ))
colnames(splitcols) <-  c( "strain", "cit", "clone" )

flow_map <- cbind( flow_map, splitcols )

gated_data <- cbind( gated_data,  flow_map[ gated_data$.id, ])

gated_data <- na.omit(gated_data)

#find WT background
background <- gated_data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(Citrine.H) %>%
  summarise(median = median(Citrine.H))

mchbackground <- gated_data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(mCherry.H) %>%
  summarise(median = median(mCherry.H))

#background correction
gated_data$Citrine.cor <- gated_data$Citrine.H - background$median
gated_data$mCherry.cor <- gated_data$mCherry.H - mchbackground$median

# ratio of citrine to mCherry for each data point
gated_data$ratio <- with( gated_data, Citrine.cor / mCherry.cor)
gated_data$ratio.H <- with( gated_data, Citrine.H / mCherry.H)
#gated_data$ratio.A <- with( gated_data, Citrine.A / mCherry.A)
#gated_data$FSCratio.A <- with( gated_data, Citrine.A / FSC.A)

write_csv(gated_data, "5standardized_normgated_data_bg_uncorrected.csv")
