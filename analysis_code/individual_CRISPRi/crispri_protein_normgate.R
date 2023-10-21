#load libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("flowCore")
# BiocManager::install("flowStats")
# BiocManager::install("flowViz")

library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

plots_on <- F

setwd(".../data/individual_CRISPRi/CRISPRi_protein/FCS_files/")

guide_key <- read.csv("../guide_key.csv", row.names = 1, header = F, col.names = c("strain","gene","alias"))

# Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("FSC.A","FSC.H", "SSC.A", "SSC.H", "Citrine.A", "Citrine.H", "mCherry.A", "mCherry.H", "time")

# Import and gate flow data - run on all .fcs files in directory
flowimport <- function( directory ){

  setwd(directory)
  # load the fcs files into a set (a list with many flowFrames, each from a different .fcs experiment file)
  flowData <- read.flowSet( files = list.files( pattern = "*.fcs"), transformation = F )

  flowCore::colnames(flowData) <- var.names
  
  sampleNames(flowData) <- sapply( sampleNames(flowData),
                                   function(x){ substr( x, nchar(x) - 6, nchar(x) - 4 ) } )
  
  # gating: fits 2d normal distribution, selects all events within the Mahalanobis distance
  # the sampling process to estimate covariance may explain why we get a very slightly different number of filtered events
  # each time it's run on the same data.
  normfilter <- norm2Filter( "FSC.A", "SSC.A") 
  normresults <- flowCore::filter( flowData, normfilter )
  flowDataGated <- flowCore::Subset( flowData, normresults )
  
  # diagnostic plots of flow data
  if(plots_on == TRUE) {
    # plot all FSC vs SSC data and superimpose gate boundary
    png( "fsc-ssc-normgates.png", width = 20, height = 20, units = "in", res = 300 )
    print( xyplot( `SSC.A` ~ `FSC.A`, data = flowData, smooth = FALSE, filter = normresults, ylim = c(0,20000), xlim = c(0,100000) ))
    dev.off()
    # plot gated FSC vs SSC data only
    png( "fsc-ssc-normgated-only.png", width = 20, height = 20, units = "in", res = 300 )
    print( xyplot( `SSC.A` ~ `FSC.A`, data = flowDataGated, smooth = FALSE, ylim = c(0,20000), xlim = c(0,100000) ))
    dev.off()
    # plot gated fluorescence data only
    png( "fluor-normgated-only.png", width = 20, height = 20, units = "in", res = 300 )
    print( xyplot( `Citrine.H` ~ `mCherry.H`, data = flowDataGated, smooth = FALSE, ylim = c(0,2000), xlim = c(0,10000) ))
    dev.off()
    # plot all fluorescence data
    png("fluor-ungated.png", width = 20, height = 20, units = "in", res = 300 )
    print( xyplot( `Citrine.H` ~ `mCherry.H`, data = flowData, smooth = FALSE, ylim = c(0,2000), xlim = c(0,10000) ))
    dev.off()
  }
  
  # convert the flowSet data to matrices then dataframes and combine them into one big dataframe
  gated_data <- fsApply( flowDataGated, exprs, simplify = F, use.exprs = F )
  gated_data <- lapply( gated_data, data.frame )
  gated_data <- ldply( gated_data )
  
  # map well to sample
  flow_map <- read.csv( "flow_map.csv", header = T, row.names = 1)
  colnames(flow_map) <- sprintf("%02d", 1:ncol(flow_map)) 
  wells <- apply( expand.grid(row.names(flow_map), colnames(flow_map)), 1, function(x){ paste0(x[1], x[2]) })
  
  flow_map <- data.frame( sample = c(as.matrix(flow_map)), row.names = wells )

  # to avoid problems with empty wells. might not work great if the legit well values don't all contain -
  split_line <- function( line ) {
    if ( is.na(line) | line == "" ) { NA }
    else { strsplit( line, "-") }
  }
  splitcols <- do.call(rbind, sapply( flow_map$sample, split_line ))
  colnames(splitcols) <-  c( "strain", "clone", "cit", "media" )
  
  flow_map <- cbind( flow_map, splitcols )
  flow_map <- cbind( flow_map, guide_key[flow_map$strain,] )
  
  gated_data <- cbind( gated_data,  flow_map[ gated_data$.id, ])

  setwd("..")
  return(gated_data)
}

# list data folders from different dates (10_28_21, 11_1_21, 11_24_21, 12_3_21)
dirs <- list.files(getwd(), pattern="_21")

# run flowimport on all data!
data <- do.call(rbind, lapply( dirs, flowimport ))

data <- na.omit(data)

#find WT background
background <- data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(Citrine.H) %>%
  summarise(median = median(Citrine.H))

mchbackground <- data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(mCherry.H) %>%
  summarise(median = median(mCherry.H))

#background correction
data$Citrine.cor <- data$Citrine.H - background$median
data$mCherry.cor <- data$mCherry.H - mchbackground$median

data$ratio <- data$Citrine.cor / data$mCherry.cor

write_csv(data, "../crispri_normgated_data_bg_corrected.csv")
