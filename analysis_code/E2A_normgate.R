##E2A reporter flow cytometry file analysis

library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)


#Load Data
datadir <- "data/E2A/"
flow_map_file = file.path( datadir, "flow_map.csv")

PrimaryDirectory <- getwd()
plots_on <- FALSE

#Import and gate flow data - run on all .fcs files in the working directory folder 
#Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("FSC.A","FSC.H", "SSC.A", "SSC.H", "Citrine.A", "Citrine.H", 
               "mCherry.A", "mCherry.H", "iRFP.A", "iRFP.H", "time")

#load the fcs files into a set
flowData <- read.flowSet( files = list.files(path = datadir, pattern = "*.fcs", full.names = T), transformation = F )
#list with 57 flowFrames, each from a different .fcs experiment file

flowCore::colnames(flowData) <- var.names
sampleNames(flowData) <- sapply( list.files( path = datadir, pattern = ".fcs", full.names = T),
                                 function(x){ substr( x, nchar(x) - 6, nchar(x) - 4 ) })

# gating: fits 2d normal distribution, selects all events within the Mahalanobis distance
normfilter <- norm2Filter( "FSC.A", "SSC.A") 
normresults <- flowCore::filter( flowData, normfilter )
flowDataGated <- flowCore::Subset( flowData, normresults )

# diagnostic plots of flow data
if(plots_on == TRUE) {
  #plot data with gates
  png( file.path( datadir, "fsc-ssc-normgates.png" ), width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `SSC.A` ~ `FSC.A`, data = flowData, smooth = FALSE, filter = normresults, ylim = c(0,20000), xlim = c(0,100000) ))
  dev.off()
  # plot gated data only
  png( file.path( datadir, "fsc-ssc-normgated-only.png" ), width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `SSC.A` ~ `FSC.A`, data = flowDataGated, smooth = FALSE, ylim = c(0,20000), xlim = c(0,100000) ))
  dev.off()
  # plot gated fluorescence data only
  png( file.path( datadir, "fluor-normgated-only.png" ), width = 20, height = 20, units = "in", res = 300 )
  print( xyplot( `Citrine.H` ~ `mCherry.H`, data = flowDataGated, smooth = FALSE, ylim = c(0,2000), xlim = c(0,10000) ))
  dev.off()
  png( file.path( datadir, "fluor-all.png" ), width = 20, height = 20, units = "in", res = 300 )
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
citbackground <- gated_data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(Citrine.A) %>%
  summarise(median = median(Citrine.A))

mchbackground <- gated_data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(mCherry.A) %>%
  summarise(median = median(mCherry.A))

iRFPbackground <- gated_data %>%
  dplyr::filter(strain == 'BY4743') %>%
  select(iRFP.A) %>%
  summarise(median = median(iRFP.A))

#background correction
gated_data$Citrine.cor <- gated_data$Citrine.A - citbackground$median
gated_data$mCherry.cor <- gated_data$mCherry.A - mchbackground$median
gated_data$iRFP.cor <- gated_data$iRFP.A - iRFPbackground$median

# ratio of citrine to mCherry for each data point
gated_data$cit.iRFP.ratio.A <- with( gated_data, Citrine.cor / iRFP.cor)
gated_data$mCh.iRFP.ratio.A <- with( gated_data, mCherry.cor / iRFP.cor)
gated_data$mCh.cit.ratio.A <- with( gated_data, mCh.iRFP.ratio.A / cit.iRFP.ratio.A)

gated_data$cit.mCh.ratio.A.raw <- with( gated_data, Citrine.cor / mCherry.cor)

gated_data <- na.omit(gated_data)

write_csv(gated_data, file.path( datadir, "E2A_normgated_data_corrected.csv" ))

# ci_ratio <- gated_data %>%
#   filter(cit == "citMin") %>%
#   select(mCh.cit.ratio.A) %>%
#   summarise(median = median(mCh.cit.ratio.A))
# 
# c9_ratio <- gated_data %>%
#   filter(cit == "cit9") %>%
#   select(mCh.cit.ratio.A) %>%
#   summarise(median = median(mCh.cit.ratio.A))
# 
# all_data_medians <- gated_data %>%
#   group_by(cit, clone) %>%
#   filter(strain != "BY4743") %>%
#   summarise(mCh.cit.med = median(mCh.cit.ratio.A),
#             cit.med = median(cit.iRFP.ratio.A),
#             mCh.med = median(mCh.iRFP.ratio.A),
#             cit.raw = median(Citrine.cor),
#             mCh.raw = median(mCherry.cor),
#             cit.mCh.raw = median(cit.mCh.ratio.A.raw))
# 
# write_csv(all_data_medians, "E2A_medians.csv")
# 
# 
