install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("flowStats")
BiocManager::install("flowViz")

library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)

setwd("../../data/RQC_knockouts_and_chimeras/")
PrimaryDirectory <- getwd()

plots_on <- FALSE

#Import OD600 data
tname <- "OD600.xlsx"
tecanimport <- function(tname){
  tmap <- read_excel(tname, sheet = "tecan_map", col_names = TRUE)
  names(tmap)[1] <- "letter"
  #pivot plate map to tidy format
  tmap <- tmap %>%
    pivot_longer(cols = 2:ncol(tmap), names_to = "number") %>%
    transmute(tecan_well = paste(letter, number, sep = ""), name = value)
  #remove empty wells
  tmap <- subset(tmap, is.na(name) == FALSE)
  #parse conditions from names
  for (i in (1:nrow(tmap))) {
    temp_name <- str_split_fixed(tmap[i,2],"-", n = 4)
    final_name <- as_tibble(temp_name, .name_repair = "minimal")
    tmap[i,3] <- final_name[1,1]
    tmap[i,4] <- final_name[1,2]
    tmap[i,5] <- final_name[1,3]
    tmap[i,6] <- final_name[1,4]
  }
  #import OD600 data from Tecan Data
  OD600_data <- read_excel(tname, sheet = "OD600", col_names = FALSE)  
  names(OD600_data)[1] <- "tecan_well"
  names(OD600_data)[2] <- "OD600"
  #remove empty wells
  OD600_data <- subset(OD600_data, is.na(OD600)==FALSE)
  #assign sample names to OD600 values
  tdata <- left_join(OD600_data, tmap, by = "tecan_well")
  tdata_all <<- tdata
  names(tdata_all)[4] <<- "strain"
  names(tdata_all)[5] <<- "clone"
  names(tdata_all)[6] <<- "cit"
}
tecanimport(tname)

#Import and gate flow data - this function runs on all .fcs files in the working directory folder 
#Set variable names for flow data files - must correspond to .fcs file parameters
#Requires "flow_map.xlsx"; column and row 1 indicate well letter and number, rows 2+ are sample names identical to tecan_map
  #flow_map wells must match the wells used on the physical flow plate; indicated in the last 3 characters of each .fcs filename
var.names = c("FSC.A","FSC.H", "SSC.A", "SSC.H", "Citrine.A", "Citrine.H", "mCherry.A", "mCherry.H", "time")
#This paramater changes gating: "The bandwidth factor used for smoothing of the density estimate."
bwfacValue <- 0.5
plots_on <- FALSE
flowimport <- function(itsflowtime){
  #load the fcs files into a set
  flowData <- read.flowSet(
    files = list.files(pattern = "*.fcs"),
    transformation = FALSE
  )
  flowCore::colnames(flowData) <- var.names
  FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
  FileNames <- as.matrix(FileNames)
  sampleNames(flowData) <- FileNames
  truncTrans <- truncateTransform(
    transformationId = "Truncate-transformation",
    a = 1
  )
  myTrans <- transformList(
    c(
      var.names[6],
      var.names[8]
    ),
    truncTrans
  )
  flowData.tt <- flowCore::transform(
    flowData,
    myTrans
  )
  
  #Numerical gating
  c2f <- curv2Filter(
    "FSC.A",
    "SSC.A",
    bwFac = bwfacValue
  )
  c2f.results <- flowCore::filter(
    flowData.tt,
    c2f
  )
  c2f.split <- flowCore::split(
    flowData.tt,
    c2f.results
  )
  foo <- base::Filter(
    function(x) !(nrow(x) < 2),
    mclapply(
      mclapply(
        mclapply(
          c2f.split[2:length(c2f.split)],
          fsApply,
          exprs,
          simplify = FALSE,
          use.exprs = FALSE
        ),
        lapply,
        data.frame
      ),
      ldply
    )
  )
  bar <- ldply(
    Map(
      function(x, y) dplyr::mutate(x, area = y),
      foo,
      names(foo)
    )
  )
  bar$.id <- as.factor(bar$.id)
  bar$area <- as.factor(bar$area)
  #write csv of events by sample area
  write_csv(
    x = bar %>% dplyr::count(.id, area),
    file = "events-by-sample-by-area.csv"
  )
  #gate data; generates a dataframe of all datapoints from all samples within each sample's gate
  gated_data <- ldply(
    Map(
      function(x, y) dplyr::filter(x, area == y),
      dlply(
        bar,
        .(.id)
      ),
      (bar %>%
         dplyr::count(.id, area) %>%
         dplyr::group_by(.id) %>%
         dplyr::filter(n==max(n)))$area
    )
  )
  ########This parser pulls the first 3 characters preceding the ".fcs" in the filename to get the well id
  well_list <- unique(gated_data[,1])
  well_list <- as.matrix(well_list)
  flow_well <- as_tibble(well_list)

  for (i in 1:(nrow(well_list))) {
    temp_name <- str_sub(well_list[i,1], -7, -5)
    flow_well[i,2] <- temp_name[1]
  }
  
  for (i in 1:(nrow(well_list))) {
    if(str_sub(flow_well[i,2], start = -2) != 10) {
      flow_well[i,2] <- gsub('*0*', '', flow_well[i,2])
    }else{
      i <- i+1
    }
  }
  
  names(flow_well)[1] <- ".id"
  names(flow_well)[2] <- "flow_well"
  #join well names to gated flow data
  gated_data <- inner_join(
    gated_data, flow_well,
    by = ".id"
  )
  #Import flow map
  flow_map_file = "../flow_map.csv"
  flow_map <- flow_map_file
  names(flow_map)[1] <- "letter"
  #pivot plate map to tidy format
  flow_map_long <- flow_map %>%
    pivot_longer(cols = 2:ncol(flow_map), names_to = "number") %>%
    transmute(flow_well = paste(letter, number, sep = ""), name = value)
  
  #join tecan data and flow data by flow_well
  flow_data_full <- full_join(
    flow_map_long, gated_data,
    by = "flow_well"
  )
  flow_data_all <<- flow_data_full
  
  #make FSC vs SSC plots
  if(plots_on == TRUE){
  pdf("fsc-ssc-c2f.pdf",
    width = 10,
    height = 30
  )
  print(
    xyplot(
      `SSC.A` ~ `FSC.A`,
      data = flowData.tt,
      smooth = FALSE,
      filter = c2f,
      ylim = c(0,10000),
      xlim = c(0,100000)
    )
  )
  dev.off()
  }
}
flowimport(itsflowtime)

#join flow data and OD600 data by sample name
final_dataset <- full_join(
  flow_data_all, tdata_all,
  by = "name"
)

final_dataset <- na.omit(final_dataset)

write_excel_csv(final_dataset, "gated_data.csv")

datacheck <- final_dataset %>%
  distinct(strain, cit, clone)

final_dataset <- read_csv("gated_data.csv")
final_dataset <- final_dataset[,1:18]

#Plot FSC vs SSC
#load data to plot:

#Fluorescence Analysis
#find background autofluorescence
background <- final_dataset %>%
  dplyr::filter(strain == '4743') %>%
  dplyr::filter(cit == '4743') %>%
  select(Citrine.H) %>%
  summarise(median = median(Citrine.H))

names(background)[1] <- 'citrine'

background[1,2] <- final_dataset %>%
  dplyr::filter(strain == '4743') %>%
  dplyr::filter(cit == '4743') %>%
  select(mCherry.H) %>%
  summarise(median = median(mCherry.H))

names(background)[2] <- 'mCherry'


#subtract background fluorescence
final_dataset_bkgd <- final_dataset %>%
  mutate(Citrine.H = Citrine.H - 70) %>%
  mutate(mCherry.H = mCherry.H - 72)

#add median to flow data
flow_data_ratios <- final_dataset_bkgd %>%
  mutate(f_ratio = Citrine.H/mCherry.H)

#get median f_ratio
f_medians <- distinct(
  flow_data_ratios[c("name", "strain", "clone", "cit")]
)

for (i in 1:nrow(f_medians)) {
  f_medians[i,5] <- flow_data_ratios %>%
    dplyr::filter(name == as.character(f_medians[i,1])) %>%
    summarise(median = median(f_ratio))
}

f_med_raw <- f_medians
write_excel_csv(f_med_chimera, "chimera medians.csv")

f_med_chimera <- f_med_raw %>%
  dplyr::filter(strain == 'Lion' | strain == 'Sphinx' | strain == 'WT')

write_csv(gated_data, "../chimeras_normgated_data_bg_corrected.csv")


