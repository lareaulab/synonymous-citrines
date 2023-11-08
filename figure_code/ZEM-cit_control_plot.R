##plotting the ZEM TF - citrine fusion protein control showing that 
##citmin and cit9 still fluoresce to different degrees when fused 
##to the ZEM transcription factor
library(plyr)
library(tidyverse)

setwd("~/data/ciBERseq/ZEM-citrine_fusion_control")
data <- read.csv("ZEM-cit_controls_normgated_data.csv", header = T)

#finding medians
medians <- data %>%
  group_by(cit) %>%
  summarise(median = median(Citrine.A))

se <- aggregate( Citrine.A ~  cit, data = data, function(x) sd(x)/sqrt(length(x)) )
medians$se <- se$Citrine.A
row.names(medians) <- medians$cit
medians.reordered <- arrange(medians, median)
medians_bar <- medians.reordered$median
names(medians_bar) <- medians.reordered$cit


cit.labels <- c("no citrine", "slow citrine", "fast citrine")

pdf("ZEM-cit_fusioncontrols.pdf", width = 3.34646 , height = 2.79429, pointsize = 7, useDingbats = F, bg = "white" )
#set margin stuff
par( mex = 0.65 ) 
par( mar = c(6,6.5,5,9) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

mp = barplot(medians_bar,
             beside = TRUE,
             col = "darkgrey",
             ylab = "citrine fluorescence",
             border = NA, axes = F,
             ylim = c(0, 450),
             names.arg = cit.labels)


axis(2)

arrows( mp, with( medians.reordered, median - se ), mp, with( medians.reordered, median + se ), angle=90, code=3, length = 0.02)
dev.off()
