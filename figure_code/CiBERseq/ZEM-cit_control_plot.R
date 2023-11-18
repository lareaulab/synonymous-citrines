##plotting the ZEM TF - citrine fusion protein control showing that 
##citmin and cit9 still fluoresce to different degrees when fused 
##to the ZEM transcription factor
library(dplyr)
library(tidyverse)

setwd("../../data/ciBERseq/ZEM-citrine_fusion_control/")
data <- read.csv("ZEM-cit_controls_normgated_data.csv", header = T)

#finding medians
medians <- data %>%
  group_by(cit) %>%
  summarise(median = median(Citrine.cor))
####################################
se <- aggregate( Citrine.cor ~  cit, data, function(x) sd(x)/sqrt(length(x)) )
medians$se <- se$Citrine.cor
row.names(medians) <- medians$cit

#and it must be a matrix:
medians.reordered <- arrange(medians, median)
medians.reordered<-medians.reordered[rownames(medians.reordered)!= "nocit",]

medians_bar <- medians.reordered$median
names(medians_bar) <- medians.reordered$cit

cit.labels <- c("slow cit", "fast cit")

pdf("ZEM-cit_fusioncontrols.pdf", width = 1.75, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
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
axis(2, lwd = 0.75)
arrows( mp, with( medians.reordered, median - se ),
        mp, with( medians.reordered, median + se ), 
        angle=90, code=3, length = 0.02)
dev.off()
