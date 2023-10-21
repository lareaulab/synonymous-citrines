library(plyr)
library(tidyverse)

setwd("../../data/individual_CRISPRi/growth_defect_controls/")

datafile <- ("crispri_growth_defect_controls_normgated_data.csv")
data = read.csv( datafile, header=T)

medians <- aggregate( ratio ~ clone + cit + gene + media, data, median) #medians of cit/mcherry ratios for each strain in each condition for each clone
avgs <- aggregate( ratio ~  cit + gene + media, medians, mean ) # average of the median green/red ratios of three clones 

se <- aggregate( ratio ~  cit + gene + media, medians, function(x) sd(x)/sqrt(length(x)) )
avgs$se <- se$ratio

citmins <- avgs[ with( avgs, cit == "ci"), ]
cit9s <- avgs[ with( avgs, cit == "c9"), ]


citminstet <- citmins[ which(citmins$media=='tet') ]
citminstet <- citmins[with( citmins, media == "tet")]
cit9stet <- cit9s[with( cit9s, media == "tet")]


c9ci_ratios <- as.data.frame(citmins[,c(2:3)])
c9ci_ratios$c9ci <- cit9s$ratio / citmins$ratio
c9ci_ratios$se <- c9ci_ratios$c9ci * sqrt( (with(cit9s, se/ratio)^2 + with(citmins, se/ratio)^2) )

c9ci_ratios_tet <- c9ci_ratios[c9ci_ratios$media != 'tet', ]  
c9ci_ratios_scd <- c9ci_ratios[c9ci_ratios$media != 'SCD', ]  

#c9ci_ratios$alias[ c9ci_ratios$gene == "HO" ] <- "control"
#c9ci_ratios$alias <-  paste0( "(", c9ci_ratios$alias, ")" )

row.names(c9ci_ratios) <- c9ci_ratios$gene

#push the data into the appropriate format for bar graphs
c9ci_ratios_only <- subset(c9ci_ratios, select = -se )
c9ci_ratios_wide <- c9ci_ratios_only %>%
  pivot_wider(
    names_from = c(gene),
    values_from = c(c9ci)
  )
c9ci_ratios_wide <- as.data.frame(c9ci_ratios_wide)
#columns must have just numbers, no characters:
rownames(c9ci_ratios_wide) <- c9ci_ratios_wide[,1]
c9ci_ratios_wide[,1] <- NULL
#and it must be a matrix:
c9ci_ratios_wide.matrix <- as.matrix(c9ci_ratios_wide)

c9ci_ratios_reordered <- c9ci_ratios %>% arrange(gene)

mediacols <- c("#39c0c4", "#99cccc")
gene.labels <- c("CDC6", "IRRI")

pdf("../../../figures/c9ciratio_induceduninduced_growthcontrols.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
#set margin stuff
par( mex = 0.65 ) 
par( mar = c(6,6.5,5,9) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

mp = barplot(c9ci_ratios_wide.matrix,
             beside = TRUE,
             col = mediacols,
             ylab = "slow/fast citrine fluorescence",
             border = NA, axes = F,
             ylim = c(0, 0.5),
             names.arg = gene.labels)

axis(2)

arrows( mp, with( c9ci_ratios_reordered, c9ci - se ), mp, with( c9ci_ratios_reordered, c9ci + se ), angle=90, code=3, length = 0.02)

legend( 'topright', legend = c("induced", "uninduced"), 
        fill = mediacols, border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()
