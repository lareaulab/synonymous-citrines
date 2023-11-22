library(plyr)
library(tidyverse)

datadir <- "data/individual_CRISPRi/growth_defect_controls/"
figdir <- "figures"

datafile <- file.path( datadir, "crispri_growth_defect_controls_normgated_data.csv" )
data = read.csv( datafile, header=T)

medians <- aggregate( ratio ~ clone + cit + gene + media, data, median) #medians of cit/mcherry ratios for each strain in each condition for each clone
avgs <- aggregate( ratio ~  cit + gene + media, medians, mean ) # average of the median green/red ratios of three clones 

se <- aggregate( ratio ~  cit + gene + media, medians, function(x) sd(x)/sqrt(length(x)) )
avgs$se <- se$ratio

citmins <- avgs[ with( avgs, cit == "ci"), ]
cit9s <- avgs[ with( avgs, cit == "c9"), ]

c9ci_ratios <- as.data.frame(citmins[,c(2:3)])
c9ci_ratios$c9ci <- cit9s$ratio / citmins$ratio
c9ci_ratios$se <- c9ci_ratios$c9ci * sqrt( (with(cit9s, se/ratio)^2 + with(citmins, se/ratio)^2) )

c9ci_matrix <- cbind( c9ci_ratios$c9ci[ with( c9ci_ratios, gene == "CDC6") ], 
                      c9ci_ratios$c9ci[ with( c9ci_ratios, gene == "IRRI") ] )
row.names( c9ci_matrix ) <- c9ci_ratios$media[ with( c9ci_ratios, gene == "CDC6") ]
                      
se_matrix <- cbind( c9ci_ratios$se[ with( c9ci_ratios, gene == "CDC6") ], 
                    c9ci_ratios$se[ with( c9ci_ratios, gene == "IRRI") ] )
row.names( se_matrix ) <- c9ci_ratios$media[ with( c9ci_ratios, gene == "CDC6") ]

mediacols <- c( tet = "#39c0c4", SCD = "#99cccc")
labels <- c( tet = "induced", SCD = "uninduced" )
genes <- c("CDC6", "IRRI")

pdf( file.path( figdir, "growthcontrols.pdf" ), width = 1.5, height = 1.3, pointsize = 6.5, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(4,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
mp = barplot( c9ci_matrix, 
              beside = TRUE,
              col = mediacols[ row.names(c9ci_matrix) ],
              ylab = "slow citrine as\nfraction of fast citrine", 
              xlab = NA,
              border = NA, axes = F,
              ylim = c(0, 0.5) )
axis( 1, lwd = 0.75, labels = genes, at = apply( mp, 2, mean), font = 3, tick = F )
axis( 2, lwd = 0.75 )
arrows( mp, c9ci_matrix - se_matrix, mp, c9ci_matrix + se_matrix, angle=90, code=3, length = 0.02)
legend( 'topright', legend = labels[ row.names(c9ci_matrix) ],
        fill = mediacols[ row.names(c9ci_matrix) ],
        border = NA, bty = "n", inset = c( -0.2, -0.5 ))
dev.off()