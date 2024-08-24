##plotting the ZEM TF - citrine fusion protein control showing that 
##citmin and cit9 still fluoresce to different degrees when fused 
##to the ZEM transcription factor

datadir <- "data/ciBERseq/ZEM-citrine_fusion_control/"
figdir <- "figures/"

data <- read.csv( file.path( datadir, "ZEM-cit_controls_normgated_data.csv" ), header = T)


medians <- aggregate( Citrine.cor ~ cit, data, median ) # median of the citrine measurement of all points per sample
se <- aggregate( Citrine.cor ~  cit, data, function(x) sd(x)/sqrt(length(x)) ) # standard error

medians <- merge( medians, se, by = "cit")
names(medians) <- c("cit", "median", "se")

row.names(medians) <- medians$cit

medians <- medians[c("c9", "ci"),]

cit.labels <- c("slow\ncitrine", "fast\ncitrine")

cairo_pdf( file.path( figdir, "ZEM-cit_fusioncontrols.pdf" ), width = 1.3, height = 1.3, pointsize = 6.5 )
#set margin stuff
par( mex = 0.65 )
par( mar = c(6,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )

mp = barplot(medians$median,
        #beside = TRUE,
        col = "darkgrey",
        ylab = "citrine fluorescence",
        border = NA, axes = F,
        ylim = c(0, 450),
        names.arg = cit.labels,
        mgp = c(3, 2, 0) 
        )
axis(2, lwd = 0.75)

arrows( mp, with( medians, median - se ), mp, with( medians, median + se ), angle=90, code=3, length = 0.02)

dev.off()
