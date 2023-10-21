##Calculating and plotting the translation efficiency (protein/mRNA) for 
##each SL-citrine construct.

library(dplyr)
library(plyr)
library(tidyr)

setwd(".../data/stemloop/")

#qPCR (mRNA) data:
mRNA = read.csv("SL_mRNA/J.WT.mRNA.csv", header = T)
colnames(mRNA)[4] = "Normalized.mRNA"

#fluorescence (protein) data:
gated = read.csv("SL_protein/normgated_data_bg_corrected.csv", header=T, row.names=1)
  medians <- aggregate( ratio ~ clone + cit + strain, gated, median)
  medians <- medians[ medians$strain != "BY4743", ]
  #add elongation rate for each citrine
  rates <- c( citmin = 164.6625671,
              cit0 = 234.2350459, 
              cit3 = 262.7644388, 
              cit6 = 269.5761919, 
              cit9 = 302.6383892, 
              citmax = 388.8861452 )
  medians$elongation_time = rates[medians$cit] 
  j.wt.medians <- medians[ medians$strain !="HC1g", ]
  j.wt.medians <- subset(j.wt.medians, select = -c(elongation_time))
  j.wt.medians$strain[j.wt.medians$strain == "HC1j"] <- "HC1J"
    fluorescence <- j.wt.medians
      colnames(fluorescence)[1] = "Clone"
      colnames(fluorescence)[2] = "Cit"
      colnames(fluorescence)[3] = "Strain"
      colnames(fluorescence)[4] = "Normalized.fluorescence"
        fluorescence$Cit[fluorescence$Cit == "citmax"] <- "citMax"
        fluorescence$Cit[fluorescence$Cit == "citmin"] <- "citMin"

#making the translation efficiency dataframe:
translation <- merge(fluorescence, mRNA, by = c("Strain", "Cit", "Clone"))
rates <- c( citMin = 164.6625671,
            cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919,
            cit9 = 302.6383892, 
            citMax = 388.8861452)
translation$elongation_time <- rates[translation$Cit] 

translation$TE <- translation$Normalized.fluorescence/translation$Normalized.mRNA
#averages:
wt.avg = aggregate( TE ~ Cit + Strain + elongation_time, translation[translation$Strain == "WT",], mean)
j.avg = aggregate( TE ~ Cit + Strain + elongation_time, translation[translation$Strain == "HC1J",], mean)


#basic (main figure) plot:
pdf("SLvsWT.TEaverages.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( translation$elongation_time, translation$TE, 
      col = ifelse(translation$Strain == "WT", "#3584c6", "#e95c64"),
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0,0.8),
      xlab = "",
      ylab = "translation efficiency\n(fluorescence/mRNA)"
)
axis( 1 )
axis( 2, seq(0,0.8,0.133), labels = c(0.0, NA, 0.3, NA, 0.5, NA, 0.8) )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

lines( wt.avg$elongation_time, wt.avg$TE, lwd = 1.5, col = "#3584c6")
lines( j.avg$elongation_time, j.avg$TE, lwd = 1.5, col = "#e95c64")

legend( "topright", pch = 20, 
        legend = c("no hairpin", "strong hairpin"), 
        col=c("#3584c6", "#e95c64"),
        cex=0.55,
        bty = "n")
dev.off()


#normalizing data for each SL vs WT citrine to its respective citmax average:
means <- aggregate(TE ~ Strain + Cit, translation, mean)
citMax.avg <- means[means$Cit == "citMax",]
JcitMax.avg <- citMax.avg[citMax.avg$Strain == "HC1J",]
JcitMax.avg <- as.numeric(JcitMax.avg[,3])
WTcitMax.avg <- citMax.avg[citMax.avg$Strain == "WT",]
WTcitMax.avg <- as.numeric(WTcitMax.avg[1,3])

normedTranslation <- translation[, c("Strain", "Cit", "Clone", "elongation_time")]

normedTranslation$normedTE <- ifelse(translation$Strain <= "HC1J",  translation$TE/JcitMax.avg, 
                                     ifelse(translation$Strain <= "WT", translation$TE/WTcitMax.avg, NA)
)
#averages:
normed.wt.avg = aggregate( normedTE ~ Cit + Strain + elongation_time, normedTranslation[normedTranslation$Strain == "WT",], mean)
normed.j.avg = aggregate( normedTE ~ Cit + Strain + elongation_time, normedTranslation[normedTranslation$Strain == "HC1J",], mean)

#plotting normalized TE:
pdf("SLvsWT.TEnormedtocitMax.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( normedTranslation$elongation_time, normedTranslation$normedTE, 
      col = ifelse(translation$Strain == "WT", "#3584c6", "#e95c64"),
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0,10),
      xlab = "",
      ylab = "translation efficiency\n(fluorescence/mRNA)"
)
axis( 1 )
#axis( 2 )
axis( 2, seq(0,10,1.666666666), labels = c(0.0, NA, NA, 5, NA, NA, 10))
title( main = "TE normed to citMax", xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

#lines( normed.wt.avg$elongation_time, normed.wt.avg$normedTE, lwd = 1.5, col = "#3584c6")
#lines( normed.j.avg$elongation_time, normed.j.avg$normedTE, lwd = 1.5, col = "#e95c64")

legend( "topright", pch = 20, 
        legend = c("no hairpin", "strong hairpin"), 
        col=c("#3584c6", "#e95c64"),
        cex=0.55,
        bty = "n")

dev.off()


#calculating SL/WT translation efficiency
dividedTE <- means %>%
  pivot_wider(
    id_cols = c("Cit"),
    names_from = Strain,
    values_from = c(TE),
    unused_fn = NULL
  )

dividedTE$ratio <- dividedTE$HC1J/dividedTE$WT

rates <- c( cit0 = 234.2350459, 
            cit3 = 262.7644388, 
            cit6 = 269.5761919,
            cit9 = 302.6383892, 
            citMax = 388.8861452,
            citMin = 164.6625671)

dividedTE$elongation_rate <- rates[dividedTE$Cit] 

HP.WTratiosTE <- dividedTE$ratio
trend <- lm(HP.WTratiosTE~rates)

#plotting change in TE across citrines:
pdf("SLvsWT_TE_ratios.pdf", width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot( dividedTE$elongation_rate, dividedTE$ratio, 
      col = "#aaabab",
      pch = 20,
      #cex = 0.6,
      axes = F,
      xlim = c(150,400),
      ylim = c(0, 1),
      xlab = "",
      ylab = "SL / WT\n TE ratio"
)
#abline(trend,col='#b2b2b2', lty = 2)

axis( 1 )
axis( 2 )
#axis( 2, seq(0,2,0.33), labels = c(0.0, NA, 0.6, NA, 1.3, NA, 2))
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )

dev.off()

```