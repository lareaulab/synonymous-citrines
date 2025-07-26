#EFL 7/23/25
#written with GPT
#Quantitating citrine western blot

library(plyr)
library(tidyverse)

# Set path and read file
datadir <- "data/WT/"
figdir <- "figures/"

df <- read_csv( file.path( datadir, "western_3-19-25-hi_res-lane_quant.csv" ))

# Pivot: get 700 and 800 channel signals as columns
df_wide <- df %>%
  select(Name, Channel, Signal) %>%
  pivot_wider(names_from = Channel, values_from = Signal) %>%
  mutate(
    Normalized_Citrine = `800` / `700`,
    Name = factor(Name, levels = c("citMin", setdiff(unique(Name), "citMin")))
  )

#make plot
cairo_pdf( file.path( figdir, "western_blot.pdf"), width = 1.75, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(8,6.5,3,3) )
par( mar = c(8,6.5,3,1) )
par( oma = c(0,0.5,1,0) )
par( xpd = F )
mp <- barplot( df_wide$Normalized_Citrine,
               ylim = c(0,1.6),
               space = 0.2,
               axes = F,
               ylab = "citrine / GAPDH",
               xlab = NA,
               border = NA,
               col = "darkgrey"
)
axis(1, labels = df_wide$Name, at = mp, font = 1, tick = F, las = 2, lwd = 0.75)
axis(2, at = seq(0, 1.8, by = 0.4), lwd = 0.75)
title( xlab = "citrine strain", line = 6 )
dev.off()

