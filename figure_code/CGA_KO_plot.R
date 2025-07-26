# Plot the citrine/mCherry fluorescence ratio for CGA KO and related strains
#EFL 7-23-25 written with GPT
library(plyr)
library(tidyverse)

datadir <- "data/CGA_reporter_control/"
figdir <- "figures/"

gated_data <- read_csv( file.path( datadir, "gated_data_bg.csv" ))

# Calculate median ratios
gated_data$ratio.A <- with(gated_data, Citrine.cor / mCherry.cor)

medians <- gated_data %>%
  filter(strain != "BY4743") %>%
  group_by(sample, clone, cit, strain) %>%
  summarise(ratio = median(ratio.A),
            Citrine = median(Citrine.cor),
            mCherry = median(mCherry.cor), .groups = "drop")

medians <- medians[ medians$strain != "BY4743", ]

# switch name order for double knockout
medians$strain[medians$strain == "syh1smy2"] <- "smy2syh1"

# Set plotting order
strain_order <- c("hel2", "hel2syh1", "smy2syh1", "syh1", "WT")
medians$strain <- factor(medians$strain, levels = strain_order)
strain_numeric <- as.numeric(medians$strain)

# Jitter for plotting, but base axis labels on center positions
set.seed(1)  # for reproducibility
xpoints <- strain_numeric + runif(length(strain_numeric), -0.15, 0.15)

# Colors (optional)
cols <- rep("darkgrey", nrow(medians))

# Axis labels
labels <- c(expression(paste(italic("hel2"), Delta)), 
            expression(paste(italic("hel2"), Delta, italic("syh1"), Delta)),
            expression(paste(italic("smy2"), Delta, italic("syh1"), Delta)),
            expression(paste(italic("syh1"), Delta)),
            "WT")

# Plot
cairo_pdf(file.path(figdir, "CGA_KO.pdf"), width = 1.75, height = 1.4, pointsize = 6.5)
par(mex = 0.65)
par(mar = c(9, 6.5, 1, 1))
par(oma = c(0, 0.5, 1, 0))
par(xpd = NA)

plot(xpoints, medians$ratio,
     ylim = c(0, 0.6),
     pch = 20,
     #col = cols,
     xlab = NA,
     ylab = "GFP / mCherry\nfluorescence ratio",
     axes = FALSE
)

# y-axis
axis(2, at = seq(0, 0.6, by = 0.2), lwd = 0.75)

# x-axis: perfectly aligned with mean position of each strain group
axis(1, at = 1:length(labels), labels = labels,
     col = NA, lwd = 0.75, las = 2)

dev.off()

