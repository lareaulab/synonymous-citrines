##plotting the normalized citrine read counts from RNAseq data of each polysome fraction over the potential traces from the fractionation

datadir <- "data/polysomes/"
figdir <- "figures/"

traceB <- read.delim( file.path( datadir, "fractionation_traces/i-x_B2.csv" ), header = T, sep = ",", col.names = c("time", "potential") )
traceC <- read.delim( file.path( datadir, "fractionation_traces/i-x_C2.csv" ), header = T, sep = ",", col.names = c("time", "potential") )

citB <- read.table( file.path( datadir, "RNAseq/citrine_counts_B.txt" ), header=F, col.names = c("count", "frac", "cit") )
citC <- read.table( file.path( datadir, "RNAseq/citrine_counts_C.txt" ), header=F, col.names = c("count", "frac", "cit") )

humanB <- read.table( file.path( datadir, "RNAseq/human_counts_B.txt" ), header=F, col.names = c("count", "frac") )
humanC <- read.table( file.path( datadir, "RNAseq/human_counts_C.txt" ), header=F, col.names = c("count", "frac") )



# hack to make it possible to look up reading by x axis time (avoid floating point issues)
traceB$time <- as.integer(round(1000 * traceB$time, 0))
traceC$time <- as.integer(round(1000 * traceC$time, 0))


# fragile hack to get it into the right shape
citB_counts <- as.data.frame( matrix( citB$count, nrow = 6))
colnames(citB_counts) <- unique(citB$frac)
row.names(citB_counts) <- unique(citB$cit)

citC_counts <- as.data.frame( matrix( citC$count, nrow = 6))
colnames(citC_counts) <- unique(citC$frac)
row.names(citC_counts) <- unique(citC$cit)


humanB_counts <- as.data.frame( matrix( rep(humanB$count, each=6), nrow = 6))
colnames(humanB_counts) <- humanB$frac

humanC_counts <- as.data.frame( matrix( rep(humanC$count, each=6), nrow = 6))
colnames(humanC_counts) <- humanC$frac


normB <- citB_counts/humanB_counts
percentB <- normB / rowSums(normB)

normC <- citC_counts/humanC_counts
percentC <- normC / rowSums(normC)


fracB <- c( 2.25, 3.3, 3.966667, 4.4, 4.766667, 5.05, 5.266667, 5.483333 )
fracB <- as.integer(round(1000 * fracB, 0))

# 0.1 minutes is the dead time between the spectrophotometer and the end of the outflow tube
fracB_shift <- fracB - (0.1 * 1000)
# y values at those timepoints
fracB_val_shift <- traceB$potential[ sapply( fracB_shift, function(x) which(traceB$time == x) ) ]

# add endpoints to set up fraction boundaries for the gradient plot
#divB <- c( 0, fracB_shift, max(traceB$time) )
divB <- c( 1100, fracB_shift, 6000 )
midpointsB <- sapply( 1:9, function(i)( (divB[i] + divB[i+1])/2 ) )

fracC = c( 2.45, 3.5, 4.15, 4.583333, 4.933333, 5.233333, 5.433333, 5.6 )
fracC <- as.integer(round(1000 * fracC, 0))

fracC_shift <- fracC - (0.1 * 1000)
# y values at those timepoints
fracC_val_shift <- traceC$potential[ sapply( fracC_shift, function(x) which(traceC$time == x) ) ]

#divC <- c( 0, fracC_shift, max(traceC$time) )
divC <- c( 1100, fracC_shift, 6000 )
midpointsC <- sapply( 1:9, function(i)( (divC[i] + divC[i+1])/2 ) )

yticks = seq(0, 0.2, by = 0.1)

cols <- c( citmin = "magenta3", cit0 = "royalblue2", cit3 = "green3", cit6 = "gold1", cit9 = "darkorange2", citmax = "red2")

# scaling to get the line graphs underneath the polysome trace
mult = 5.5 
offset = -1.5


# fix the gain problem for trace B
# measurements before ~2.5 minutes were taken with gain of 0.05 instead of 0.1
# rescale to put it all on the same axis
changepoint = 2400
c <- which( traceB$time == changepoint)
n <- length(traceB$time)
adj_val <- fracB_val_shift
adj_val[ which(fracB_shift < changepoint) ] = adj_val[ which(fracB_shift < changepoint) ]/2



cairo_pdf( file.path( figdir, "polysomes.pdf" ), width = 4.5, height = 1.3, pointsize = 6.5)
par( mex = 0.65 ) # sets margin stuff
par( oma = c(1,1,1,8) )
par( mar = c(2,1,1,1))
par( xpd = F )
par( mfrow = c(1, 2))

plot( traceC$time, traceC$potential, 
      type = "l",
      ylim = c(offset, 1),
      xlim = c(1100, 6000),
      axes = F,
      lwd = 1.5,
      col = "lightgrey",
      xlab = "polysome fraction",
      ylab = NA,
      mgp = c(1, 1, 0)
)
segments(fracC_shift, 0, fracC_shift, fracC_val_shift, lwd = 1, col = "lightgrey")
axis( 1, at = midpointsC[seq(1,length(midpointsC),by=2)], labels = c(0,1,2,3,4,5,6,7,"8+")[seq(1,length(midpointsC),by=2)], tick = F, line = -1 ) 
axis( 1, at = midpointsC[seq(2,length(midpointsC),by=2)], labels = c(0,1,2,3,4,5,6,7,"8+")[seq(2,length(midpointsC),by=2)], tick = F, line = -1 ) 

lapply( 1:6, function(x){ lines( midpointsC, mult * percentC[x,] + offset, col=cols[ row.names(percentC)[x] ], lwd = 1.5) } )
mtext("replicate 1", 3, adj = 0)

plot( traceB$time[c+1:n], traceB$potential[c+1:n], 
      type = "l",
      ylim = c(offset, 1),
      xlim = c(1000, 6000),
      axes = F,
      lwd = 1.5,
      col = "lightgrey",
      xlab = "polysome fraction",
      ylab = NA,
      mgp = c(1, 1, 0)
)
lines( traceB$time[1:c], traceB$potential[1:c] / 2,
       lwd = 1.5,
       col = "lightgrey")
segments(fracB_shift, 0, fracB_shift, adj_val, lwd = 1, col = "lightgrey")
axis( 1, at = midpointsB[seq(1,length(midpointsB),by=2)], labels = c(0,1,2,3,4,5,6,7,"8+")[seq(1,length(midpointsB),by=2)], tick = F, line = -1 ) 
axis( 1, at = midpointsB[seq(2,length(midpointsB),by=2)], labels = c(0,1,2,3,4,5,6,7,"8+")[seq(2,length(midpointsB),by=2)], tick = F, line = -1 ) 

lapply( 1:6, function(x){ lines( midpointsB, mult * percentB[x,] + offset, col=cols[ row.names(normB)[x] ], lwd = 1.5) } )
mtext("replicate 2", 3, adj = 0)

axis( 4, at = mult * yticks + offset, labels = NA, lwd = 0, lwd.ticks = 0.75, line = 0, tck = 0.025)
axis( 4, at = (mult * yticks + offset)[c(1,3)], labels = c(0, 10, 20)[c(1,3)], lwd = 0, las = 0, line = -0.75)
axis( 4, at = (mult * yticks + offset)[2], labels = c(0, 10, 20)[2], lwd = 0, las = 0, line = -0.75)
mtext("% of total", side = 4, line = 1.5, adj = 0.04)

par(xpd = NA)
legend( "bottomright", col = cols, pch = 20, legend = c("fast", NA, NA, NA, "slow", "slowest"), bty = "n" , inset = c(-0.5, -0.3))
dev.off()
