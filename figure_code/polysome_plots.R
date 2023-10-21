##plotting the normalized citrine read counts from RNAseq data of each polysome fraction over the potential traces from the fractionation

setwd( "../../data/polysomes/" )

traceB <- read.delim( "fractionation_traces/WTi9-B-i.csv", header = T, sep = ",", col.names = c("time", "potential") )
traceC <- read.delim( "fractionation_traces/WTi9-C-i.csv", header = T, sep = ",", col.names = c("time", "potential") )

cit9 <- read.delim( "RNAseq/cit9_counts.txt", header=F, col.names = c("fraction", "slow"), row.names = 1 )
citmin <- read.delim( "RNAseq/citmin_counts.txt", header=F, col.names = c("fraction", "fast"), row.names = 1 )

human <- read.delim( "RNAseq/human_perfect_totals.txt", header=F, col.names = c("fraction", "human"), row.names = 1 )

norm <- cbind( cit9/human, citmin/human )

repB <- norm[2:6,]
repC <- norm[8:12,]

repB$slow <- repB$slow / sum(repB$slow)
repB$fast <- repB$fast / sum(repB$fast)

repC$slow <- repC$slow / sum(repC$slow)
repC$fast <- repC$fast / sum(repC$fast)

fracB <- c(2.2, 3.15, 3.75, 4.166667, 4.516667)
fracB_val <- c(0.606, 0.284, 0.248, 0.201, 0.171)

fracB_shift <- fracB - 0.1
fracB_shift_x <- sapply( round(fracB_shift, 3), function(x) which( round(traceB$time, 3) == x ))
fracB_val_shift <- traceB$potential[ fracB_shift_x ]

divB <- c(fracB_shift, traceB$time[nrow(traceB)])
midpointsB <- sapply( 1:5, function(i)( (divB[i] + divB[i+1])/2 ) )

fracC = c(2.866667, 3.733333, 4.333333, 4.8, 5.116667)
fracC_val <- c(0.784, 0.325, 0.287, 0.267, 0.151)

fracC_shift <- fracC - 0.1
fracC_shift_x <- sapply( round(fracC_shift, 3), function(x) which( round(traceC$time, 3) == x ))
fracC_val_shift <- traceC$potential[ fracC_shift_x ]

divC <- c(fracC_shift, traceC$time[nrow(traceC)])
midpointsC <- sapply( 1:5, function(i)( (divC[i] + divC[i+1])/2 ) )

yticks = seq(0, 0.6, by = 0.2)
cols = c( slow = "darkorange2", fast = "magenta3")

pdf("../../figures/polysome_rep_B.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,1,2,6) )
par( oma = c(0,0.5,1,0.5) )
par( xpd = F )
plot( traceB$time, traceB$potential, 
      type = "l",
      ylim = c(0, 1),
      xlim = c(0.5,5.5),
      axes = F,
      lwd = 2,
      col = "lightgrey",
      xlab = "polysome fraction",
      ylab = NA )
segments(fracB_shift, 0, fracB_shift, fracB_val_shift, lwd = 2, col = "lightgrey")
axis( 1, at = midpointsB, labels = c(1,2,3,4,"5+"), tick = F ) 
lines( midpointsB, 1.5 * repB$slow, lwd = 2, col = cols["slow"]) #"darkorange2" )
lines( midpointsB, 1.5 * repB$fast, lwd = 2, col = cols["fast"])# "magenta3" )
axis( 4, at = 1.5 * yticks, labels = yticks, lwd = 0, lwd.ticks = 1, las = 1)
text( midpointsB[1], 1.5*sum(repB[ 1,])/2 + c(0.06, -0.06), labels = c("slow citrine", "fast citrine"), col = cols[ c("slow", "fast") ], pos = 2)
mtext("mRNA as fraction of total", side = 4, line = 3.5)
mtext("replicate 1", 3, adj = 0)
dev.off()

pdf("../../figures/polysome_rep_C.pdf", width = 2.5, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(7,1,2,6) )
par( oma = c(0,0.5,1,0.5) )
par( xpd = F )
plot( traceC$time, traceC$potential, 
      type = "l",
      ylim = c(0, 1),
      xlim = c(1,6),
      axes = F,
      lwd = 2,
      col = "lightgrey",
      xlab = "polysome fraction",
      ylab = NA )
segments(fracC_shift, 0, fracC_shift, fracC_val_shift, lwd = 2, col = "lightgrey")
axis( 1, at = midpointsC, labels = c(1,2,3,4,"5+"), tick = F ) 
lines( midpointsC, 1.5 * repC$slow, lwd = 2, col = "darkorange2" )
lines( midpointsC, 1.5 * repC$fast, lwd = 2, col = "magenta3" )
axis( 4, at = 1.5 * yticks, labels = yticks, lwd = 0, lwd.ticks = 1, las = 1)
text( midpointsC[1], 1.5*sum(repC[ 1,])/2 + c(0.06, -0.06), labels = c("slow citrine", "fast citrine"), col = cols[ c("slow", "fast") ], pos = 2)
mtext("mRNA as fraction of total", side = 4, line = 3.5)
mtext("replicate 2", 3, adj = 0)
dev.off()
