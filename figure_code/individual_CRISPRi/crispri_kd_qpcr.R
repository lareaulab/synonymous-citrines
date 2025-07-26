figdir <- "figures"
datadir <- "data/individual_CRISPRi"

cq <- read.delim( file.path( datadir, "qpcr_20250514.csv"), sep=",")

# A1,2,3 are biological replicates - A1,B1,C1 are technical replicates

cq$Biorep <- rep(1:3, 21)
cq$Abundance <- 2^-cq$Cq

name <- c(Fun12 = "FUN12", Rpg1 = "RPG1", HO = "HO")

# averages of three technical replicates for each target
techmeans <-  aggregate( Abundance ~ Biorep + Target + Sample, cq, mean)

actin_target <- techmeans[ techmeans$Target == "Actin", ]
rpg1_target <- techmeans[ techmeans$Target == "Rpg1", ]
fun12_target <- techmeans[ techmeans$Target == "Fun12", ]

fun12 <- merge( fun12_target, actin_target, by = c("Sample", "Biorep"))
rpg1 <- merge( rpg1_target, actin_target, by = c("Sample", "Biorep"))

fun12 <- fun12[order(fun12$Sample, decreasing=T),]
rpg1 <- rpg1[order(rpg1$Sample, decreasing=F),]

# ratio of target to actin
fun12$Relative <- fun12$Abundance.x / fun12$Abundance.y
rpg1$Relative <- rpg1$Abundance.x / rpg1$Abundance.y

# average of three biological replicates of ratio of target to actin
fun12.avg <- aggregate( Relative ~ Sample, fun12, mean )
rpg1.avg <- aggregate( Relative ~ Sample, rpg1, mean )

fun12.avg <- fun12.avg[order(fun12.avg$Sample, decreasing=T),]
rpg1.avg <- rpg1.avg[order(rpg1.avg$Sample, decreasing=F),]

fun12.p <- t.test( fun12$Relative[fun12$Sample == "HO"], fun12$Relative[fun12$Sample == "Fun12"], var.equal = F )$p.value
fun12.p <- round( fun12.p, 3 )
# 0.012
rpg1.p <- t.test( rpg1$Relative[rpg1$Sample == "HO"], rpg1$Relative[rpg1$Sample == "Rpg1"], var.equal = F )$p.value
rpg1.p <- round( rpg1.p, 3 )
# 0.166

cairo_pdf( file.path( figdir, "fun12_qpcr.pdf"), width = 1.3, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,2,2) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )
# plot the averages as a barplot
max <- max(fun12$Relative)
bp <- barplot( fun12.avg$Relative, names.arg = name[fun12.avg$Sample], 
               xlab = "CRISPRi target",
               ylab = "relative mRNA\nabundance, FUN12",
               ylim = c(0, max)
)
mps <- bp[,1]
names(mps) <- fun12.avg$Sample
# add the individual ratios
points( mps[ fun12$Sample ], fun12$Relative, pch = 19 )
# p value
segments(mps[1], max*1.1, mps[2])
text( mean(mps), max*1.1, labels = paste0("p = ", fun12.p), pos = 3)
dev.off()



cairo_pdf( file.path( figdir, "rpg1_qpcr.pdf"), width = 1.3, height = 1.3, pointsize = 6.5 )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,2,2) )
par( oma = c(0,0.5,2,0) )
par( xpd = NA )
# plot the averages as a barplot
max <- max(rpg1$Relative)
bp <- barplot( rpg1.avg$Relative, names.arg = name[rpg1.avg$Sample], 
               xlab = "CRISPRi target",
               ylab = "relative mRNA\nabundance, RPG1",
               ylim = c(0, max)
)
mps <- bp[,1]
names(mps) <- rpg1.avg$Sample
# add the individual ratios
points( mps[ rpg1$Sample ], rpg1$Relative, pch = 19 )
# p value
segments(mps[1], max*1.1, mps[2])
text( mean(mps), max*1.1, labels = paste0("p = ", rpg1.p), pos = 3)
dev.off()

