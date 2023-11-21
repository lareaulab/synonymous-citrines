##re-analysis of citrine mRNA qPCR data from Tunney, McGlincey, et al 2018
## uses delta delta Cq method (efficiency = 2) 

setwd("../data/WT/WT_mRNA/")

citrine_construct_scores_fname <- "../../codon_scores/citrine_scores_full_model.tsv"
cit <- read.delim(citrine_construct_scores_fname, header=F, row.names = 1)
names(cit) = c("time")
namemap <- c( Y00 = "cit0", Y333 = "cit3", Y66 = "cit6", Y99 = "cit9", MAX = "citmax", MIN = "citmin")


qpcr1 <- read.csv("codopt_Ingolia Lab_2017-09-13 10-16-37_CT009900_sample-info.csv", header = T)
qpcr1.cq <- read.csv("codopt_Ingolia Lab_2017-09-13 10-16-37_CT009900 -  Quantification Summary_0.csv", header = T)

qpcr2 <- read.csv("codopt_MIN069_Ingolia Lab_2017-10-24 14-15-53_CT009900_sample-info.csv", header = T)
qpcr2.cq <- read.csv("codopt_MIN069_Ingolia Lab_2017-10-24 14-15-53_CT009900 -  Quantification Summary_0.csv", header = T)

qpcr3 <- read.csv("codopt_MIN6_redux_Ingolia Lab_2017-11-01 11-34-53_CT009900_sample-info.csv", header = T)
qpcr3.cq <- read.csv("codopt_MIN6_redux_Ingolia Lab_2017-11-01 11-34-53_CT009900 -  Quantification Summary_0.csv", header = T)

# strains min, 333, and max 
qpcr1 <- qpcr1[ with( qpcr1, order( row, column )),]
qpcr1$Cq <- qpcr1.cq$Cq 
qpcr1 <- qpcr1[ with( qpcr1, which(!(Strain == "MAX" & Isolate == 3)) ), ] # remove isolate 3 of max from qpcr1 - pellet was aspirated
qpcr1 <- qpcr1[ with( qpcr1, Strain %in% c("MIN", "Y333", "MAX")), ]

# strains min, 000, 999 (the 666 isolates didn't grow, so we will discard them from this set and re-measure in the next dataset)
qpcr2 <- qpcr2[ with( qpcr2, order( row, column )),]
qpcr2$Cq <- qpcr2.cq$Cq 
qpcr2 <- qpcr2[ with( qpcr2, which(!(strain == "Y99" & isolate == 3)) ), ] # remove isolate 3 of 999 from qpcr2 - discarded from project because of copy number
qpcr2 <- qpcr2[ with( qpcr2, strain %in% c("MIN", "Y00", "Y99")), ]

# strains min and 666
qpcr3.cq <- qpcr3.cq[ with(qpcr3.cq, !is.nan(Cq)), ]
qpcr3 <- qpcr3[ with(qpcr3, order( row, column )), ]
qpcr3$Cq <- qpcr3.cq$Cq 
qpcr3 <- qpcr3[ with( qpcr3, strain %in% c("MIN", "Y66")), ]

# average the Cqs of technical replicates
qpcr1.avgs <- aggregate(Cq ~ Isolate + Strain + amplicon2, qpcr1, mean)
qpcr2.avgs <- aggregate(Cq ~ isolate + strain + amplicon2, qpcr2, mean)
qpcr3.avgs <- aggregate(Cq ~ isolate + strain + amplicon2, qpcr3, mean)
names(qpcr1.avgs) <- names(qpcr2.avgs)

# citrine to mcherry delta Cq; then delta delta Cq against average citmin from each plate
qpcr3.ratios <- qpcr3.avgs[1:6, c("strain","isolate")]
qpcr3.ratios$dcq <- qpcr3.avgs$Cq[1:6] - qpcr3.avgs$Cq[7:12] # citrine - mCherry delta Cq
min3 <- mean( qpcr3.ratios$dcq[ with( qpcr3.ratios, strain == "MIN" ) ] )
qpcr3.ratios$ddcq <- qpcr3.ratios$dcq - min3 # normalize to average delta Cq for citMin for each plate so we can compare between plates
qpcr3.ratios$strain <- namemap[qpcr3.ratios$strain] # standardize the strain naming 

qpcr2.ratios <- qpcr2.avgs[1:8, c("strain","isolate")]
qpcr2.ratios$dcq <- qpcr2.avgs$Cq[1:8] - qpcr2.avgs$Cq[9:16]
min2 <- mean( qpcr2.ratios$dcq[ with( qpcr2.ratios, strain == "MIN" ) ] )
qpcr2.ratios$ddcq <- qpcr2.ratios$dcq - min2
qpcr2.ratios$strain <- namemap[qpcr2.ratios$strain]

qpcr1.ratios <- qpcr1.avgs[1:8, c("strain","isolate")]
qpcr1.ratios$dcq <- qpcr1.avgs$Cq[1:8] - qpcr1.avgs$Cq[9:16]
min1 <- mean( qpcr1.ratios$dcq[ with( qpcr1.ratios, strain == "MIN" ) ] )
qpcr1.ratios$ddcq <- qpcr1.ratios$dcq - min1
qpcr1.ratios$strain <- namemap[qpcr1.ratios$strain]

# combine three plates and compute ratio of citrine to mCherry; keep the MIN data from the first plate
ratios <- rbind( qpcr1.ratios, 
                 qpcr2.ratios[ qpcr2.ratios$strain != "citmin", ],
                 qpcr3.ratios[ qpcr3.ratios$strain != "citmin", ] )
ratios <- within( ratios, ratio <- 2^-ddcq )
#ratios <- within( ratios, citscore <- cit.scores[strain] )
ratios <- within( ratios, citscore <- cit[strain,] )

write.csv(ratios, "WT_mRNA_ratios.csv", row.names=FALSE)
