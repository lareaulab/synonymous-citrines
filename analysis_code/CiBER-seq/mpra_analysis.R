options(stringsAsFactors = FALSE)
library("mpra")

datadir <- "data/ciBERseq/ciBERseq/"

# data
counts <- read.delim( file.path( datadir, "counts_all.txt"), row.names = 1, header = T )

# other information
grna <- read.delim( file.path( datadir, "grna-assign-barcode-grna-good.txt" ), row.names = 1, header = T)
target <- read.delim( file.path( datadir, "guide-good-targets-171010.txt" ), header = T)
empty <- read.delim( file.path( datadir, "mod-probably-empty.txt" ), header = FALSE)

# if (!file.exists("SGD_features.tab")) {
#   sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile = "../SGD_features.tab")
# }
sgd <- read.delim( file.path( datadir, "SGD_features.tab" ), header = FALSE, quote = "",
                  col.names = c("sgdid", "type", "qual", "name", "gene", "alias",
                                "parent", "sgdid2", "chrom", "start", "end",
                                "strand", "genpos", "cver", "sver", "desc"))

cutoff <- 24

# fast and slow
# barcodes with enough reads in both DNA samples only
counts <- counts[ counts$fast_DNA_pre_L > cutoff & counts$fast_DNA_pre_R > cutoff & counts$slow_DNA_pre_L > cutoff & counts$slow_DNA_pre_R > cutoff, ]
# match barcodes from input counts to guide names (genes or empty controls)
counts <- counts[row.names(counts) %in% c(row.names(grna), empty$V1) ,]

# mpra package data for link between barcode and guide
guide_eid <- grna[match(row.names(counts), row.names(grna)),]
guide_eid <- ifelse(is.na(guide_eid), 
                    paste0("Empty", seq(1, length(guide_eid))),
                    guide_eid)

run_mpra <- function( data, sample ) {
  dna <- data[ , grep( paste0(sample, "_DNA"), colnames(data) )]
  rna <- data[ , grep( paste0(sample, "_RNA"), colnames(data) )]
  
  colnames(dna) <- sub( paste0(sample, "_DNA_"), "", colnames(dna))
  colnames(rna) <- sub( paste0(sample, "_RNA_"), "", colnames(rna))
  
  mpraset <- MPRASet(DNA = dna, RNA = rna, 
                         eid = guide_eid, eseq = NULL, 
                         barcode = row.names(counts))
  design <- data.frame(intcpt = 1, 
                           induced = grepl("post", colnames(dna)),
                           tstat = grepl("L", colnames(dna)))
  mpralm_fit <- mpralm(object = mpraset, design = design,
                           aggregate = "sum",
                           normalize = TRUE, model_type = "indep_groups", plot = TRUE)

  out <- topTable(mpralm_fit, coef = 2, number = Inf)
  
  out$yorf <- target[match(row.names(out), target$Guide), "Yorf1"]
  out$yorfs <- target[match(row.names(out), target$Guide), "Yorfs"]
  out$gene <- sgd[match(out$yorf, sgd$name), "gene"]
  out$desc <- strtrim(sgd[match(out$yorf, sgd$name), "desc"], 160)
  
  out <- out[!is.na(out$yorf),]
  return(out)
}

fast <- run_mpra( counts, "fast")
slow <- run_mpra( counts, "slow")


# put the two in the same order (tested that these lists are identical)
fast = fast[order(row.names(fast)),]
slow = slow[order(row.names(slow)),]

# if there's no gene name, give it the yorf name as a gene name
fast$gene[which(fast$gene == "")] = fast$yorf[which(fast$gene == "")]
slow$gene[which(slow$gene == "")] = slow$yorf[which(slow$gene == "")]

#Write output
write.table( x = fast[,!grepl("desc", colnames(fast))], 
             file = file.path( datadir, "fast_mpra_results_byguide.txt" ),
             row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table( x = slow[,!grepl("desc", colnames(slow))], 
             file = file.path( datadir, "slow_mpra_results_byguide.txt" ),
             row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

up_cutoff <- log2(2)
diff_cutoff <- log2(1.5)

slowup = slow$logFC > up_cutoff & (slow$logFC - fast$logFC) > diff_cutoff
slowup_strict = slow$logFC > log2(3) & (slow$logFC - fast$logFC) > log2(1.5)

# background set for GO analysis
write.table( unique( slow$yorf ), 
             file = file.path( datadir, "all_yorfs.txt" ), 
             quote = F, row.names = F, col.names = F)
# genes with direction of effect we want
write.table( unique( slow$yorf[slowup] ), 
             file = file.path( datadir, "slow_up.txt" ), 
             quote = F, row.names = F, col.names = F)
write.table( unique( slow$yorf[slowup_strict] ),
             file = file.path( datadir, "slow_up_strict.txt" ),
             quote = F, row.names = F, col.names = F)

# most-up set: 
# "FUN12" "OAF1"  "KIP1"  "RPG1"  "TRP1"  "SEC3"  "SEC27" "MPT5"  "TPN1"  "NOP7"  "NET1" 
# "SMF3"  "NIP1"  "TOP2" 
# 	formation of cytoplasmic translation initiation complex (3)	3.84e-3
#   translational initiation (4)	9.16e-2
# fun12, rpg1, nip1, [mpt5]
