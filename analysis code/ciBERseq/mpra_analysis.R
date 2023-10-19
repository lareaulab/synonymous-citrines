options(stringsAsFactors = FALSE)
library("mpra")

# data
counts <- read.delim( "counts_all.txt", stringsAsFactors = FALSE, row.names = 1, header = T )

# other information
grna <- read.delim("grna-assign-barcode-grna-good.txt", stringsAsFactors = FALSE, row.names = 1, header = T)
target <- read.delim("guide-good-targets-171010.txt", stringsAsFactors = FALSE, header = T)
empty <- read.delim("mod-probably-empty.txt", stringsAsFactors = FALSE, header = FALSE)

if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile = "../SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header = FALSE, quote = "",
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
  out$Unind <- out$AveExpr - out$logFC/2 # RNA to DNA ratio in the uninduced sample, for quality control
  
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
write.table(x = fast[,!grepl("desc", colnames(fast))], file = "fast_mpra_results_byguide.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = slow[,!grepl("desc", colnames(slow))], file = "slow_mpra_results_byguide.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


# all genes within a strict set of bounds: slow foldchange > 1 or 1.5, slow foldchange - fast foldchange > 0.5
slowupfc1 = slow$logFC > 1 & (slow$logFC - fast$logFC) > 0.5
slowupfc1.5 = slow$logFC > 1.5 & (slow$logFC - fast$logFC) > 0.5

write.table(slow$gene[slowupfc1], file = "slow_all_targets_up_fc1.txt", quote = F, row.names = F, col.names = F)
write.table(slow$gene[slowupfc1.5], file = "slow_all_targets_up_fc1.5.txt", quote = F, row.names = F, col.names = F)

