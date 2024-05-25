library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

setwd(getSrcDirectory(function(){})[1])

data <- read.csv('/data/TCGA-LUAD.mutect2_sbs.tsv', sep='\t')
head(data)
sigs.input <- mut.to.sigs.input(mut.ref = data, 
                                sample.id = "Sample_ID", 
                                chr = "chrom", 
                                pos = "start", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

samples = rownames(as.data.frame(sigs.input))
signatures_all = sapply(samples, function(x) whichSignatures(tumor.ref = as.data.frame(sigs.input),
                                                             signatures.ref=signatures.cosmic,
                                                             sample.id = x,
                                                             contexts.needed = TRUE,
))

colnames <- colnames(as.data.frame(signatures_all[,"TCGA-97-7938-01A"]$weights))

pivot_table <- data.frame(matrix(ncol = length(colnames), nrow = 0))
colnames(pivot_table) <- colnames

for (i in colnames(signatures_all)) {
  pivot_table <- rbind(pivot_table, as.data.frame(signatures_all[,i]$weights))
}

signature_groups <- list(
  "Spontaneous_deamination_of_5-methylcytosine" = c("Signature.1"),
  "AID/APOBEC_family_of_cytidine_deaminases" = c("Signature.2", "Signature.13"),
  "Failure_of_DNA_double-strand_break-repair_by_homologous_recombination" = c("Signature.3"),
  "Smoking" = c("Signature.4"),
  "Unknown" = c("Signature.5", "Signature.8", "Signature.12", "Signature.17", "Signature.18", "Signature.19", "Signature.21", "Signature.23", "Signature.24", "Signature.25", "Signature.27", "Signature.28", "Signature.30"),
  "Defective_DNA_mismatch_repair" = c("Signature.6", "Signature.15", "Signature.20", "Signature.26"),
  "Ultraviolet_light_exposure" = c("Signature.7"),
  "Polymerase_activity" = c("Signature.9", "Signature.10"),
  "Alkylating_agents" = c("Signature.11"),
  "Aristolochic_acid_exposure" = c("Signature.22"),
  "Tobacco_chewing" = c("Signature.29")
)

combine_signatures <- function(signatures_table, signature_groups) {
  combined_table <- data.frame(matrix(ncol = length(signature_groups), nrow = nrow(signatures_table)))
  colnames(combined_table) <- names(signature_groups)
  
  for (group_name in names(signature_groups)) {
    combined_table[, group_name] <- rowSums(signatures_table[, signature_groups[[group_name]], drop = FALSE])
  }
  rownames(combined_table) <- rownames(signatures_table)
  
  return(combined_table)
}

combined_table <- combine_signatures(pivot_table, signature_groups)
write.table(combined_table, '/data/combined_table_signatures.tsv', sep='\t')

data1 <- read.csv('/data/TCGA-LUAD.GDC_phenotype.tsv', sep='\t', row.names = 1)
data2 <- data.frame("tobacco_smoking_history" = data1$tobacco_smoking_history)
row.names(data2) <- row.names(data1)

combined_table$ID <- rownames(combined_table)
combined_table <- combined_table[, c("ID", names(combined_table)[-ncol(combined_table)])]

data2$ID <- rownames(data2)
data2 <- data2[, c("ID", names(data2)[-ncol(data2)])]

combined_table_tobaco <- merge(combined_table, data2, by = "ID")

combined_table_tobaco_without_na <- combined_table_tobaco %>% drop_na()
write.table(combined_table_tobaco_without_na, '/data/combined_table_signatures_tobacco.tsv', sep='\t')