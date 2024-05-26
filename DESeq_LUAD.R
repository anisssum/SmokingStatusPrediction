# Load required libraries
library(dplyr)
library(DESeq2)
library(tidyverse)
library(ggrepel)

setwd(getSrcDirectory(function(){})[1])

# Read and prepare index data
index <- read.csv("data/transpose_ID_SmokingHistory.csv", sep = " ", header = FALSE)
colnames(index) <- head(index, 1)
index <- index[-1, ]

# Read and prepare sample type data
sample_type <- head(read.csv("data/transpose_ID_SampleType.csv", sep = " ", header = FALSE), 2)
cols <- head(sample_type, 1)
colnames(sample_type) <- cols
sample_type <- sample_type[-1, ]
sample_type[sample_type == 'Primary'] <- 1
sample_type[sample_type == 'Recurrent'] <- 1
sample_type[sample_type == 'Solid'] <- 0
sample_type <- lapply(sample_type, as.numeric)
sample_type <- as.data.frame(sample_type)
colnames(sample_type) <- cols

# Read and prepare gene data
data <- read.csv("data/tcga.gene_sums.LUAD.R109", sep = "\t", header = FALSE)
cols <- head(data, 1)
colnames(data) <- cols
data <- data[-1, ]
data <- lapply(data, as.numeric)
data <- as.data.frame(data)
colnames(data) <- cols

# Read and prepare gene annotation data
gene_refseq <- read.csv("data/gene_annotation_table.txt", sep = "\t", header = FALSE)
colnames(gene_refseq) <- head(gene_refseq, 1)
gene_refseq <- gene_refseq[-1, ]

# Merge gene data with annotation
data_refseq <- merge(data, gene_refseq, by = 'gene_id', all = TRUE)
data_refseq <- data_refseq %>% select(gene_name, everything())
data <- subset(data_refseq, select = -c(gene_id))

# Combine data with sample type information
data <- bind_rows(data, sample_type)
data_t <- t(data)
colnames(data_t) <- data_t[1, ]
data_t <- data_t[-1, ]

# Filter to keep only tumor samples
data_t1 <- data_t[!(data_t[, 54043] == 0), ]
data_1 <- t(data_t1)
data_1 <- as.data.frame(data_1[-nrow(data_1), ])
data_1[] <- lapply(data_1, function(x) {
  if (is.character(x)) as.numeric(as.character(x)) else x
})

# Combine data with index
data_2 <- bind_rows(data_1, index)

# Filter out columns with NA values
data_2 <- data_2[, colSums(is.na(data_2)) == 0]

# Filter data based on training indices
train_index <- read.csv("data/train_id.csv", header = FALSE)
data_2 <- data_2[, colnames(data_2) %in% train_index$V1]

# Sort data
data_t <- t(data_2)
data_t_sort <- data_t[order(data_t[, 54043]), ]
data_sort <- t(data_t_sort)
data_sort_1 <- data_sort[-nrow(data_sort), ]
# write.table(data_t_sort, file = "/data/data_t_sort_all.txt", sep = "\t")

# Prepare metadata for DESeq2
metaF <- data.frame(samples = colnames(data_sort_1),
                    type = c(rep('13', times = 58), rep('245', times = 97), rep('13', times = 129), rep('245', times = 148), rep('245', times = 3)))

# Create DESeq2 dataset
counts_mtx <- data_sort_1 %>% as.matrix()
dds <- DESeqDataSetFromMatrix(countData = counts_mtx, colData = metaF, design = ~ type)

# Normalize data
dds <- estimateSizeFactors(dds)

# Perform variance stabilizing transformation and plot PCA
transformed_dds <- vst(dds)
plotPCA(transformed_dds, intgroup = c('type'))

# Perform DESeq analysis
dds$type <- relevel(dds$type, ref = '13')
dds <- DESeq(dds)
res <- results(dds, name = resultsNames(dds)[2])
write.table(res, file = "/data/DESeq_res.csv")

# Process DESeq results
res2 <- res %>% as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  mutate(FC = 2^log2FoldChange,
         FC_sig = case_when(log2FoldChange > 0 ~ FC,
                            log2FoldChange < 0 ~ -(1 / FC))) %>%
  dplyr::select(gene, log2FoldChange, padj, FC, FC_sig, pvalue) %>%
  drop_na(padj) %>%
  arrange(desc(abs(log2FoldChange)))

# Write top DEGs to files
DEG_up <- res2 %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0)
write.table(head(DEG_up$gene, 250), file = "/data/DEG_up.csv", row.names = FALSE)

DEG_down <- res2 %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0)
write.table(head(DEG_down$gene, 250), file = "/data/DEG_down.csv", row.names = FALSE)

# Create volcano plot
res_data <- data.frame(res)
res_data$gene_symbol <- row.names(res_data)
res_data$diffexpressed <- "NO"
res_data$diffexpressed[res_data$log2FoldChange > 0.6 & res_data$pvalue < 0.05] <- "UP"
res_data$diffexpressed[res_data$log2FoldChange < -0.6 & res_data$pvalue < 0.05] <- "DOWN"
res_data$delabel <- ifelse(res_data$gene_symbol %in% head(res_data[order(res_data$padj), "gene_symbol"], 20), res_data$gene_symbol, NA)
myvolcanoplot <- ggplot(data = res_data, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 75), xlim = c(-8, 8)) +
  labs(color = 'Expression Change',
       x = expression("log"[2]*"Fold Change"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  ggtitle('Lifelong Non-smoker and Current reformed smoker for > 15 years vs Current smoker and Current reformed smoker for ≤15 years') +
  geom_text_repel(max.overlaps = Inf)
myvolcanoplot
