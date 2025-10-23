#!/usr/bin/env Rscript
library(DESeq2)
library(tidyverse)

# read config
config <- yaml::read_yaml("config.yaml")
samples <- config$samples
conds <- unlist(config$conditions)

# read counts into matrix
count_list <- lapply(samples, function(s){
    df <- read_tsv(paste0("counts/", s, ".tsv"), col_names = c("gene", s))
    df
})
counts_df <- reduce(count_list, full_join, by = "gene")
gene <- counts_df$gene
counts_mat <- as.matrix(counts_df[,-1])
rownames(counts_mat) <- gene
colnames(counts_mat) <- samples

coldata <- data.frame(
    sample = samples,
    condition = factor(unlist(config$conditions)[samples])
)
rownames(coldata) <- samples

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdf <- as.data.frame(res)
resdf$gene <- rownames(resdf)
write_tsv(resdf, "de/deseq2_results.tsv")
