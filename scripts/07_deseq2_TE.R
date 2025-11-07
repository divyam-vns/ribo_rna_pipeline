#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
cfg <- yaml::yaml.load_file(args[1])
library(DESeq2); library(data.table); library(tidyverse)
proj <- cfg$project_dir
rna_f <- file.path(proj,"step06_counts","rna_counts.txt")
ribo_f<- file.path(proj,"step06_counts","ribo_counts.txt")
samples <- read_tsv(cfg$sample_table)
read_counts <- function(f){
  df <- fread(f, skip=1, header=TRUE)
  colnames(df)[1] <- "GeneID"
  return(df)
}
rna_df <- read_counts(rna_f); ribo_df <- read_counts(ribo_f)
mat_rna <- as.matrix(rna_df[,7:ncol(rna_df)])
rownames(mat_rna) <- rna_df$GeneID
mat_ribo<- as.matrix(ribo_df[,7:ncol(ribo_df)])
rownames(mat_ribo) <- ribo_df$GeneID
common <- intersect(rownames(mat_ribo), rownames(mat_rna))
mat_rna <- mat_rna[common,]; mat_ribo<- mat_ribo[common,]
colnames(mat_rna) <- paste0(samples$sample_id, "_RNA")
colnames(mat_ribo)<- paste0(samples$sample_id, "_Ribo")
counts_comb <- cbind(mat_ribo, mat_rna)
coldata <- data.frame(sample=colnames(counts_comb),
                      condition=rep(samples$condition,2),
                      assay=rep(c("Ribo","RNA"), each=nrow(samples)))
rownames(coldata) <- coldata$sample
dds <- DESeqDataSetFromMatrix(counts_comb, colData=coldata, design=~condition + assay + condition:assay)
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)
int_name <- grep(":", resultsNames(dds), value=TRUE)
res_int <- results(dds, name=int_name)
outdir <- file.path(proj,"step07_deseq2")
dir.create(outdir, showWarnings=FALSE)
write_tsv(as.data.frame(res_int) %>% rownames_to_column("GeneID"), file.path(outdir,"differential_TE_results.tsv"))
# compute TE per sample
norm <- counts(dds, normalized=TRUE)
ribo_norm <- norm[, grepl("_Ribo", colnames(norm))]
rna_norm  <- norm[, grepl("_RNA", colnames(norm))]
# ensure matched order
samples_ids <- samples$sample_id
ribo_cols <- paste0(samples_ids,"_Ribo")
rna_cols  <- paste0(samples_ids,"_RNA")
te_mat <- log2(ribo_norm[,ribo_cols] + 1) - log2(rna_norm[,rna_cols] + 1)
write_tsv(as.data.frame(te_mat) %>% rownames_to_column("GeneID"), file.path(outdir,"TE_log2_per_sample.tsv"))
cat("DESeq2 TE analysis done. Outputs in", outdir, "\n")
