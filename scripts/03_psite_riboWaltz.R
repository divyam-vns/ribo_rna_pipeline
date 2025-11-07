#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
cfg <- yaml::yaml.load_file(args[1])
library(riboWaltz); library(GenomicAlignments); library(GenomicFeatures); library(tidyverse)
samples <- read_tsv(cfg$sample_table)
proj <- cfg$project_dir
bams <- list()
for(i in 1:nrow(samples)){
  base <- sub(".fastq.gz$","", basename(samples$ribo_fastq[i]))
  bams[[samples$sample_id[i]]] <- file.path(proj,"step02_star", base, "Aligned.sortedByCoord.out.bam")
}
txdb <- makeTxDbFromGFF(cfg$gtf, format="gtf")
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
reads_list <- list()
for(s in names(bams)){
  gal <- readGAlignments(bams[[s]])
  df <- data.frame(seqnames=as.character(seqnames(gal)), start=start(gal), end=end(gal), strand=as.character(strand(gal)))
  df$length <- width(gal)
  reads_list[[s]] <- df
}
psites <- psite(reads_list, cds)
pdf(file.path(proj,"step03_riboQC","ribo_qc.pdf"))
length_dist(reads_list); frame_psite(psites)
dev.off()
psite_df <- psite_df(psites)
dir.create(file.path(proj,"step03_riboQC"), recursive=TRUE, showWarnings=FALSE)
write_tsv(psite_df, file.path(proj,"step03_riboQC","psite_table.tsv"))
cat("P-site assignment done\n")
