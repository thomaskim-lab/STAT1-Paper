############################################################
### 00_bulk_setup.R
### Shared setup for STAT1 paper bulk RNA-seq analyses
############################################################

### Clear workspace if desired
# rm(list = ls())

############################################################
### 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(openxlsx)
  library(readr)
  library(readxl)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(VennDiagram)
  library(ggpubr)

  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)

  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(fgsea)
  library(enrichplot)
})

############################################################
### 2. Define project paths
############################################################

# This assumes scripts are run from the Bulk_RNA-Seq folder.
project_dir <- normalizePath("..")

metadata_dir <- file.path(project_dir, "metadata")
gene_sets_dir <- file.path(project_dir, "gene_sets")
results_dir <- file.path(project_dir, "results", "bulk_rnaseq")
session_info_dir <- file.path(project_dir, "session_info")

dir.create(metadata_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gene_sets_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(session_info_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
### 3. Helper functions
############################################################

strip_ensembl_version <- function(x) {
  sub("\\..*$", "", x)
}

make_integer_counts <- function(count_mat) {
  mode(count_mat) <- "numeric"
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
  count_mat
}

filter_low_count_genes <- function(count_mat, min_count = 10, min_samples = 2) {
  count_mat[rowSums(count_mat >= min_count) >= min_samples, , drop = FALSE]
}

collapse_to_gene_symbol <- function(vsd_mat, gene_symbols) {
  stopifnot(length(gene_symbols) == nrow(vsd_mat))

  df <- as.data.frame(vsd_mat)
  df$gene_symbol <- gene_symbols
  df$mean_expr <- rowMeans(vsd_mat, na.rm = TRUE)

  df <- df %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    arrange(gene_symbol, desc(mean_expr)) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  rownames(df) <- df$gene_symbol
  df$gene_symbol <- NULL
  df$mean_expr <- NULL

  as.matrix(df)
}

save_csv <- function(x, filename) {
  write.csv(x, file.path(results_dir, filename), row.names = FALSE)
}

save_table <- function(x, filename) {
  write.table(
    x,
    file = file.path(results_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}