############################################################
### INPUT BLOCK FOR MATRIX WITH *_cpm AND *_count COLUMNS
############################################################

library(tidyverse)
library(DESeq2)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggrepel)

# Input files
setwd("C:/Users/au731006/Downloads")

base_dir <- "C:/Users/au731006/Downloads/Plasmidsaurus"

matrix_files <- c(
  "5SZCYD" = file.path(base_dir, "Pilot", "5SZCYD_results", "5SZCYD-expression-matrix.tsv"),
  "BN4GZZ" = file.path(base_dir, "Mar23", "BN4GZZ_results", "BN4GZZ-expression-matrix.tsv")
)

sample_map_files <- c(
  "5SZCYD" = file.path(base_dir, "Pilot", "5SZCYD_results", "5SZCYD.csv"),
  "BN4GZZ" = file.path(base_dir, "Mar23", "BN4GZZ_results", "BN4GZZ.csv")
)

outdir <- file.path(base_dir, "Merged_STAT1_siRNA_5SZCYD_BN4GZZ")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
### 1-5. Read, parse, rename, and merge 5SZCYD + BN4GZZ
############################################################

# Input files
matrix_files <- c(
  "5SZCYD" = "C:/Users/au731006/Downloads/Plasmidsaurus/Pilot/5SZCYD_results/5SZCYD-expression-matrix.tsv",
  "BN4GZZ" = "C:/Users/au731006/Downloads/Plasmidsaurus/Mar23/BN4GZZ_results/BN4GZZ-expression-matrix.tsv"
)

sample_map_files <- c(
  "5SZCYD" = "C:/Users/au731006/Downloads/Plasmidsaurus/Pilot/5SZCYD_results/5SZCYD.csv",
  "BN4GZZ" = "C:/Users/au731006/Downloads/Plasmidsaurus/Mar23/BN4GZZ_results/BN4GZZ.csv"
)

read_one_run <- function(matrix_file, sample_map_file, batch_name) {
  
  ############################################################
  ### Read matrix and sample map
  ############################################################
  
  expr_df <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
  sample_map <- read.csv(sample_map_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  # expected columns in sample_map:
  # run_id, Sample_name
  stopifnot(all(c("run_id", "Sample_name") %in% colnames(sample_map)))
  stopifnot(all(c("gene_id", "gene_name", "gene_biotype") %in% colnames(expr_df)))
  
  ############################################################
  ### Extract count and CPM columns separately
  ############################################################
  
  count_cols <- grep("_count$", colnames(expr_df), value = TRUE)
  cpm_cols   <- grep("_cpm$", colnames(expr_df), value = TRUE)
  
  if (length(count_cols) == 0) {
    stop("No *_count columns found in the expression matrix for batch: ", batch_name)
  }
  
  # function to strip suffix
  strip_suffix <- function(x) {
    sub("_(count|cpm)$", "", x)
  }
  
  count_run_ids <- strip_suffix(count_cols)
  cpm_run_ids   <- strip_suffix(cpm_cols)
  
  ############################################################
  ### Build count matrix
  ############################################################
  
  count_mat <- as.matrix(expr_df[, count_cols, drop = FALSE])
  rownames(count_mat) <- sub("\\..*$", "", expr_df$gene_id)
  
  # rename columns from 5SZCYD_x_count / BN4GZZ_x_count to 5SZCYD_x / BN4GZZ_x
  colnames(count_mat) <- count_run_ids
  
  # remove duplicated genes if any
  count_mat <- count_mat[!duplicated(rownames(count_mat)), , drop = FALSE]
  
  # DESeq2 requires integer counts
  mode(count_mat) <- "numeric"
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
  
  ############################################################
  ### Optional: build CPM matrix too
  ############################################################
  
  cpm_mat <- NULL
  
  if (length(cpm_cols) > 0) {
    cpm_mat <- as.matrix(expr_df[, cpm_cols, drop = FALSE])
    rownames(cpm_mat) <- sub("\\..*$", "", expr_df$gene_id)
    colnames(cpm_mat) <- cpm_run_ids
    cpm_mat <- cpm_mat[!duplicated(rownames(cpm_mat)), , drop = FALSE]
    mode(cpm_mat) <- "numeric"
  }
  
  ############################################################
  ### Match run IDs to real sample names
  ############################################################
  
  # Check all count columns have mapping
  missing_map <- setdiff(colnames(count_mat), sample_map$run_id)
  if (length(missing_map) > 0) {
    stop(
      "These run IDs are missing from sample_map for batch ",
      batch_name,
      ": ",
      paste(missing_map, collapse = ", ")
    )
  }
  
  # Reorder sample_map to match matrix columns
  idx <- match(colnames(count_mat), sample_map$run_id)
  
  if (any(is.na(idx))) {
    stop(
      "Failed to match these count_mat columns for batch ",
      batch_name,
      ": ",
      paste(colnames(count_mat)[is.na(idx)], collapse = ", ")
    )
  }
  
  sample_map <- sample_map[idx, , drop = FALSE]
  
  stopifnot(all(sample_map$run_id == colnames(count_mat)))
  
  # Rename matrices to real sample names
  colnames(count_mat) <- sample_map$Sample_name
  
  if (!is.null(cpm_mat)) {
    missing_cpm_map <- setdiff(colnames(cpm_mat), sample_map$run_id)
    
    if (length(missing_cpm_map) == 0) {
      cpm_idx <- match(sample_map$run_id, colnames(cpm_mat))
      cpm_mat <- cpm_mat[, cpm_idx, drop = FALSE]
      colnames(cpm_mat) <- sample_map$Sample_name
    }
  }
  
  # Add batch/run information
  sample_map$Batch <- batch_name
  
  rownames(sample_map) <- sample_map$Sample_name
  
  return(list(
    count_mat = count_mat,
    cpm_mat = cpm_mat,
    sample_map = sample_map
  ))
}

############################################################
### Read both runs
############################################################

run_5SZCYD <- read_one_run(
  matrix_file = matrix_files["5SZCYD"],
  sample_map_file = sample_map_files["5SZCYD"],
  batch_name = "5SZCYD"
)

run_BN4GZZ <- read_one_run(
  matrix_file = matrix_files["BN4GZZ"],
  sample_map_file = sample_map_files["BN4GZZ"],
  batch_name = "BN4GZZ"
)

############################################################
### Merge count matrices
############################################################

common_genes <- Reduce(
  intersect,
  list(
    rownames(run_5SZCYD$count_mat),
    rownames(run_BN4GZZ$count_mat)
  )
)

count_mat <- cbind(
  run_5SZCYD$count_mat[common_genes, , drop = FALSE],
  run_BN4GZZ$count_mat[common_genes, , drop = FALSE]
)

############################################################
### Merge CPM matrices, if present
############################################################

cpm_mat <- NULL

if (!is.null(run_5SZCYD$cpm_mat) && !is.null(run_BN4GZZ$cpm_mat)) {
  common_cpm_genes <- Reduce(
    intersect,
    list(
      rownames(run_5SZCYD$cpm_mat),
      rownames(run_BN4GZZ$cpm_mat)
    )
  )
  
  cpm_mat <- cbind(
    run_5SZCYD$cpm_mat[common_cpm_genes, , drop = FALSE],
    run_BN4GZZ$cpm_mat[common_cpm_genes, , drop = FALSE]
  )
}

############################################################
### Merge sample maps
############################################################

sample_map <- bind_rows(
  run_5SZCYD$sample_map,
  run_BN4GZZ$sample_map
)

sample_map <- sample_map[match(colnames(count_mat), sample_map$Sample_name), , drop = FALSE]
rownames(sample_map) <- sample_map$Sample_name

stopifnot(all(colnames(count_mat) == rownames(sample_map)))

############################################################
### Optional: exclude problematic sample 5SZCYD_1 / HMC3_Lipo_24h_R1
############################################################

exclude_samples <- c("HMC3_Lipo_24h_R1")

keep_samples <- setdiff(colnames(count_mat), exclude_samples)

count_mat <- count_mat[, keep_samples, drop = FALSE]
sample_map <- sample_map[keep_samples, , drop = FALSE]

if (!is.null(cpm_mat)) {
  cpm_mat <- cpm_mat[, keep_samples, drop = FALSE]
}

stopifnot(all(colnames(count_mat) == rownames(sample_map)))


############################################################
### 6. Parse metadata from sample names
############################################################


parse_sample_metadata <- function(sample_names, sample_map) {
  tibble(Sample_name = sample_names) %>%
    left_join(
      sample_map %>% 
        dplyr::select(Sample_name, Batch),
      by = "Sample_name"
    ) %>%
    mutate(
      Condition = case_when(
        grepl("^HMC3_Lipo_24h_R[0-9]+$", Sample_name) ~ "Lipo_only",
        grepl("^HMC3_siSTAT1_1_24h_R[0-9]+$", Sample_name) ~ "siSTAT1_1",
        grepl("^HMC3_siSTAT1_3_24h_R[0-9]+$", Sample_name) ~ "siSTAT1_3",
        TRUE ~ "Other"
      ),
      siRNA = case_when(
        Condition == "Lipo_only" ~ "Control",
        Condition == "siSTAT1_1" ~ "siSTAT1",
        Condition == "siSTAT1_3" ~ "siSTAT1",
        TRUE ~ "Other"
      ),
      Duplex = case_when(
        Condition == "Lipo_only" ~ "None",
        Condition == "siSTAT1_1" ~ "DsiRNA_1",
        Condition == "siSTAT1_3" ~ "DsiRNA_3",
        TRUE ~ "Other"
      ),
      Replicate = sub(".*_R([0-9]+)$", "R\\1", Sample_name)
    ) %>%
    mutate(
      Batch = factor(Batch, levels = c("5SZCYD", "BN4GZZ")),
      Condition = factor(
        Condition,
        levels = c("Lipo_only", "siSTAT1_1", "siSTAT1_3")
      ),
      siRNA = factor(siRNA, levels = c("Control", "siSTAT1")),
      Duplex = factor(Duplex, levels = c("None", "DsiRNA_1", "DsiRNA_3")),
      Replicate = factor(Replicate)
    )
}

meta <- parse_sample_metadata(colnames(count_mat), sample_map)
rownames(meta) <- meta$Sample_name

stopifnot(all(colnames(count_mat) == rownames(meta)))

print(meta)
print(table(meta$Batch, meta$Condition))

if (any(is.na(meta$Condition)) || any(meta$Condition == "Other")) {
  stop("Some samples were parsed incorrectly. Check sample names.")
}


############################################################
### 7. Save parsed metadata
############################################################

write.xlsx(
  meta,
  file = file.path(outdir, "Parsed_sample_metadata.xlsx"),
  overwrite = TRUE
)


############################################################
### 8. Build DESeq2 object - batch-adjusted 3-condition model
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ Batch + Condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)

# Inspect available coefficients
resultsNames(dds)


############################################################
### 9. Annotation
############################################################

gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(dds),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

############################################################
### 10. Save count/CPM matrices with real sample names
############################################################

write.xlsx(
  as.data.frame(count_mat) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "Counts_matrix_renamed.xlsx"),
  overwrite = TRUE
)

if (!is.null(cpm_mat)) {
  write.xlsx(
    as.data.frame(cpm_mat) %>% rownames_to_column("ENSEMBL"),
    file = file.path(outdir, "CPM_matrix_renamed.xlsx"),
    overwrite = TRUE
  )
}

############################################################
### 11. Batch-adjusted DE analysis: STAT1 siRNA trial
############################################################

make_res_table <- function(res_obj, annot_tbl) {
  as.data.frame(res_obj) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(annot_tbl, by = "ENSEMBL") %>%
    arrange(padj)
}

# Main contrasts
res_siSTAT1_1_vs_Lipo <- results(
  dds,
  contrast = c("Condition", "siSTAT1_1", "Lipo_only")
)

res_siSTAT1_3_vs_Lipo <- results(
  dds,
  contrast = c("Condition", "siSTAT1_3", "Lipo_only")
)

res_siSTAT1_3_vs_1 <- results(
  dds,
  contrast = c("Condition", "siSTAT1_3", "siSTAT1_1")
)

# Annotated result tables
res_siSTAT1_1_vs_Lipo_df <- make_res_table(res_siSTAT1_1_vs_Lipo, gene_annot)
res_siSTAT1_3_vs_Lipo_df <- make_res_table(res_siSTAT1_3_vs_Lipo, gene_annot)
res_siSTAT1_3_vs_1_df    <- make_res_table(res_siSTAT1_3_vs_1, gene_annot)

# Save DE results
write.xlsx(
  list(
    siSTAT1_1_vs_Lipo = res_siSTAT1_1_vs_Lipo_df,
    siSTAT1_3_vs_Lipo = res_siSTAT1_3_vs_Lipo_df,
    siSTAT1_3_vs_siSTAT1_1 = res_siSTAT1_3_vs_1_df
  ),
  file = file.path(outdir, "DE_results_STAT1_siRNA_batch_adjusted.xlsx"),
  overwrite = TRUE
)


############################################################
### 12. Optional pooled siSTAT1 analysis
############################################################

# This tests both STAT1 siRNAs together against Lipo-only.
# Useful as a broad knockdown signature, but keep individual siRNA results primary.

meta_pooled <- meta
meta_pooled$siRNA <- factor(meta_pooled$siRNA, levels = c("Control", "siSTAT1"))

dds_pooled <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_pooled,
  design = ~ Batch + siRNA
)

keep_pooled <- rowSums(counts(dds_pooled) >= 10) >= 2
dds_pooled <- dds_pooled[keep_pooled, ]

dds_pooled <- DESeq(dds_pooled)

res_siSTAT1_pooled_vs_Lipo <- results(
  dds_pooled,
  contrast = c("siRNA", "siSTAT1", "Control")
)

res_siSTAT1_pooled_vs_Lipo_df <- make_res_table(
  res_siSTAT1_pooled_vs_Lipo,
  gene_annot
)

write.xlsx(
  list(
    siSTAT1_pooled_vs_Lipo = res_siSTAT1_pooled_vs_Lipo_df
  ),
  file = file.path(outdir, "DE_results_STAT1_siRNA_pooled_batch_adjusted.xlsx"),
  overwrite = TRUE
)

############################################################
### 13. PCA plot
############################################################

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c("Condition", "Batch"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    color = Condition,
    shape = Batch,
    label = name
  )
) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - HMC3 STAT1 siRNA 24h, batch-adjusted design") +
  theme_bw()

ggsave(
  file.path(outdir, "PCA_HMC3_STAT1_siRNA_24h.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
### 14. Session summary
############################################################

message("Analysis completed successfully.")
message("Output directory: ", outdir)
message("Samples processed: ", ncol(count_mat))
message("Genes retained in batch-adjusted model: ", nrow(dds))

# Required libraries for downstream analysis
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(VennDiagram)
library(grid)
library(ReactomePA)
library(msigdbr)
library(fgsea)
library(enrichplot)
library(patchwork)

############################################################
### Define contrast list for this experiment
############################################################

# These result data frames should already exist from section 11/12:
# res_siSTAT1_1_vs_Lipo_df
# res_siSTAT1_3_vs_Lipo_df
# res_siSTAT1_3_vs_1_df
# Optional from pooled analysis:
# res_siSTAT1_pooled_vs_Lipo_df

contrast_results <- list(
  siSTAT1_1_vs_Lipo = res_siSTAT1_1_vs_Lipo_df,
  siSTAT1_3_vs_Lipo = res_siSTAT1_3_vs_Lipo_df,
  siSTAT1_3_vs_siSTAT1_1 = res_siSTAT1_3_vs_1_df
)

if (exists("res_siSTAT1_pooled_vs_Lipo_df")) {
  contrast_results[["siSTAT1_pooled_vs_Lipo"]] <- res_siSTAT1_pooled_vs_Lipo_df
}

############################################################
### 15. Variance stabilized matrix and normalized counts
############################################################

# Use the batch-adjusted dds object created with:
# design = ~ Batch + Condition
vsd <- vst(dds, blind = FALSE)

norm_counts <- counts(dds, normalized = TRUE)

write.xlsx(
  as.data.frame(norm_counts) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "Normalized_counts.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  as.data.frame(assay(vsd)) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "VST_matrix.xlsx"),
  overwrite = TRUE
)

############################################################
### 16. Sample distance heatmap
############################################################

sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

annotation_df <- as.data.frame(colData(dds)[, c("Condition", "Batch"), drop = FALSE])
rownames(annotation_df) <- colnames(dds)

png(file.path(outdir, "Sample_distance_heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sample_dist_mat,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "Sample-to-sample distances"
)
dev.off()

############################################################
### 17. Volcano plot function
############################################################

plot_volcano <- function(res_df, title_text, label_top_n = 15) {
  df <- res_df %>%
    mutate(
      sig = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 0.5  ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -0.5 ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(ifelse(is.na(pvalue) | pvalue <= 0, NA, pvalue))
    )
  
  top_labels <- df %>%
    filter(sig != "NS", !is.na(SYMBOL)) %>%
    arrange(padj) %>%
    slice_head(n = label_top_n)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10p, color = sig)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = top_labels,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      x = "log2 fold change",
      y = "-log10(p-value)",
      color = ""
    )
  
  return(p)
}

volcano_titles <- c(
  siSTAT1_1_vs_Lipo = "STAT1 DsiRNA #1 vs Lipofectamine-only",
  siSTAT1_3_vs_Lipo = "STAT1 DsiRNA #3 vs Lipofectamine-only",
  siSTAT1_3_vs_siSTAT1_1 = "STAT1 DsiRNA #3 vs STAT1 DsiRNA #1",
  siSTAT1_pooled_vs_Lipo = "Pooled STAT1 siRNA vs Lipofectamine-only"
)

for (nm in names(contrast_results)) {
  title_text <- ifelse(nm %in% names(volcano_titles), volcano_titles[[nm]], nm)
  p <- plot_volcano(contrast_results[[nm]], title_text)
  
  ggsave(
    file.path(outdir, paste0("Volcano_", nm, ".png")),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    file.path(outdir, paste0("Volcano_", nm, ".eps")),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

############################################################
### 18. Top DEGs
############################################################

get_top_genes <- function(res_df, n = 10) {
  top_up <- res_df %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = n)
  
  top_down <- res_df %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < -1) %>%
    arrange(log2FoldChange) %>%
    slice_head(n = n)
  
  list(up = top_up, down = top_down)
}

top_deg_list <- list()

for (nm in names(contrast_results)) {
  tops <- get_top_genes(contrast_results[[nm]], n = 10)
  top_deg_list[[paste0(nm, "_up")]] <- tops$up
  top_deg_list[[paste0(nm, "_down")]] <- tops$down
}

write.xlsx(
  top_deg_list,
  file = file.path(outdir, "Top_DEGs.xlsx"),
  overwrite = TRUE
)

############################################################
### 19. Heatmap of selected STAT1 / IFN / inflammatory genes
############################################################

selected_genes <- c(
  "STAT1", "STAT2", "IRF1", "IRF7",
  "ISG15", "IFIT1", "IFIT3", "MX1", "MX2",
  "NFKBIA", "CXCL8", "CCL2", "TNFAIP3"
)

sel_annot <- gene_annot %>%
  filter(SYMBOL %in% selected_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(sel_annot) > 0) {
  keep_genes <- intersect(sel_annot$ENSEMBL, rownames(assay(vsd)))
  sel_mat <- assay(vsd)[keep_genes, , drop = FALSE]
  rownames(sel_mat) <- sel_annot$SYMBOL[match(rownames(sel_mat), sel_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Condition, Batch, Replicate) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(sel_mat))
  sel_mat <- sel_mat[, ordered_samples, drop = FALSE]
  
  annotation_df2 <- meta[, c("Sample_name", "Condition", "Batch"), drop = FALSE] %>%
    as.data.frame()
  rownames(annotation_df2) <- annotation_df2$Sample_name
  annotation_df2$Sample_name <- NULL
  annotation_df2 <- annotation_df2[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Selected_STAT1_IFN_genes_heatmap.png"), width = 1200, height = 1000, res = 150)
  pheatmap(
    sel_mat,
    scale = "row",
    annotation_col = annotation_df2,
    cluster_cols = FALSE,
    main = "Selected STAT1 / IFN / inflammatory genes"
  )
  dev.off()
}

############################################################
### 20. GO enrichment
############################################################

run_go <- function(res_df, direction = c("up", "down")) {
  direction <- match.arg(direction)
  
  if (direction == "up") {
    genes <- res_df %>%
      filter(!is.na(padj), padj < 0.05, log2FoldChange > 0.5, !is.na(ENTREZID)) %>%
      pull(ENTREZID) %>%
      unique()
  } else {
    genes <- res_df %>%
      filter(!is.na(padj), padj < 0.05, log2FoldChange < -0.5, !is.na(ENTREZID)) %>%
      pull(ENTREZID) %>%
      unique()
  }
  
  if (length(genes) < 5) return(NULL)
  
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

go_list <- list()
go_objects <- list()

for (nm in names(contrast_results)) {
  ego_up <- run_go(contrast_results[[nm]], "up")
  ego_down <- run_go(contrast_results[[nm]], "down")
  
  go_objects[[paste0(nm, "_up")]] <- ego_up
  go_objects[[paste0(nm, "_down")]] <- ego_down
  
  if (!is.null(ego_up)) {
    go_list[[paste0("GO_", nm, "_up")]] <- as.data.frame(ego_up)
  }
  if (!is.null(ego_down)) {
    go_list[[paste0("GO_", nm, "_down")]] <- as.data.frame(ego_down)
  }
}

if (length(go_list) > 0) {
  write.xlsx(go_list, file = file.path(outdir, "GO_enrichment_results.xlsx"), overwrite = TRUE)
}

save_go_plot <- function(go_obj, filename, title_text) {
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj)) > 0) {
    p <- dotplot(go_obj, showCategory = 10) + ggtitle(title_text)
    ggsave(file.path(outdir, filename), p, width = 8, height = 6, dpi = 300)
  }
}

for (nm in names(go_objects)) {
  save_go_plot(
    go_objects[[nm]],
    paste0("GO_", nm, ".png"),
    paste0("GO BP: ", nm)
  )
  save_go_plot(
    go_objects[[nm]],
    paste0("GO_", nm, ".eps"),
    paste0("GO BP: ", nm)
  )
}

############################################################
### 21. Venn diagram of significant DEGs
############################################################

sig_gene_sets <- lapply(contrast_results, function(res_df) {
  res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
    pull(SYMBOL) %>%
    unique()
})

# VennDiagram supports 2-5 sets. This will run for 3 or 4 contrasts.
if (length(sig_gene_sets) >= 2 && length(sig_gene_sets) <= 5) {
  venn.plot <- venn.diagram(
    x = sig_gene_sets,
    filename = NULL,
    fill = c("purple", "green3", "red", "orange", "skyblue")[seq_along(sig_gene_sets)],
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.1,
    main = "DEG overlap across STAT1 siRNA contrasts"
  )
  
  png(file.path(outdir, "DEG_venn.png"), width = 1400, height = 1200, res = 150)
  grid.draw(venn.plot)
  dev.off()
}

common_genes_all <- if (length(sig_gene_sets) >= 2) Reduce(intersect, sig_gene_sets) else character(0)

overlap_export <- c(
  list(common_genes_all = data.frame(Gene = common_genes_all)),
  lapply(sig_gene_sets, function(x) data.frame(Gene = x))
)

write.xlsx(
  overlap_export,
  file = file.path(outdir, "DEG_overlap_lists.xlsx"),
  overwrite = TRUE
)

############################################################
### 22. Focused lipid / LD / IFN modules
############################################################

lipid_modules <- list(
  LD_structure = c("PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5", "CIDEC", "CIDEA"),
  Lipolysis = c("PNPLA2", "ABHD5", "LIPE", "MGLL"),
  Esterification = c("ACAT1", "ACAT2", "SOAT1", "SOAT2", "DGAT1", "DGAT2", "ACSL3"),
  FA_activation_FAO = c("ACSL1", "ACSL3", "CPT1A", "CPT2", "ACADM", "ACADVL", "PPARGC1A"),
  Lysosome_lipophagy = c("LIPA", "LAMP1", "LAMP2", "CTSD", "RAB7A", "ATG5", "ATG7", "SQSTM1"),
  ER_stress = c("DDIT3", "ATF4", "XBP1", "HSPA5", "ERN1"),
  IFN_STAT_IRF = c("STAT1", "STAT2", "IRF1", "IRF7", "ISG15", "IFIT1", "IFIT3", "MX1", "MX2"),
  NFkB_inflammation = c("NFKBIA", "CXCL8", "CCL2", "TNFAIP3", "RELB", "JUNB")
)

extract_module_results <- function(res_df, module_list, contrast_label) {
  bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    res_df %>%
      filter(SYMBOL %in% genes) %>%
      mutate(Module = mod, Contrast = contrast_label)
  }))
}

focused_module_results <- lapply(names(contrast_results), function(nm) {
  extract_module_results(contrast_results[[nm]], lipid_modules, nm)
})
names(focused_module_results) <- names(contrast_results)

write.xlsx(
  focused_module_results,
  file = file.path(outdir, "Focused_lipid_STAT1_IFN_gene_results.xlsx"),
  overwrite = TRUE
)

############################################################
### 23. Focused heatmap of lipid / LD / stress / IFN genes
############################################################

all_lipid_genes <- unique(unlist(lipid_modules))

lipid_annot <- gene_annot %>%
  filter(SYMBOL %in% all_lipid_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(lipid_annot) > 0) {
  keep_genes <- intersect(lipid_annot$ENSEMBL, rownames(assay(vsd)))
  lipid_mat <- assay(vsd)[keep_genes, , drop = FALSE]
  rownames(lipid_mat) <- lipid_annot$SYMBOL[match(rownames(lipid_mat), lipid_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Condition, Batch, Replicate) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(lipid_mat))
  lipid_mat <- lipid_mat[, ordered_samples, drop = FALSE]
  
  anno_lipid <- meta[, c("Sample_name", "Condition", "Batch"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_lipid) <- anno_lipid$Sample_name
  anno_lipid$Sample_name <- NULL
  anno_lipid <- anno_lipid[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Focused_lipid_STAT1_IFN_heatmap.png"), width = 1600, height = 1300, res = 160)
  pheatmap(
    lipid_mat,
    scale = "row",
    annotation_col = anno_lipid,
    cluster_cols = FALSE,
    main = "Focused lipid / LD / lysosome / stress / IFN genes"
  )
  dev.off()
  
  cairo_ps(
    file = file.path(outdir, "Focused_lipid_STAT1_IFN_heatmap.eps"),
    width = 10,
    height = 8,
    fallback_resolution = 300
  )
  pheatmap(
    lipid_mat,
    scale = "row",
    annotation_col = anno_lipid,
    cluster_cols = FALSE,
    main = "Focused lipid / LD / lysosome / stress / IFN genes"
  )
  dev.off()
}

############################################################
### 24. Module-level summary
############################################################

summarize_modules <- function(res_df, module_list, contrast_label) {
  bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    sub <- res_df %>% filter(SYMBOL %in% genes)
    
    data.frame(
      Module = mod,
      Contrast = contrast_label,
      n_genes_found = nrow(sub),
      mean_log2FC = ifelse(nrow(sub) > 0, mean(sub$log2FoldChange, na.rm = TRUE), NA),
      median_log2FC = ifelse(nrow(sub) > 0, median(sub$log2FoldChange, na.rm = TRUE), NA),
      n_up_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange > 1), 0),
      n_down_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange < -1), 0)
    )
  }))
}

module_summary_all <- bind_rows(lapply(names(contrast_results), function(nm) {
  summarize_modules(contrast_results[[nm]], lipid_modules, nm)
}))

write.xlsx(
  module_summary_all,
  file = file.path(outdir, "Lipid_STAT1_IFN_module_summary.xlsx"),
  overwrite = TRUE
)

p_module_mean <- ggplot(module_summary_all, aes(x = Module, y = mean_log2FC, fill = Contrast)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Lipid / LD / IFN module summary", x = "", y = "Mean log2FC")

ggsave(
  file.path(outdir, "Lipid_STAT1_IFN_module_summary_barplot.png"),
  p_module_mean,
  width = 10,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Lipid_STAT1_IFN_module_summary_barplot.eps"),
  p_module_mean,
  width = 10,
  height = 6.5,
  device = cairo_ps
)

p_module_sig_up <- ggplot(module_summary_all, aes(x = Module, y = n_up_sig, fill = Contrast)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Significant upregulated genes per module", x = "", y = "Count")

ggsave(
  file.path(outdir, "Lipid_STAT1_IFN_module_sig_up_barplot.png"),
  p_module_sig_up,
  width = 10,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Lipid_STAT1_IFN_module_sig_up_barplot.eps"),
  p_module_sig_up,
  width = 10,
  height = 6.5,
  device = cairo_ps
)

############################################################
### 25. STAT1 / lipid axis focused export
############################################################

custom_qpcr_panel <- c(
  "PLIN2", "PLIN3", "ACAT1", "ACAT2", "LPIN1", "PNPLA2", "LIPE",
  "ABHD5", "LIPA", "CPT1A", "ACSL1", "ACSL3", "PPARGC1A",
  "DDIT3", "LAMP1", "STAT1", "STAT2", "IRF1", "IRF7",
  "ISG15", "IFIT1", "IFIT3", "MX1", "MX2", "NFKBIA", "CXCL8"
)

extract_qpcr_panel <- function(res_df, genes, contrast_label) {
  res_df %>%
    filter(SYMBOL %in% genes) %>%
    dplyr::select(ENSEMBL, SYMBOL, GENENAME, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(match(SYMBOL, genes)) %>%
    mutate(Contrast = contrast_label)
}

qpcr_panel_list <- lapply(names(contrast_results), function(nm) {
  extract_qpcr_panel(contrast_results[[nm]], custom_qpcr_panel, nm)
})
names(qpcr_panel_list) <- names(contrast_results)

write.xlsx(
  qpcr_panel_list,
  file = file.path(outdir, "Custom_lipid_STAT1_IFN_qPCR_panel.xlsx"),
  overwrite = TRUE
)

############################################################
### 26. Targeted panel heatmap
############################################################

panel_annot <- gene_annot %>%
  filter(SYMBOL %in% custom_qpcr_panel) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(panel_annot) > 0) {
  keep_genes <- intersect(panel_annot$ENSEMBL, rownames(assay(vsd)))
  panel_mat <- assay(vsd)[keep_genes, , drop = FALSE]
  rownames(panel_mat) <- panel_annot$SYMBOL[match(rownames(panel_mat), panel_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Condition, Batch, Replicate) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(panel_mat))
  panel_mat <- panel_mat[, ordered_samples, drop = FALSE]
  
  anno_panel <- meta[, c("Sample_name", "Condition", "Batch"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_panel) <- anno_panel$Sample_name
  anno_panel$Sample_name <- NULL
  anno_panel <- anno_panel[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Custom_lipid_STAT1_IFN_qPCR_panel_heatmap.png"), width = 1600, height = 1100, res = 160)
  pheatmap(
    panel_mat,
    scale = "row",
    annotation_col = anno_panel,
    cluster_cols = FALSE,
    main = "Custom lipid / STAT1 / IFN axis panel"
  )
  dev.off()
}

############################################################
### 27. Quick interpretation helper table
############################################################

interpret_panel_direction <- function(panel_df) {
  panel_df %>%
    mutate(
      Direction = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 1  ~ "Up_sig",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down_sig",
        !is.na(log2FoldChange) & log2FoldChange > 0 ~ "Up_nonsig",
        !is.na(log2FoldChange) & log2FoldChange < 0 ~ "Down_nonsig",
        TRUE ~ "NA"
      )
    )
}

qpcr_panel_interpretation_list <- lapply(qpcr_panel_list, interpret_panel_direction)

write.xlsx(
  qpcr_panel_interpretation_list,
  file = file.path(outdir, "Custom_lipid_STAT1_IFN_qPCR_panel_interpretation.xlsx"),
  overwrite = TRUE
)

############################################################
### 28. Reactome GSEA setup
############################################################

make_ranked_vector <- function(res_df, use_stat = TRUE) {
  if (use_stat && "stat" %in% colnames(res_df)) {
    ranked <- res_df %>%
      filter(!is.na(ENTREZID), !is.na(stat)) %>%
      group_by(ENTREZID) %>%
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(stat)) %>%
      dplyr::select(ENTREZID, stat)
    
    geneList <- ranked$stat
    names(geneList) <- ranked$ENTREZID
  } else {
    ranked <- res_df %>%
      filter(!is.na(ENTREZID), !is.na(log2FoldChange), !is.na(pvalue), pvalue > 0) %>%
      mutate(rank_score = sign(log2FoldChange) * -log10(pvalue)) %>%
      group_by(ENTREZID) %>%
      slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(rank_score)) %>%
      dplyr::select(ENTREZID, rank_score)
    
    geneList <- ranked$rank_score
    names(geneList) <- ranked$ENTREZID
  }
  
  geneList <- sort(geneList, decreasing = TRUE)
  return(geneList)
}

run_reactome_gsea <- function(geneList) {
  if (length(geneList) < 10) return(NULL)
  
  ReactomePA::gsePathway(
    geneList = geneList,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

gene_lists <- lapply(contrast_results, make_ranked_vector, use_stat = TRUE)
gsea_reactome_results <- lapply(gene_lists, run_reactome_gsea)

gsea_export <- list()
for (nm in names(gsea_reactome_results)) {
  obj <- gsea_reactome_results[[nm]]
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    gsea_export[[nm]] <- as.data.frame(obj)
  }
}

if (length(gsea_export) > 0) {
  write.xlsx(
    gsea_export,
    file = file.path(outdir, "Reactome_GSEA_results.xlsx"),
    overwrite = TRUE
  )
}

############################################################
### 29. Selected Reactome GSEA enrichment plots
############################################################

plot_selected_gseapaths <- function(gsea_obj, pattern_vec, prefix, outdir) {
  if (is.null(gsea_obj)) return(NULL)
  
  gsea_df <- as.data.frame(gsea_obj)
  if (nrow(gsea_df) == 0) return(NULL)
  
  for (pat in pattern_vec) {
    idx <- grep(pat, gsea_df$Description, ignore.case = TRUE)
    
    if (length(idx) > 0) {
      first_id <- gsea_df$ID[idx[1]]
      plot_title <- gsea_df$Description[idx[1]]
      
      p <- enrichplot::gseaplot2(
        gsea_obj,
        geneSetID = first_id,
        title = plot_title
      )
      
      safe_pat <- gsub("[^A-Za-z0-9]+", "_", pat)
      
      ggsave(
        filename = file.path(outdir, paste0(prefix, "_", safe_pat, ".png")),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(outdir, paste0(prefix, "_", safe_pat, ".eps")),
        plot = p,
        width = 8,
        height = 6,
        device = cairo_ps
      )
    }
  }
}

pattern_vec <- c(
  "interferon",
  "immune",
  "cytokine",
  "JAK",
  "STAT",
  "cholesterol",
  "fatty acid",
  "lysosome",
  "autophagy",
  "ER"
)

for (nm in names(gsea_reactome_results)) {
  plot_selected_gseapaths(
    gsea_obj = gsea_reactome_results[[nm]],
    pattern_vec = pattern_vec,
    prefix = paste0("Reactome_GSEA_", nm, "_selected"),
    outdir = outdir
  )
}

############################################################
### 30. Final add-on summary
############################################################

message("Extended downstream analysis completed successfully.")
message("Output directory: ", outdir)
message("Contrasts processed: ", paste(names(contrast_results), collapse = ", "))




