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
matrix_file <- "C:/Users/au731006/Downloads/Plasmidsaurus/R4NCSL_results/R4NCSL-expression-matrix.tsv"
sample_map_file <- "C:/Users/au731006/Downloads/Plasmidsaurus/R4NCSL_results/R4NCSL.csv"
outdir <- "C:/Users/au731006/Downloads/Plasmidsaurus/Output/R4NCSL"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
### 1. Read matrix and sample map
############################################################

# Matrix is TSV, so use read.delim
expr_df <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
sample_map <- read.csv(sample_map_file, check.names = FALSE, stringsAsFactors = FALSE)

# expected columns in sample_map:
# run_id, Sample_name
stopifnot(all(c("run_id", "Sample_name") %in% colnames(sample_map)))
stopifnot(all(c("gene_id", "gene_name", "gene_biotype") %in% colnames(expr_df)))

############################################################
### 2. Extract count and CPM columns separately
############################################################

count_cols <- grep("_count$", colnames(expr_df), value = TRUE)
cpm_cols   <- grep("_cpm$", colnames(expr_df), value = TRUE)

if (length(count_cols) == 0) {
  stop("No *_count columns found in the expression matrix.")
}

# function to strip suffix
strip_suffix <- function(x) {
  sub("_(count|cpm)$", "", x)
}

count_run_ids <- strip_suffix(count_cols)
cpm_run_ids   <- strip_suffix(cpm_cols)

############################################################
### 3. Build count matrix
############################################################

count_mat <- as.matrix(expr_df[, count_cols, drop = FALSE])
rownames(count_mat) <- sub("\\..*$", "", expr_df$gene_id)

# rename columns from R4NCSL_x_count -> R4NCSL_x
colnames(count_mat) <- count_run_ids

# remove duplicated genes if any
count_mat <- count_mat[!duplicated(rownames(count_mat)), , drop = FALSE]

# DESeq2 requires integer counts
mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

############################################################
### 4. Optional: build CPM matrix too
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
### 5. Match run IDs to real sample names
############################################################

# Check all count columns have mapping
missing_map <- setdiff(colnames(count_mat), sample_map$run_id)
if (length(missing_map) > 0) {
  stop("These run IDs are missing from sample_map: ", paste(missing_map, collapse = ", "))
}

# Reorder sample_map to match matrix columns using base R
idx <- match(colnames(count_mat), sample_map$run_id)

if (any(is.na(idx))) {
  stop("Failed to match these count_mat columns: ",
       paste(colnames(count_mat)[is.na(idx)], collapse = ", "))
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

############################################################
### 6. Parse metadata from sample names
############################################################

# Example names:
# HMC3_Lipo48h_R1
# HMC3_siSTAT1_3_48h_R1
# HMC3_Lipo48h_LPS24h_R1
# HMC3_siSTAT1_3_48h_LPS24h_R3

parse_sample_metadata <- function(sample_names) {
  tibble(Sample_name = sample_names) %>%
    mutate(
      Knockdown = ifelse(grepl("siSTAT1", Sample_name), "siSTAT1", "Control"),
      Lipid = ifelse(grepl("Lipo48h", Sample_name), "Lipo48h", "Other"),
      LPS = ifelse(grepl("LPS24h", Sample_name), "LPS", "No_LPS"),
      Replicate = sub(".*_R([0-9]+)$", "R\\1", Sample_name)
    ) %>%
    mutate(
      Condition = case_when(
        Knockdown == "Control" & LPS == "No_LPS" ~ "Lipo_only",
        Knockdown == "siSTAT1" & LPS == "No_LPS" ~ "siSTAT1_Lipo",
        Knockdown == "Control" & LPS == "LPS" ~ "Lipo_LPS",
        Knockdown == "siSTAT1" & LPS == "LPS" ~ "siSTAT1_Lipo_LPS",
        TRUE ~ "Other"
      )
    ) %>%
    mutate(
      Knockdown = factor(Knockdown, levels = c("Control", "siSTAT1")),
      LPS = factor(LPS, levels = c("No_LPS", "LPS")),
      Condition = factor(
        Condition,
        levels = c("Lipo_only", "siSTAT1_Lipo", "Lipo_LPS", "siSTAT1_Lipo_LPS")
      )
    )
}

meta <- parse_sample_metadata(colnames(count_mat))
rownames(meta) <- meta$Sample_name

############################################################
### 7. Save parsed metadata
############################################################

write.xlsx(meta, file = file.path(outdir, "Parsed_sample_metadata.xlsx"), overwrite = TRUE)

############################################################
### 8. Build DESeq2 object - factorial model
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ Knockdown + LPS + Knockdown:LPS
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)

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
### 11. Group-based DE analysis
############################################################

meta$Group <- factor(
  meta$Condition,
  levels = c("Lipo_only", "siSTAT1_Lipo", "Lipo_LPS", "siSTAT1_Lipo_LPS")
)

dds_group <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ Group
)

keep2 <- rowSums(counts(dds_group) >= 10) >= 2
dds_group <- dds_group[keep2, ]
dds_group <- DESeq(dds_group)

vsd_group <- vst(dds_group, blind = FALSE)

make_res_table <- function(res_obj, annot_tbl) {
  as.data.frame(res_obj) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(annot_tbl, by = "ENSEMBL") %>%
    arrange(padj)
}

# Main group contrasts
res_siSTAT1_noLPS <- results(dds_group, contrast = c("Group", "siSTAT1_Lipo", "Lipo_only"))
res_LPS_control   <- results(dds_group, contrast = c("Group", "Lipo_LPS", "Lipo_only"))
res_siSTAT1_LPS   <- results(dds_group, contrast = c("Group", "siSTAT1_Lipo_LPS", "Lipo_LPS"))
res_LPS_siSTAT1   <- results(dds_group, contrast = c("Group", "siSTAT1_Lipo_LPS", "siSTAT1_Lipo"))

res_siSTAT1_noLPS_df <- make_res_table(res_siSTAT1_noLPS, gene_annot)
res_LPS_control_df   <- make_res_table(res_LPS_control, gene_annot)
res_siSTAT1_LPS_df   <- make_res_table(res_siSTAT1_LPS, gene_annot)
res_LPS_siSTAT1_df   <- make_res_table(res_LPS_siSTAT1, gene_annot)

write.xlsx(
  list(
    siSTAT1_noLPS = res_siSTAT1_noLPS_df,
    LPS_in_control = res_LPS_control_df,
    siSTAT1_in_LPS = res_siSTAT1_LPS_df,
    LPS_in_siSTAT1 = res_LPS_siSTAT1_df
  ),
  file = file.path(outdir, "DE_results_group_comparisons.xlsx"),
  overwrite = TRUE
)

############################################################
### 12. Optional factorial-model results
############################################################

# Inspect coefficient names first:
coef_names <- resultsNames(dds)
writeLines(coef_names)

# Usually these will correspond to:
# Knockdown_siSTAT1_vs_Control
# LPS_LPS_vs_No_LPS
# KnockdownsiSTAT1.LPSLPS

res_factorial_list <- list()

if ("Knockdown_siSTAT1_vs_Control" %in% coef_names) {
  res_factorial_list[["siSTAT1_main_effect"]] <-
    make_res_table(results(dds, name = "Knockdown_siSTAT1_vs_Control"), gene_annot)
}

if ("LPS_LPS_vs_No_LPS" %in% coef_names) {
  res_factorial_list[["LPS_main_effect"]] <-
    make_res_table(results(dds, name = "LPS_LPS_vs_No_LPS"), gene_annot)
}

# interaction term name can vary slightly, so detect it
interaction_term <- grep("Knockdown.*LPS|LPS.*Knockdown", coef_names, value = TRUE)

if (length(interaction_term) >= 1) {
  res_factorial_list[["Interaction_effect"]] <-
    make_res_table(results(dds, name = interaction_term[1]), gene_annot)
}

if (length(res_factorial_list) > 0) {
  write.xlsx(
    res_factorial_list,
    file = file.path(outdir, "DE_results_factorial_model.xlsx"),
    overwrite = TRUE
  )
}

############################################################
### 13. PCA plot
############################################################

pcaData <- plotPCA(vsd_group, intgroup = c("Group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - HMC3 siSTAT1 / LPS design") +
  theme_bw()

ggsave(
  file.path(outdir, "PCA_HMC3_siSTAT1_LPS.png"),
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
message("Genes retained in factorial model: ", nrow(dds))
message("Genes retained in group model: ", nrow(dds_group))

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
### 15. Variance stabilized matrix and normalized counts
############################################################

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

sample_dists <- dist(t(assay(vsd_group)))
sample_dist_mat <- as.matrix(sample_dists)

annotation_df <- as.data.frame(colData(dds_group)[, c("Group"), drop = FALSE])
rownames(annotation_df) <- colnames(dds_group)

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

p_siSTAT1_noLPS <- plot_volcano(res_siSTAT1_noLPS_df, "siSTAT1 vs Control (no LPS)")
p_LPS_control   <- plot_volcano(res_LPS_control_df, "LPS vs no LPS (Control)")
p_siSTAT1_LPS   <- plot_volcano(res_siSTAT1_LPS_df, "siSTAT1 vs Control (with LPS)")
p_LPS_siSTAT1   <- plot_volcano(res_LPS_siSTAT1_df, "LPS vs no LPS (siSTAT1)")

ggsave(file.path(outdir, "Volcano_siSTAT1_noLPS.png"), p_siSTAT1_noLPS, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_LPS_control.png"), p_LPS_control, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_siSTAT1_LPS.png"), p_siSTAT1_LPS, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_LPS_siSTAT1.png"), p_LPS_siSTAT1, width = 7, height = 5, dpi = 300)

ggsave(file.path(outdir, "Volcano_siSTAT1_noLPS.eps"), p_siSTAT1_noLPS, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_LPS_control.eps"), p_LPS_control, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_siSTAT1_LPS.eps"), p_siSTAT1_LPS, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_LPS_siSTAT1.eps"), p_LPS_siSTAT1, width = 7, height = 5, dpi = 300)

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

top_siSTAT1_noLPS <- get_top_genes(res_siSTAT1_noLPS_df)
top_LPS_control   <- get_top_genes(res_LPS_control_df)
top_siSTAT1_LPS   <- get_top_genes(res_siSTAT1_LPS_df)
top_LPS_siSTAT1   <- get_top_genes(res_LPS_siSTAT1_df)

write.xlsx(
  list(
    siSTAT1_noLPS_up = top_siSTAT1_noLPS$up,
    siSTAT1_noLPS_down = top_siSTAT1_noLPS$down,
    LPS_control_up = top_LPS_control$up,
    LPS_control_down = top_LPS_control$down,
    siSTAT1_LPS_up = top_siSTAT1_LPS$up,
    siSTAT1_LPS_down = top_siSTAT1_LPS$down,
    LPS_siSTAT1_up = top_LPS_siSTAT1$up,
    LPS_siSTAT1_down = top_LPS_siSTAT1$down
  ),
  file = file.path(outdir, "Top_DEGs.xlsx"),
  overwrite = TRUE
)

############################################################
### 19. Heatmap of selected genes
############################################################

selected_genes <- c("STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3", "MX1", "MX2", "NFKBIA", "CXCL8", "CCL2")

sel_annot <- gene_annot %>%
  filter(SYMBOL %in% selected_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(sel_annot) > 0) {
  sel_mat <- assay(vsd_group)[sel_annot$ENSEMBL, , drop = FALSE]
  rownames(sel_mat) <- sel_annot$SYMBOL[match(rownames(sel_mat), sel_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Group) %>%
    pull(Sample_name)
  
  sel_mat <- sel_mat[, ordered_samples, drop = FALSE]
  
  annotation_df2 <- meta[, c("Sample_name", "Group"), drop = FALSE]
  annotation_df2 <- as.data.frame(annotation_df2)
  rownames(annotation_df2) <- annotation_df2$Sample_name
  annotation_df2$Sample_name <- NULL
  annotation_df2 <- annotation_df2[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Selected_genes_heatmap.png"), width = 1200, height = 1000, res = 150)
  pheatmap(
    sel_mat,
    scale = "row",
    annotation_col = annotation_df2,
    cluster_cols = FALSE,
    main = "Selected genes heatmap"
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

ego_siSTAT1_noLPS_up <- run_go(res_siSTAT1_noLPS_df, "up")
ego_LPS_control_up   <- run_go(res_LPS_control_df, "up")
ego_siSTAT1_LPS_up   <- run_go(res_siSTAT1_LPS_df, "up")
ego_LPS_siSTAT1_up   <- run_go(res_LPS_siSTAT1_df, "up")

go_list <- list()
if (!is.null(ego_siSTAT1_noLPS_up)) go_list[["GO_siSTAT1_noLPS_up"]] <- as.data.frame(ego_siSTAT1_noLPS_up)
if (!is.null(ego_LPS_control_up))   go_list[["GO_LPS_control_up"]]   <- as.data.frame(ego_LPS_control_up)
if (!is.null(ego_siSTAT1_LPS_up))   go_list[["GO_siSTAT1_LPS_up"]]   <- as.data.frame(ego_siSTAT1_LPS_up)
if (!is.null(ego_LPS_siSTAT1_up))   go_list[["GO_LPS_siSTAT1_up"]]   <- as.data.frame(ego_LPS_siSTAT1_up)

if (length(go_list) > 0) {
  write.xlsx(go_list, file = file.path(outdir, "GO_enrichment_results.xlsx"), overwrite = TRUE)
}

save_go_plot <- function(go_obj, filename, title_text) {
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj)) > 0) {
    p <- dotplot(go_obj, showCategory = 10) + ggtitle(title_text)
    ggsave(file.path(outdir, filename), p, width = 8, height = 6, dpi = 300)
  }
}

save_go_plot(ego_siSTAT1_noLPS_up, "GO_siSTAT1_noLPS_up.png", "GO BP: siSTAT1 vs Control (no LPS)")
save_go_plot(ego_LPS_control_up,   "GO_LPS_control_up.png",   "GO BP: LPS in Control")
save_go_plot(ego_siSTAT1_LPS_up,   "GO_siSTAT1_LPS_up.png",   "GO BP: siSTAT1 vs Control (with LPS)")
save_go_plot(ego_LPS_siSTAT1_up,   "GO_LPS_siSTAT1_up.png",   "GO BP: LPS in siSTAT1")

save_go_plot(ego_siSTAT1_noLPS_up, "GO_siSTAT1_noLPS_up.eps", "GO BP: siSTAT1 vs Control (no LPS)")
save_go_plot(ego_LPS_control_up,   "GO_LPS_control_up.eps",   "GO BP: LPS in Control")
save_go_plot(ego_siSTAT1_LPS_up,   "GO_siSTAT1_LPS_up.eps",   "GO BP: siSTAT1 vs Control (with LPS)")
save_go_plot(ego_LPS_siSTAT1_up,   "GO_LPS_siSTAT1_up.eps",   "GO BP: LPS in siSTAT1")

############################################################
### 21. Venn diagram of significant DEGs
############################################################

sig_siSTAT1_noLPS <- res_siSTAT1_noLPS_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_LPS_control <- res_LPS_control_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_siSTAT1_LPS <- res_siSTAT1_LPS_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_LPS_siSTAT1 <- res_LPS_siSTAT1_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

# 4-way Venn can get crowded, but this will still run
gene_sets <- list(
  siSTAT1_noLPS = sig_siSTAT1_noLPS,
  LPS_control   = sig_LPS_control,
  siSTAT1_LPS   = sig_siSTAT1_LPS,
  LPS_siSTAT1   = sig_LPS_siSTAT1
)

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  fill = c("purple", "green3", "red", "orange"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.1,
  main = "DEG overlap across contrasts"
)

png(file.path(outdir, "DEG_venn.png"), width = 1400, height = 1200, res = 150)
grid.draw(venn.plot)
dev.off()

common_genes_all <- Reduce(intersect, gene_sets)

write.xlsx(
  list(
    common_genes_all = data.frame(Gene = common_genes_all),
    siSTAT1_noLPS = data.frame(Gene = sig_siSTAT1_noLPS),
    LPS_control = data.frame(Gene = sig_LPS_control),
    siSTAT1_LPS = data.frame(Gene = sig_siSTAT1_LPS),
    LPS_siSTAT1 = data.frame(Gene = sig_LPS_siSTAT1)
  ),
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

lipid_res_siSTAT1_noLPS <- extract_module_results(res_siSTAT1_noLPS_df, lipid_modules, "siSTAT1_noLPS")
lipid_res_LPS_control   <- extract_module_results(res_LPS_control_df, lipid_modules, "LPS_control")
lipid_res_siSTAT1_LPS   <- extract_module_results(res_siSTAT1_LPS_df, lipid_modules, "siSTAT1_LPS")
lipid_res_LPS_siSTAT1   <- extract_module_results(res_LPS_siSTAT1_df, lipid_modules, "LPS_siSTAT1")

write.xlsx(
  list(
    siSTAT1_noLPS = lipid_res_siSTAT1_noLPS,
    LPS_control = lipid_res_LPS_control,
    siSTAT1_LPS = lipid_res_siSTAT1_LPS,
    LPS_siSTAT1 = lipid_res_LPS_siSTAT1
  ),
  file = file.path(outdir, "Focused_lipid_gene_results.xlsx"),
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
  keep_genes <- intersect(lipid_annot$ENSEMBL, rownames(assay(vsd_group)))
  lipid_mat <- assay(vsd_group)[keep_genes, , drop = FALSE]
  rownames(lipid_mat) <- lipid_annot$SYMBOL[match(rownames(lipid_mat), lipid_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Group) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(lipid_mat))
  lipid_mat <- lipid_mat[, ordered_samples, drop = FALSE]
  
  anno_lipid <- meta[, c("Sample_name", "Group"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_lipid) <- anno_lipid$Sample_name
  anno_lipid$Sample_name <- NULL
  anno_lipid <- anno_lipid[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Focused_lipid_heatmap.png"), width = 1600, height = 1300, res = 160)
  pheatmap(
    lipid_mat,
    scale = "row",
    annotation_col = anno_lipid,
    cluster_cols = FALSE,
    main = "Focused lipid / LD / lysosome / stress / IFN genes"
  )
  dev.off()
}

#eps
all_lipid_genes <- unique(unlist(lipid_modules))

lipid_annot <- gene_annot %>%
  filter(SYMBOL %in% all_lipid_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(lipid_annot) > 0) {
  keep_genes <- intersect(lipid_annot$ENSEMBL, rownames(assay(vsd_group)))
  lipid_mat <- assay(vsd_group)[keep_genes, , drop = FALSE]
  rownames(lipid_mat) <- lipid_annot$SYMBOL[match(rownames(lipid_mat), lipid_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Group) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(lipid_mat))
  lipid_mat <- lipid_mat[, ordered_samples, drop = FALSE]
  
  anno_lipid <- meta[, c("Sample_name", "Group"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_lipid) <- anno_lipid$Sample_name
  anno_lipid$Sample_name <- NULL
  anno_lipid <- anno_lipid[ordered_samples, , drop = FALSE]
  
  # PNG export
  png(file.path(outdir, "Focused_lipid_heatmap.png"), width = 1600, height = 1300, res = 160)
  pheatmap(
    lipid_mat,
    scale = "row",
    annotation_col = anno_lipid,
    cluster_cols = FALSE,
    main = "Focused lipid / LD / lysosome / stress / IFN genes"
  )
  dev.off()
  
  # EPS export
  cairo_ps(
    file = file.path(outdir, "Focused_lipid_heatmap.eps"),
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

module_summary_siSTAT1_noLPS <- summarize_modules(res_siSTAT1_noLPS_df, lipid_modules, "siSTAT1_noLPS")
module_summary_LPS_control   <- summarize_modules(res_LPS_control_df, lipid_modules, "LPS_control")
module_summary_siSTAT1_LPS   <- summarize_modules(res_siSTAT1_LPS_df, lipid_modules, "siSTAT1_LPS")
module_summary_LPS_siSTAT1   <- summarize_modules(res_LPS_siSTAT1_df, lipid_modules, "LPS_siSTAT1")

module_summary_all <- bind_rows(
  module_summary_siSTAT1_noLPS,
  module_summary_LPS_control,
  module_summary_siSTAT1_LPS,
  module_summary_LPS_siSTAT1
)

write.xlsx(
  module_summary_all,
  file = file.path(outdir, "Lipid_module_summary.xlsx"),
  overwrite = TRUE
)

p_module_mean <- ggplot(module_summary_all, aes(x = Module, y = mean_log2FC, fill = Contrast)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Lipid / LD module summary", x = "", y = "Mean log2FC")

ggsave(file.path(outdir, "Lipid_module_summary_barplot.png"),
       p_module_mean, width = 10, height = 6.5, dpi = 300)

ggsave(
  file.path(outdir, "Lipid_module_summary_barplot.eps"),
  p_module_mean,
  width = 10,
  height = 6.5,
  device = cairo_ps
)

p_module_sig <- ggplot(module_summary_all, aes(x = Module, y = n_up_sig, fill = Contrast)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Significant upregulated genes per module", x = "", y = "Count")

ggsave(file.path(outdir, "Lipid_module_sig_up_barplot.png"),
       p_module_sig, width = 10, height = 6.5, dpi = 300)

ggsave(
  file.path(outdir, "Lipid_module_sig_up_barplot.eps"),
  p_module_mean,
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
  "DDIT3", "LAMP1", "STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3",
  "MX1", "MX2", "NFKBIA", "CXCL8"
)

extract_qpcr_panel <- function(res_df, genes, contrast_label) {
  res_df %>%
    filter(SYMBOL %in% genes) %>%
    select(ENSEMBL, SYMBOL, GENENAME, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(match(SYMBOL, genes)) %>%
    mutate(Contrast = contrast_label)
}

qpcr_panel_siSTAT1_noLPS <- extract_qpcr_panel(res_siSTAT1_noLPS_df, custom_qpcr_panel, "siSTAT1_noLPS")
qpcr_panel_LPS_control   <- extract_qpcr_panel(res_LPS_control_df, custom_qpcr_panel, "LPS_control")
qpcr_panel_siSTAT1_LPS   <- extract_qpcr_panel(res_siSTAT1_LPS_df, custom_qpcr_panel, "siSTAT1_LPS")
qpcr_panel_LPS_siSTAT1   <- extract_qpcr_panel(res_LPS_siSTAT1_df, custom_qpcr_panel, "LPS_siSTAT1")

write.xlsx(
  list(
    siSTAT1_noLPS = qpcr_panel_siSTAT1_noLPS,
    LPS_control = qpcr_panel_LPS_control,
    siSTAT1_LPS = qpcr_panel_siSTAT1_LPS,
    LPS_siSTAT1 = qpcr_panel_LPS_siSTAT1
  ),
  file = file.path(outdir, "Custom_lipid_qPCR_panel.xlsx"),
  overwrite = TRUE
)

############################################################
### 26. Targeted panel heatmap
############################################################

panel_annot <- gene_annot %>%
  filter(SYMBOL %in% custom_qpcr_panel) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(panel_annot) > 0) {
  keep_genes <- intersect(panel_annot$ENSEMBL, rownames(assay(vsd_group)))
  panel_mat <- assay(vsd_group)[keep_genes, , drop = FALSE]
  rownames(panel_mat) <- panel_annot$SYMBOL[match(rownames(panel_mat), panel_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    arrange(Group) %>%
    pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(panel_mat))
  panel_mat <- panel_mat[, ordered_samples, drop = FALSE]
  
  anno_panel <- meta[, c("Sample_name", "Group"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_panel) <- anno_panel$Sample_name
  anno_panel$Sample_name <- NULL
  anno_panel <- anno_panel[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Custom_lipid_qPCR_panel_heatmap.png"), width = 1600, height = 1100, res = 160)
  pheatmap(
    panel_mat,
    scale = "row",
    annotation_col = anno_panel,
    cluster_cols = FALSE,
    main = "Custom lipid / STAT1 axis panel"
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

qpcr_panel_siSTAT1_noLPS_interpret <- interpret_panel_direction(qpcr_panel_siSTAT1_noLPS)
qpcr_panel_LPS_control_interpret   <- interpret_panel_direction(qpcr_panel_LPS_control)
qpcr_panel_siSTAT1_LPS_interpret   <- interpret_panel_direction(qpcr_panel_siSTAT1_LPS)
qpcr_panel_LPS_siSTAT1_interpret   <- interpret_panel_direction(qpcr_panel_LPS_siSTAT1)

write.xlsx(
  list(
    siSTAT1_noLPS = qpcr_panel_siSTAT1_noLPS_interpret,
    LPS_control = qpcr_panel_LPS_control_interpret,
    siSTAT1_LPS = qpcr_panel_siSTAT1_LPS_interpret,
    LPS_siSTAT1 = qpcr_panel_LPS_siSTAT1_interpret
  ),
  file = file.path(outdir, "Custom_lipid_qPCR_panel_interpretation.xlsx"),
  overwrite = TRUE
)

############################################################
### Selected Reactome GSEA enrichment plots
############################################################

############################################################
### Reactome GSEA setup for your 4 contrasts
############################################################

library(ReactomePA)
library(enrichplot)

make_ranked_vector <- function(res_df, use_stat = TRUE) {
  if (use_stat && "stat" %in% colnames(res_df)) {
    ranked <- res_df %>%
      filter(!is.na(ENTREZID), !is.na(stat)) %>%
      group_by(ENTREZID) %>%
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(stat)) %>%
      select(ENTREZID, stat)
    
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
      select(ENTREZID, rank_score)
    
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

geneList_siSTAT1_noLPS <- make_ranked_vector(res_siSTAT1_noLPS_df, use_stat = TRUE)
geneList_LPS_control   <- make_ranked_vector(res_LPS_control_df, use_stat = TRUE)
geneList_siSTAT1_LPS   <- make_ranked_vector(res_siSTAT1_LPS_df, use_stat = TRUE)
geneList_LPS_siSTAT1   <- make_ranked_vector(res_LPS_siSTAT1_df, use_stat = TRUE)

gsea_react_siSTAT1_noLPS <- run_reactome_gsea(geneList_siSTAT1_noLPS)
gsea_react_LPS_control   <- run_reactome_gsea(geneList_LPS_control)
gsea_react_siSTAT1_LPS   <- run_reactome_gsea(geneList_siSTAT1_LPS)
gsea_react_LPS_siSTAT1   <- run_reactome_gsea(geneList_LPS_siSTAT1)

class(gsea_react_siSTAT1_noLPS)
class(gsea_react_LPS_control)
class(gsea_react_siSTAT1_LPS)
class(gsea_react_LPS_siSTAT1)
if (!is.null(gsea_react_siSTAT1_noLPS)) nrow(as.data.frame(gsea_react_siSTAT1_noLPS))
if (!is.null(gsea_react_LPS_control)) nrow(as.data.frame(gsea_react_LPS_control))
if (!is.null(gsea_react_siSTAT1_LPS)) nrow(as.data.frame(gsea_react_siSTAT1_LPS))
if (!is.null(gsea_react_LPS_siSTAT1)) nrow(as.data.frame(gsea_react_LPS_siSTAT1))

############################################################
### Selected Reactome GSEA enrichment plots
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
  "cholesterol",
  "fatty acid",
  "lysosome",
  "autophagy",
  "ER"
)

plot_selected_gseapaths(
  gsea_obj = gsea_react_siSTAT1_noLPS,
  pattern_vec = pattern_vec,
  prefix = "Reactome_GSEA_siSTAT1_noLPS_selected",
  outdir = outdir
)

plot_selected_gseapaths(
  gsea_obj = gsea_react_LPS_control,
  pattern_vec = pattern_vec,
  prefix = "Reactome_GSEA_LPS_control_selected",
  outdir = outdir
)

plot_selected_gseapaths(
  gsea_obj = gsea_react_siSTAT1_LPS,
  pattern_vec = pattern_vec,
  prefix = "Reactome_GSEA_siSTAT1_LPS_selected",
  outdir = outdir
)

plot_selected_gseapaths(
  gsea_obj = gsea_react_LPS_siSTAT1,
  pattern_vec = pattern_vec,
  prefix = "Reactome_GSEA_LPS_siSTAT1_selected",
  outdir = outdir
)

############################################################
### 28. Final add-on summary
############################################################

message("Extended downstream analysis completed successfully.")



