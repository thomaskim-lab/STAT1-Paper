############################################
### BULK RNA-SEQ ANALYSIS PIPELINE
### Thomas dataset
############################################

###############################
### 0. Install packages if needed
###############################
# install.packages(c("tidyverse", "pheatmap", "ggrepel", "openxlsx", "VennDiagram"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "apeglm"))

###############################
### 1. Load libraries
###############################
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

###############################
### 2. Set input/output paths
###############################
setwd("C:/Users/au731006/Downloads")
counts_file <- "C:/Users/au731006/Downloads/RNASeq/Thomas.rawcounts.csv"
meta_file   <- "C:/Users/au731006/Downloads/RNASeq/Thomas_metadata.csv"
outdir      <- "C:/Users/au731006/Downloads/Plasmidsaurus/Output/LPS"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###############################
### 3. Read counts and metadata
###############################
counts_df <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)

# Standardize first column name in counts
colnames(counts_df)[1] <- "gene_id"

# Standardize metadata sample column name
if ("Sample name" %in% colnames(meta)) {
  colnames(meta)[colnames(meta) == "Sample name"] <- "Sample_name"
}

# Check sample names
sample_cols <- colnames(counts_df)[-1]

# Keep only metadata rows corresponding to count columns
meta <- meta %>% filter(Sample_name %in% sample_cols)

# Reorder metadata to match counts matrix column order
meta <- meta %>% slice(match(sample_cols, Sample_name))

if (!all(meta$Sample_name == sample_cols)) {
  stop("Metadata sample order does not match counts matrix columns.")
}

###############################
### 4. Prepare count matrix
###############################
count_mat <- as.matrix(counts_df[, -1])
rownames(count_mat) <- counts_df$gene_id

# Remove Ensembl version if present
rownames(count_mat) <- sub("\\..*$", "", rownames(count_mat))

# Remove duplicated Ensembl IDs
count_mat <- count_mat[!duplicated(rownames(count_mat)), ]

# Force numeric/integer
mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

###############################
### 5. Prepare metadata
###############################
meta <- meta %>%
  mutate(
    Time = factor(Time, levels = c("24h", "48h", "72h")),
    Treatment = factor(Treatment, levels = c("Control", "LPS")),
    Time_Treatment = factor(Time_Treatment)
  )

rownames(meta) <- meta$Sample_name

###############################
### 6. Create DESeq2 object
###############################
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ Time + Treatment + Time:Treatment
)

# Filter low-count genes before DESeq
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

###############################
### 7. Annotation
###############################
gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(dds),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

###############################
### 8. Variance stabilized data
###############################
vsd <- vst(dds, blind = FALSE)

###############################
### 9. PCA plot
###############################
pcaData <- plotPCA(vsd, intgroup = c("Time", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = Time, shape = Treatment, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of bulk RNA-seq samples") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "PCA_bulkRNAseq.png"), p_pca, width = 7, height = 5, dpi = 300)

###############################
### 10. Sample distance heatmap
###############################
sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

annotation_df <- as.data.frame(colData(dds)[, c("Time", "Treatment")])

png(file.path(outdir, "Sample_distance_heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sample_dist_mat,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "Sample-to-sample distances"
)
dev.off()

###############################
### 11. Save normalized matrices
###############################
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

###############################
### 12. Per-timepoint DE function
###############################
run_timepoint_deseq <- function(dds_obj, timepoint) {
  dds_sub <- dds_obj[, dds_obj$Time == timepoint]
  dds_sub$Time <- droplevels(dds_sub$Time)
  dds_sub$Treatment <- droplevels(dds_sub$Treatment)
  design(dds_sub) <- ~ Treatment
  dds_sub <- DESeq(dds_sub)
  
  res <- results(dds_sub, contrast = c("Treatment", "LPS", "Control"))
  
  # shrink LFC
  res_shrunk <- lfcShrink(
    dds_sub,
    contrast = c("Treatment", "LPS", "Control"),
    res = res,
    type = "normal"
  )
  
  return(list(dds = dds_sub, res = res_shrunk))
}

out_24 <- run_timepoint_deseq(dds, "24h")
out_48 <- run_timepoint_deseq(dds, "48h")
out_72 <- run_timepoint_deseq(dds, "72h")

res_24 <- out_24$res
res_48 <- out_48$res
res_72 <- out_72$res

###############################
### 13. Convert DESeq results to tables
###############################
make_res_table <- function(res_obj, annot_tbl) {
  as.data.frame(res_obj) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(annot_tbl, by = "ENSEMBL") %>%
    arrange(padj)
}

res_24_df <- make_res_table(res_24, gene_annot)
res_48_df <- make_res_table(res_48, gene_annot)
res_72_df <- make_res_table(res_72, gene_annot)

# Save full DE tables
write.xlsx(
  list(
    DE_24h = res_24_df,
    DE_48h = res_48_df,
    DE_72h = res_72_df
  ),
  file = file.path(outdir, "DE_results_by_timepoint.xlsx"),
  overwrite = TRUE
)

###############################
### 14. Volcano plot function
###############################
plot_volcano <- function(res_df, title_text, label_top_n = 15) {
  df <- res_df %>%
    mutate(
      sig = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 0.5  ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -0.5 ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(pvalue)
    )
  
  top_labels <- df %>%
    filter(sig != "NS") %>%
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

p24 <- plot_volcano(res_24_df, "24h: LPS vs Control")
p48 <- plot_volcano(res_48_df, "48h: LPS vs Control")
p72 <- plot_volcano(res_72_df, "72h: LPS vs Control")

ggsave(file.path(outdir, "Volcano_24h.png"), p24, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_48h.png"), p48, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "Volcano_72h.png"), p72, width = 7, height = 5, dpi = 300)

###############################
### 15. Top DEGs
###############################
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

top24 <- get_top_genes(res_24_df)
top48 <- get_top_genes(res_48_df)
top72 <- get_top_genes(res_72_df)

write.xlsx(
  list(
    top24_up = top24$up,
    top24_down = top24$down,
    top48_up = top48$up,
    top48_down = top48$down,
    top72_up = top72$up,
    top72_down = top72$down
  ),
  file = file.path(outdir, "Top_DEGs.xlsx"),
  overwrite = TRUE
)

###############################
### 16. Heatmap of selected genes
###############################
selected_genes <- c("MX1", "MX2", "ISG20", "CXCL8", "IFIT1", "IFIT3", "IFI16", "NFKBIA")

sel_annot <- gene_annot %>%
  filter(SYMBOL %in% selected_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

sel_mat <- assay(vsd)[sel_annot$ENSEMBL, , drop = FALSE]
rownames(sel_mat) <- sel_annot$SYMBOL[match(rownames(sel_mat), sel_annot$ENSEMBL)]

# reorder columns by Time then Treatment
ordered_samples <- meta %>%
  arrange(Time, Treatment) %>%
  pull(Sample_name)

sel_mat <- sel_mat[, ordered_samples, drop = FALSE]
annotation_df2 <- meta[ordered_samples, c("Time", "Treatment"), drop = FALSE]

png(file.path(outdir, "Selected_genes_heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sel_mat,
  scale = "row",
  annotation_col = annotation_df2,
  cluster_cols = FALSE,
  main = "Selected genes heatmap"
)
dev.off()

###############################
### 17. GO enrichment
###############################
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
  
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  return(ego)
}

ego_24_up <- run_go(res_24_df, "up")
ego_48_up <- run_go(res_48_df, "up")
ego_72_up <- run_go(res_72_df, "up")

if (!is.null(ego_24_up)) {
  p <- dotplot(ego_24_up, showCategory = 10) + ggtitle("GO BP enrichment: 24h upregulated")
  ggsave(file.path(outdir, "GO_24h_up.png"), p, width = 8, height = 6, dpi = 300)
}

if (!is.null(ego_48_up)) {
  p <- dotplot(ego_48_up, showCategory = 10) + ggtitle("GO BP enrichment: 48h upregulated")
  ggsave(file.path(outdir, "GO_48h_up.png"), p, width = 8, height = 6, dpi = 300)
}

if (!is.null(ego_72_up)) {
  p <- dotplot(ego_72_up, showCategory = 10) + ggtitle("GO BP enrichment: 72h upregulated")
  ggsave(file.path(outdir, "GO_72h_up.png"), p, width = 8, height = 6, dpi = 300)
}

# Save GO tables
go_list <- list()

if (!is.null(ego_24_up)) go_list[["GO_24h_up"]] <- as.data.frame(ego_24_up)
if (!is.null(ego_48_up)) go_list[["GO_48h_up"]] <- as.data.frame(ego_48_up)
if (!is.null(ego_72_up)) go_list[["GO_72h_up"]] <- as.data.frame(ego_72_up)

if (length(go_list) > 0) {
  write.xlsx(go_list, file = file.path(outdir, "GO_enrichment_results.xlsx"), overwrite = TRUE)
}

###############################
### 18. Venn diagram of significant DEGs
###############################
sig_24 <- res_24_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_48 <- res_48_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_72 <- res_72_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

gene_sets <- list(
  "24h" = sig_24,
  "48h" = sig_48,
  "72h" = sig_72
)

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  fill = c("purple", "green3", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20, 180),
  main = "DEG overlap across time points"
)

png(file.path(outdir, "DEG_venn.png"), width = 1200, height = 1000, res = 150)
grid.draw(venn.plot)
dev.off()

common_genes <- Reduce(intersect, gene_sets)

unique_24 <- setdiff(sig_24, union(sig_48, sig_72))
unique_48 <- setdiff(sig_48, union(sig_24, sig_72))
unique_72 <- setdiff(sig_72, union(sig_24, sig_48))

write.xlsx(
  list(
    common_genes = data.frame(Gene = common_genes),
    unique_24h = data.frame(Gene = unique_24),
    unique_48h = data.frame(Gene = unique_48),
    unique_72h = data.frame(Gene = unique_72)
  ),
  file = file.path(outdir, "DEG_overlap_lists.xlsx"),
  overwrite = TRUE
)

###############################
### 19. Optional: STAT1-like exploratory analysis
###############################
find_stat1_like <- function(dds_sub, res_df, gene_symbol = "STAT1", expr_similarity = 0.4, fc_similarity = 0.1) {
  norm_counts <- counts(dds_sub, normalized = TRUE)
  
  annot_sub <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = rownames(dds_sub),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  expr_df <- as.data.frame(norm_counts) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(annot_sub, by = "ENSEMBL")
  
  expr_df <- expr_df %>%
    mutate(
      mean_Control = rowMeans(select(., colnames(norm_counts)[dds_sub$Treatment == "Control"]), na.rm = TRUE),
      mean_LPS     = rowMeans(select(., colnames(norm_counts)[dds_sub$Treatment == "LPS"]), na.rm = TRUE)
    ) %>%
    select(ENSEMBL, SYMBOL, mean_Control, mean_LPS)
  
  merged <- res_df %>%
    left_join(expr_df, by = c("ENSEMBL", "SYMBOL"))
  
  stat1_row <- merged %>% filter(SYMBOL == gene_symbol) %>% slice(1)
  
  if (nrow(stat1_row) == 0) return(NULL)
  
  fc_low <- stat1_row$log2FoldChange * (1 - fc_similarity)
  fc_high <- stat1_row$log2FoldChange * (1 + fc_similarity)
  
  c_low <- stat1_row$mean_Control * (1 - expr_similarity)
  c_high <- stat1_row$mean_Control * (1 + expr_similarity)
  l_low <- stat1_row$mean_LPS * (1 - expr_similarity)
  l_high <- stat1_row$mean_LPS * (1 + expr_similarity)
  
  out <- merged %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(mean_Control),
      !is.na(mean_LPS),
      log2FoldChange >= fc_low,
      log2FoldChange <= fc_high,
      mean_Control >= c_low,
      mean_Control <= c_high,
      mean_LPS >= l_low,
      mean_LPS <= l_high
    ) %>%
    arrange(padj)
  
  return(out)
}

stat1_like_24 <- find_stat1_like(out_24$dds, res_24_df)
stat1_like_48 <- find_stat1_like(out_48$dds, res_48_df)
stat1_like_72 <- find_stat1_like(out_72$dds, res_72_df)

stat1_list <- list()
if (!is.null(stat1_like_24)) stat1_list[["STAT1_like_24h"]] <- stat1_like_24
if (!is.null(stat1_like_48)) stat1_list[["STAT1_like_48h"]] <- stat1_like_48
if (!is.null(stat1_like_72)) stat1_list[["STAT1_like_72h"]] <- stat1_like_72

if (length(stat1_list) > 0) {
  write.xlsx(stat1_list, file = file.path(outdir, "STAT1_like_genes.xlsx"), overwrite = TRUE)
}

###############################
### 20. Session summary
###############################
cat("Analysis complete.\n")
cat("Results saved to:\n", outdir, "\n")

############################################################
### 21. ADD-ON: Reactome / GSEA / lipid-focused analysis
############################################################

# Required packages
# install.packages(c("msigdbr", "enrichplot", "patchwork"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("ReactomePA", "fgsea"))

library(ReactomePA)
library(msigdbr)
library(fgsea)
library(enrichplot)
library(patchwork)

############################################################
### 21A. Reactome ORA from significant DEGs
############################################################

run_reactome_ora <- function(res_df, direction = c("up", "down"),
                             padj_cutoff = 0.05, lfc_cutoff = 0.5) {
  direction <- match.arg(direction)
  
  if (direction == "up") {
    genes <- res_df %>%
      dplyr::filter(!is.na(padj), padj < padj_cutoff,
                    log2FoldChange > lfc_cutoff,
                    !is.na(ENTREZID)) %>%
      dplyr::pull(ENTREZID) %>%
      unique()
  } else {
    genes <- res_df %>%
      dplyr::filter(!is.na(padj), padj < padj_cutoff,
                    log2FoldChange < -lfc_cutoff,
                    !is.na(ENTREZID)) %>%
      dplyr::pull(ENTREZID) %>%
      unique()
  }
  
  if (length(genes) < 5) return(NULL)
  
  ReactomePA::enrichPathway(
    gene = genes,
    organism = "human",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

react_ora_24_up   <- run_reactome_ora(res_24_df, "up")
react_ora_24_down <- run_reactome_ora(res_24_df, "down")
react_ora_48_up   <- run_reactome_ora(res_48_df, "up")
react_ora_48_down <- run_reactome_ora(res_48_df, "down")
react_ora_72_up   <- run_reactome_ora(res_72_df, "up")
react_ora_72_down <- run_reactome_ora(res_72_df, "down")

reactome_ora_list <- list()
if (!is.null(react_ora_24_up))   reactome_ora_list[["Reactome_ORA_24h_up"]]   <- as.data.frame(react_ora_24_up)
if (!is.null(react_ora_24_down)) reactome_ora_list[["Reactome_ORA_24h_down"]] <- as.data.frame(react_ora_24_down)
if (!is.null(react_ora_48_up))   reactome_ora_list[["Reactome_ORA_48h_up"]]   <- as.data.frame(react_ora_48_up)
if (!is.null(react_ora_48_down)) reactome_ora_list[["Reactome_ORA_48h_down"]] <- as.data.frame(react_ora_48_down)
if (!is.null(react_ora_72_up))   reactome_ora_list[["Reactome_ORA_72h_up"]]   <- as.data.frame(react_ora_72_up)
if (!is.null(react_ora_72_down)) reactome_ora_list[["Reactome_ORA_72h_down"]] <- as.data.frame(react_ora_72_down)

if (length(reactome_ora_list) > 0) {
  write.xlsx(
    reactome_ora_list,
    file = file.path(outdir, "Reactome_ORA_results.xlsx"),
    overwrite = TRUE
  )
}

save_dotplot_if_valid <- function(enrich_obj, filename, title_text, ncat = 15,
                                  width = 9, height = 6) {
  if (!is.null(enrich_obj)) {
    df <- as.data.frame(enrich_obj)
    if (nrow(df) > 0) {
      p <- enrichplot::dotplot(enrich_obj, showCategory = ncat) +
        ggtitle(title_text)
      ggsave(file.path(outdir, filename), p, width = width, height = height, dpi = 300)
    }
  }
}

save_dotplot_if_valid(react_ora_24_up,   "Reactome_ORA_24h_up.png",   "Reactome ORA - 24h Up")
save_dotplot_if_valid(react_ora_24_down, "Reactome_ORA_24h_down.png", "Reactome ORA - 24h Down")
save_dotplot_if_valid(react_ora_48_up,   "Reactome_ORA_48h_up.png",   "Reactome ORA - 48h Up")
save_dotplot_if_valid(react_ora_48_down, "Reactome_ORA_48h_down.png", "Reactome ORA - 48h Down")
save_dotplot_if_valid(react_ora_72_up,   "Reactome_ORA_72h_up.png",   "Reactome ORA - 72h Up")
save_dotplot_if_valid(react_ora_72_down, "Reactome_ORA_72h_down.png", "Reactome ORA - 72h Down")

save_dotplot_if_valid(react_ora_24_up,   "Reactome_ORA_24h_up.eps",   "Reactome ORA - 24h Up")
save_dotplot_if_valid(react_ora_24_down, "Reactome_ORA_24h_down.eps", "Reactome ORA - 24h Down")
save_dotplot_if_valid(react_ora_48_up,   "Reactome_ORA_48h_up.eps",   "Reactome ORA - 48h Up")
save_dotplot_if_valid(react_ora_48_down, "Reactome_ORA_48h_down.eps", "Reactome ORA - 48h Down")
save_dotplot_if_valid(react_ora_72_up,   "Reactome_ORA_72h_up.eps",   "Reactome ORA - 72h Up")
save_dotplot_if_valid(react_ora_72_down, "Reactome_ORA_72h_down.eps", "Reactome ORA - 72h Down")

############################################################
### 21B. Ranked GSEA helpers
############################################################

make_ranked_vector <- function(res_df, use_stat = TRUE) {
  if (use_stat && "stat" %in% colnames(res_df)) {
    ranked <- res_df %>%
      dplyr::filter(!is.na(ENTREZID), !is.na(stat)) %>%
      dplyr::group_by(ENTREZID) %>%
      dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(stat)) %>%
      dplyr::select(ENTREZID, stat)
    
    geneList <- ranked$stat
    names(geneList) <- ranked$ENTREZID
  } else {
    ranked <- res_df %>%
      dplyr::filter(!is.na(ENTREZID), !is.na(log2FoldChange), !is.na(pvalue), pvalue > 0) %>%
      dplyr::mutate(rank_score = sign(log2FoldChange) * -log10(pvalue)) %>%
      dplyr::group_by(ENTREZID) %>%
      dplyr::slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(rank_score)) %>%
      dplyr::select(ENTREZID, rank_score)
    
    geneList <- ranked$rank_score
    names(geneList) <- ranked$ENTREZID
  }
  
  geneList <- sort(geneList, decreasing = TRUE)
  return(geneList)
}

geneList_24 <- make_ranked_vector(res_24_df, use_stat = TRUE)
geneList_48 <- make_ranked_vector(res_48_df, use_stat = TRUE)
geneList_72 <- make_ranked_vector(res_72_df, use_stat = TRUE)

############################################################
### 21C. Reactome GSEA
############################################################

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

gsea_react_24 <- run_reactome_gsea(geneList_24)
gsea_react_48 <- run_reactome_gsea(geneList_48)
gsea_react_72 <- run_reactome_gsea(geneList_72)

reactome_gsea_list <- list()
if (!is.null(gsea_react_24) && nrow(as.data.frame(gsea_react_24)) > 0) reactome_gsea_list[["Reactome_GSEA_24h"]] <- as.data.frame(gsea_react_24)
if (!is.null(gsea_react_48) && nrow(as.data.frame(gsea_react_48)) > 0) reactome_gsea_list[["Reactome_GSEA_48h"]] <- as.data.frame(gsea_react_48)
if (!is.null(gsea_react_72) && nrow(as.data.frame(gsea_react_72)) > 0) reactome_gsea_list[["Reactome_GSEA_72h"]] <- as.data.frame(gsea_react_72)

if (length(reactome_gsea_list) > 0) {
  write.xlsx(
    reactome_gsea_list,
    file = file.path(outdir, "Reactome_GSEA_results.xlsx"),
    overwrite = TRUE
  )
}

save_dotplot_if_valid(gsea_react_24, "Reactome_GSEA_24h.png", "Reactome GSEA - 24h")
save_dotplot_if_valid(gsea_react_48, "Reactome_GSEA_48h.png", "Reactome GSEA - 48h")
save_dotplot_if_valid(gsea_react_72, "Reactome_GSEA_72h.png", "Reactome GSEA - 72h")

save_dotplot_if_valid(gsea_react_24, "Reactome_GSEA_24h.eps", "Reactome GSEA - 24h")
save_dotplot_if_valid(gsea_react_48, "Reactome_GSEA_48h.eps", "Reactome GSEA - 48h")
save_dotplot_if_valid(gsea_react_72, "Reactome_GSEA_72h.eps", "Reactome GSEA - 72h")

############################################################
### 21D. Hallmark GSEA using fgsea
############################################################

msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

hallmark_list <- split(msig_h$entrez_gene, msig_h$gs_name)

run_fgsea_simple <- function(geneList, pathways, minSize = 10, maxSize = 500, nperm = 10000) {
  fgsea::fgsea(
    pathways = pathways,
    stats = geneList,
    minSize = minSize,
    maxSize = maxSize,
    nperm = nperm
  ) %>%
    dplyr::arrange(padj)
}

fgsea_h_24 <- run_fgsea_simple(geneList_24, hallmark_list)
fgsea_h_48 <- run_fgsea_simple(geneList_48, hallmark_list)
fgsea_h_72 <- run_fgsea_simple(geneList_72, hallmark_list)

write.xlsx(
  list(
    Hallmark_GSEA_24h = fgsea_h_24,
    Hallmark_GSEA_48h = fgsea_h_48,
    Hallmark_GSEA_72h = fgsea_h_72
  ),
  file = file.path(outdir, "Hallmark_GSEA_results.xlsx"),
  overwrite = TRUE
)

plot_fgsea_top <- function(fgsea_res, title_text, filename, top_n = 15) {
  df <- fgsea_res %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::slice_min(order_by = padj, n = top_n) %>%
    dplyr::mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  if (nrow(df) == 0) return(NULL)
  
  p <- ggplot(df, aes(x = NES, y = pathway, size = size, color = -log10(padj))) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(title = title_text, y = "", x = "NES", color = "-log10(adj.P)", size = "Geneset size")
  
  ggsave(file.path(outdir, filename), p, width = 9, height = 6, dpi = 300)
  return(p)
}

plot_fgsea_top(fgsea_h_24, "Hallmark GSEA - 24h", "Hallmark_GSEA_24h.png")
plot_fgsea_top(fgsea_h_48, "Hallmark GSEA - 48h", "Hallmark_GSEA_48h.png")
plot_fgsea_top(fgsea_h_72, "Hallmark GSEA - 72h", "Hallmark_GSEA_72h.png")

plot_fgsea_top(fgsea_h_24, "Hallmark GSEA - 24h", "Hallmark_GSEA_24h.eps")
plot_fgsea_top(fgsea_h_48, "Hallmark GSEA - 48h", "Hallmark_GSEA_48h.eps")
plot_fgsea_top(fgsea_h_72, "Hallmark GSEA - 72h", "Hallmark_GSEA_72h.eps")

############################################################
### 21E. Custom lipid / LD / IFN modules
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

extract_module_results <- function(res_df, module_list, time_label) {
  dplyr::bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    res_df %>%
      dplyr::filter(SYMBOL %in% genes) %>%
      dplyr::mutate(Module = mod, Time = time_label)
  }))
}

lipid_res_24 <- extract_module_results(res_24_df, lipid_modules, "24h")
lipid_res_48 <- extract_module_results(res_48_df, lipid_modules, "48h")
lipid_res_72 <- extract_module_results(res_72_df, lipid_modules, "72h")

write.xlsx(
  list(
    Lipid_focus_24h = lipid_res_24,
    Lipid_focus_48h = lipid_res_48,
    Lipid_focus_72h = lipid_res_72
  ),
  file = file.path(outdir, "Focused_lipid_gene_results.xlsx"),
  overwrite = TRUE
)

############################################################
### 21F. Focused heatmap of lipid / LD / stress genes
############################################################

all_lipid_genes <- unique(unlist(lipid_modules))

lipid_annot <- gene_annot %>%
  dplyr::filter(SYMBOL %in% all_lipid_genes) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

if (nrow(lipid_annot) > 0) {
  lipid_mat <- assay(vsd)[lipid_annot$ENSEMBL, , drop = FALSE]
  rownames(lipid_mat) <- lipid_annot$SYMBOL[match(rownames(lipid_mat), lipid_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    dplyr::arrange(Time, Treatment) %>%
    dplyr::pull(Sample_name)
  
  lipid_mat <- lipid_mat[, ordered_samples, drop = FALSE]
  anno_lipid <- meta[ordered_samples, c("Time", "Treatment"), drop = FALSE]
  
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

#just 24 and 48 hours
all_lipid_genes <- unique(unlist(lipid_modules))

lipid_annot <- gene_annot %>%
  dplyr::filter(SYMBOL %in% all_lipid_genes) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

if (nrow(lipid_annot) > 0) {
  
  meta_subset <- meta %>%
    dplyr::filter(Time %in% c("24h", "48h"))
  
  keep_genes <- intersect(lipid_annot$ENSEMBL, rownames(assay(vsd)))
  ordered_samples <- meta_subset %>%
    dplyr::arrange(Time, Treatment) %>%
    dplyr::pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(assay(vsd)))
  
  lipid_mat <- assay(vsd)[keep_genes, ordered_samples, drop = FALSE]
  
  rownames(lipid_mat) <- lipid_annot$SYMBOL[
    match(rownames(lipid_mat), lipid_annot$ENSEMBL)
  ]
  
  anno_lipid <- meta_subset %>%
    dplyr::select(Sample_name, Time, Treatment) %>%
    as.data.frame()
  
  rownames(anno_lipid) <- anno_lipid$Sample_name
  anno_lipid$Sample_name <- NULL
  anno_lipid <- anno_lipid[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Focused_lipid_heatmap_24_48h.png"),
      width = 1600, height = 1300, res = 160)
  
  pheatmap(
    lipid_mat,
    scale = "row",
    annotation_col = anno_lipid,
    cluster_cols = FALSE,
    main = "Focused lipid / LD / lysosome / stress / IFN genes (24h, 48h)"
  )
  
  dev.off()
}

cairo_ps(
  file.path(outdir, "Focused_lipid_heatmap_24_48h.eps"),
  width = 10,
  height = 8
)

pheatmap(
  lipid_mat,
  scale = "row",
  annotation_col = anno_lipid,
  cluster_cols = FALSE,
  main = "Focused lipid / LD / lysosome / stress / IFN genes (24h, 48h)"
)

dev.off()

############################################################
### 21G. Module-level summary
############################################################

summarize_modules <- function(res_df, module_list, time_label) {
  dplyr::bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    sub <- res_df %>% dplyr::filter(SYMBOL %in% genes)
    
    data.frame(
      Module = mod,
      Time = time_label,
      n_genes_found = nrow(sub),
      mean_log2FC = ifelse(nrow(sub) > 0, mean(sub$log2FoldChange, na.rm = TRUE), NA),
      median_log2FC = ifelse(nrow(sub) > 0, median(sub$log2FoldChange, na.rm = TRUE), NA),
      n_up_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange > 1), 0),
      n_down_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange < -1), 0)
    )
  }))
}

module_summary_24 <- summarize_modules(res_24_df, lipid_modules, "24h")
module_summary_48 <- summarize_modules(res_48_df, lipid_modules, "48h")
module_summary_72 <- summarize_modules(res_72_df, lipid_modules, "72h")

module_summary_all <- dplyr::bind_rows(module_summary_24, module_summary_48, module_summary_72)

write.xlsx(
  module_summary_all,
  file = file.path(outdir, "Lipid_module_summary.xlsx"),
  overwrite = TRUE
)

p_module_mean <- ggplot(module_summary_all, aes(x = Module, y = mean_log2FC, fill = Time)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Lipid / LD module summary", x = "", y = "Mean log2FC")

ggsave(file.path(outdir, "Lipid_module_summary_barplot.png"),
       p_module_mean, width = 10, height = 6.5, dpi = 300)

p_module_sig <- ggplot(module_summary_all, aes(x = Module, y = n_up_sig, fill = Time)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Significant upregulated genes per module", x = "", y = "Count")

ggsave(file.path(outdir, "Lipid_module_sig_up_barplot.png"),
       p_module_sig, width = 10, height = 6.5, dpi = 300)

#just 24 and 48 hours
summarize_modules <- function(res_df, module_list, time_label) {
  dplyr::bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    sub <- res_df %>% dplyr::filter(SYMBOL %in% genes)
    
    data.frame(
      Module = mod,
      Time = time_label,
      n_genes_found = nrow(sub),
      mean_log2FC = ifelse(nrow(sub) > 0, mean(sub$log2FoldChange, na.rm = TRUE), NA),
      median_log2FC = ifelse(nrow(sub) > 0, median(sub$log2FoldChange, na.rm = TRUE), NA),
      n_up_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange > 1), 0),
      n_down_sig = ifelse(nrow(sub) > 0, sum(!is.na(sub$padj) & sub$padj < 0.05 & sub$log2FoldChange < -1), 0)
    )
  }))
}

module_summary_24 <- summarize_modules(res_24_df, lipid_modules, "24h")
module_summary_48 <- summarize_modules(res_48_df, lipid_modules, "48h")

module_summary_24_48 <- dplyr::bind_rows(module_summary_24, module_summary_48)

write.xlsx(
  module_summary_24_48,
  file = file.path(outdir, "Lipid_module_summary_24_48h.xlsx"),
  overwrite = TRUE
)

p_module_mean <- ggplot(module_summary_24_48, aes(x = Module, y = mean_log2FC, fill = Time)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Lipid / LD module summary (24h, 48h)", x = "", y = "Mean log2FC")

ggsave(file.path(outdir, "Lipid_module_summary_barplot_24_48h.png"),
       p_module_mean, width = 10, height = 6.5, dpi = 300)

ggsave(
  file.path(outdir, "Lipid_module_summary_barplot_24_48h.eps"),
  p_module_mean,
  width = 10,
  height = 6.5,
  device = cairo_ps
)

p_module_sig <- ggplot(module_summary_24_48, aes(x = Module, y = n_up_sig, fill = Time)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Significant upregulated genes per module (24h, 48h)", x = "", y = "Count")

ggsave(file.path(outdir, "Lipid_module_sig_up_barplot_24_48h.png"),
       p_module_sig, width = 10, height = 6.5, dpi = 300)

ggsave(
  file.path(outdir, "Lipid_module_sig_up_barplot_24_48h.eps"),
  p_module_mean,
  width = 10,
  height = 6.5,
  device = cairo_ps
)

############################################################
### 21H. STAT1-axis focused export
############################################################

stat1_axis_genes <- c(
  "STAT1", "STAT2", "IRF1", "IRF7", "ISG15", "IFIT1", "IFIT3", "MX1", "MX2",
  "PLIN1", "PLIN2", "PLIN3",
  "ACAT1", "ACAT2", "SOAT1", "SOAT2", "DGAT1", "DGAT2",
  "PNPLA2", "ABHD5", "LIPE", "MGLL",
  "ACSL1", "ACSL3", "CPT1A", "CPT2", "PPARGC1A",
  "LIPA", "LAMP1", "LAMP2",
  "DDIT3", "ATF4", "XBP1",
  "NFKBIA", "CXCL8", "CCL2"
)

extract_focus_panel <- function(res_df, panel_genes, time_label) {
  res_df %>%
    dplyr::filter(SYMBOL %in% panel_genes) %>%
    dplyr::select(ENSEMBL, SYMBOL, GENENAME, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    dplyr::arrange(match(SYMBOL, panel_genes)) %>%
    dplyr::mutate(Time = time_label)
}

focus_24 <- extract_focus_panel(res_24_df, stat1_axis_genes, "24h")
focus_48 <- extract_focus_panel(res_48_df, stat1_axis_genes, "48h")
focus_72 <- extract_focus_panel(res_72_df, stat1_axis_genes, "72h")

write.xlsx(
  list(
    STAT1_focus_24h = focus_24,
    STAT1_focus_48h = focus_48,
    STAT1_focus_72h = focus_72
  ),
  file = file.path(outdir, "STAT1_lipid_focus_panel.xlsx"),
  overwrite = TRUE
)

############################################################
### 21I. Targeted panel using your specific genes
############################################################

custom_qpcr_panel <- c(
  "PLIN2", "PLIN3", "ACAT1", "ACAT2", "LPIN1", "PNPLA2", "LIPE",
  "ABHD5", "LIPA", "CPT1A", "ACSL1", "ACSL3", "PPARGC1A",
  "DDIT3", "LAMP1", "STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3", "MX1", "MX2", "NFKBIA", "CXCL8"
)

extract_qpcr_panel <- function(res_df, genes, time_label) {
  res_df %>%
    dplyr::filter(SYMBOL %in% genes) %>%
    dplyr::select(ENSEMBL, SYMBOL, GENENAME, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    dplyr::arrange(match(SYMBOL, genes)) %>%
    dplyr::mutate(Time = time_label)
}

qpcr_panel_24 <- extract_qpcr_panel(res_24_df, custom_qpcr_panel, "24h")
qpcr_panel_48 <- extract_qpcr_panel(res_48_df, custom_qpcr_panel, "48h")
qpcr_panel_72 <- extract_qpcr_panel(res_72_df, custom_qpcr_panel, "72h")

write.xlsx(
  list(
    qPCR_panel_24h = qpcr_panel_24,
    qPCR_panel_48h = qpcr_panel_48,
    qPCR_panel_72h = qpcr_panel_72
  ),
  file = file.path(outdir, "Custom_lipid_qPCR_panel.xlsx"),
  overwrite = TRUE
)

############################################################
### 21J. Targeted panel heatmap
############################################################

panel_annot <- gene_annot %>%
  dplyr::filter(SYMBOL %in% custom_qpcr_panel) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

if (nrow(panel_annot) > 0) {
  panel_mat <- assay(vsd)[panel_annot$ENSEMBL, , drop = FALSE]
  rownames(panel_mat) <- panel_annot$SYMBOL[match(rownames(panel_mat), panel_annot$ENSEMBL)]
  
  ordered_samples <- meta %>%
    dplyr::arrange(Time, Treatment) %>%
    dplyr::pull(Sample_name)
  
  panel_mat <- panel_mat[, ordered_samples, drop = FALSE]
  anno_panel <- meta[ordered_samples, c("Time", "Treatment"), drop = FALSE]
  
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
### 21K. Optional: GSEA enrichment plots for selected pathways
############################################################

plot_selected_gseapaths <- function(gsea_obj, pattern_vec, prefix) {
  if (is.null(gsea_obj)) return(NULL)
  gsea_df <- as.data.frame(gsea_obj)
  if (nrow(gsea_df) == 0) return(NULL)
  
  for (pat in pattern_vec) {
    idx <- grep(pat, gsea_df$Description, ignore.case = TRUE)
    if (length(idx) > 0) {
      first_id <- gsea_df$ID[idx[1]]
      p <- enrichplot::gseaplot2(gsea_obj, geneSetID = first_id, title = gsea_df$Description[idx[1]])
      ggsave(
        file.path(outdir, paste0(prefix, "_", gsub("[^A-Za-z0-9]+", "_", pat), ".png")),
        p, width = 8, height = 6, dpi = 300
      )
    }
  }
}

plot_selected_gseapaths(
  gsea_react_24,
  pattern_vec = c("interferon", "immune", "cholesterol", "fatty acid", "lysosome", "autophagy", "ER"),
  prefix = "Reactome_GSEA_24h_selected"
)

plot_selected_gseapaths(
  gsea_react_48,
  pattern_vec = c("interferon", "immune", "cholesterol", "fatty acid", "lysosome", "autophagy", "ER"),
  prefix = "Reactome_GSEA_48h_selected"
)

plot_selected_gseapaths(
  gsea_react_72,
  pattern_vec = c("interferon", "immune", "cholesterol", "fatty acid", "lysosome", "autophagy", "ER"),
  prefix = "Reactome_GSEA_72h_selected"
)

############################################################
### 21L. Quick interpretation helper table
############################################################

interpret_panel_direction <- function(panel_df) {
  panel_df %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 1  ~ "Up_sig",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down_sig",
        !is.na(log2FoldChange) & log2FoldChange > 0 ~ "Up_nonsig",
        !is.na(log2FoldChange) & log2FoldChange < 0 ~ "Down_nonsig",
        TRUE ~ "NA"
      )
    )
}

qpcr_panel_24_interpret <- interpret_panel_direction(qpcr_panel_24)
qpcr_panel_48_interpret <- interpret_panel_direction(qpcr_panel_48)
qpcr_panel_72_interpret <- interpret_panel_direction(qpcr_panel_72)

write.xlsx(
  list(
    qPCR_panel_interpret_24h = qpcr_panel_24_interpret,
    qPCR_panel_interpret_48h = qpcr_panel_48_interpret,
    qPCR_panel_interpret_72h = qpcr_panel_72_interpret
  ),
  file = file.path(outdir, "Custom_lipid_qPCR_panel_interpretation.xlsx"),
  overwrite = TRUE
)

cat("Reactome / GSEA / lipid-focused add-on complete.\n")

