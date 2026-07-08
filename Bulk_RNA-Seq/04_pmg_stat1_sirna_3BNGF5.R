############################################################
### Primary mouse microglia STAT1 siRNA RNA-seq analysis
### Plasmidsaurus order: 3BNGF5
### Species: Mouse / Mus musculus
###
### Main analysis excludes untreated from DE contrasts.
### Untreated is retained only for baseline identity/reference QC.
############################################################

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(openxlsx)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(patchwork)

############################################################
### 0. Input paths
############################################################

base_dir <- "C:/Users/au731006/Downloads/Plasmidsaurus"

matrix_file <- file.path(
  base_dir,
  "3BNGF5_results",
  "3BNGF5-expression-matrix.tsv"
)

sample_map_file <- file.path(
  base_dir,
  "3BNGF5_results",
  "3BNGF5.csv"
)

outdir <- file.path(
  base_dir,
  "3BNGF5_PrimaryMicroglia_STAT1_siRNA_analysis"
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
### 1. Read matrix and sample map
############################################################

expr_df <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)

if (file.exists(sample_map_file)) {
  sample_map <- read.csv(sample_map_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  sample_map <- data.frame(
    run_id = paste0("3BNGF5_", 1:8),
    Sample_name = c(
      "PMG_Untreated_48-72h_R1",
      "PMG_Untreated_48-72h_R2",
      "PMG_LipoOnly_48-72h_R1",
      "PMG_LipoOnly_48-72h_R2",
      "PMG_siSTAT1_1_48-72h_R1",
      "PMG_siSTAT1_1_48-72h_R2",
      "PMG_siSTAT1_3_48-72h_R1",
      "PMG_siSTAT1_3_48-72h_R2"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    sample_map,
    file = file.path(outdir, "Generated_3BNGF5_sample_map.csv"),
    row.names = FALSE
  )
}

stopifnot(all(c("run_id", "Sample_name") %in% colnames(sample_map)))
stopifnot(all(c("gene_id", "gene_name", "gene_biotype") %in% colnames(expr_df)))

############################################################
### 2. Extract count and CPM columns
############################################################

count_cols <- grep("_count$", colnames(expr_df), value = TRUE)
cpm_cols   <- grep("_cpm$", colnames(expr_df), value = TRUE)

if (length(count_cols) == 0) {
  stop("No *_count columns found in the expression matrix.")
}

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
colnames(count_mat) <- count_run_ids

count_mat <- count_mat[!duplicated(rownames(count_mat)), , drop = FALSE]

mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

############################################################
### 4. Optional CPM matrix
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
### 5. Match run IDs to biological sample names
############################################################

missing_map <- setdiff(colnames(count_mat), sample_map$run_id)

if (length(missing_map) > 0) {
  stop("These run IDs are missing from sample_map: ", paste(missing_map, collapse = ", "))
}

idx <- match(colnames(count_mat), sample_map$run_id)
sample_map <- sample_map[idx, , drop = FALSE]

stopifnot(all(sample_map$run_id == colnames(count_mat)))

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
### 6. Parse metadata
############################################################

parse_sample_metadata <- function(sample_names) {
  tibble(Sample_name = sample_names) %>%
    mutate(
      Condition = case_when(
        grepl("^PMG_Untreated_48-72h_R[0-9]+$", Sample_name) ~ "Untreated",
        grepl("^PMG_LipoOnly_48-72h_R[0-9]+$", Sample_name) ~ "LipoOnly",
        grepl("^PMG_siSTAT1_1_48-72h_R[0-9]+$", Sample_name) ~ "siSTAT1_1",
        grepl("^PMG_siSTAT1_3_48-72h_R[0-9]+$", Sample_name) ~ "siSTAT1_3",
        TRUE ~ "Other"
      ),
      TransfectionContext = case_when(
        Condition == "Untreated" ~ "Untreated",
        Condition == "LipoOnly" ~ "LipoOnly",
        Condition %in% c("siSTAT1_1", "siSTAT1_3") ~ "RNAiMAX_siRNA",
        TRUE ~ "Other"
      ),
      siRNA = case_when(
        Condition %in% c("siSTAT1_1", "siSTAT1_3") ~ "siSTAT1",
        Condition == "LipoOnly" ~ "Control",
        Condition == "Untreated" ~ "Untreated",
        TRUE ~ "Other"
      ),
      Duplex = case_when(
        Condition == "Untreated" ~ "None",
        Condition == "LipoOnly" ~ "None",
        Condition == "siSTAT1_1" ~ "DsiRNA_1",
        Condition == "siSTAT1_3" ~ "DsiRNA_3",
        TRUE ~ "Other"
      ),
      Replicate = sub(".*_R([0-9]+)$", "R\\1", Sample_name)
    ) %>%
    mutate(
      Condition = factor(
        Condition,
        levels = c("Untreated", "LipoOnly", "siSTAT1_1", "siSTAT1_3")
      ),
      TransfectionContext = factor(
        TransfectionContext,
        levels = c("Untreated", "LipoOnly", "RNAiMAX_siRNA")
      ),
      siRNA = factor(siRNA, levels = c("Untreated", "Control", "siSTAT1")),
      Duplex = factor(Duplex, levels = c("None", "DsiRNA_1", "DsiRNA_3")),
      Replicate = factor(Replicate)
    )
}

meta_all <- parse_sample_metadata(colnames(count_mat))
rownames(meta_all) <- meta_all$Sample_name

stopifnot(all(colnames(count_mat) == rownames(meta_all)))

print(meta_all)
print(table(meta_all$Condition))

if (any(is.na(meta_all$Condition)) || any(meta_all$Condition == "Other")) {
  stop("Some samples were parsed incorrectly. Check sample names.")
}

write.xlsx(
  meta_all,
  file = file.path(outdir, "Parsed_sample_metadata_all_samples.xlsx"),
  overwrite = TRUE
)

############################################################
### 7. Save matrices with all samples
############################################################

write.xlsx(
  as.data.frame(count_mat) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "Counts_matrix_all_samples.xlsx"),
  overwrite = TRUE
)

if (!is.null(cpm_mat)) {
  write.xlsx(
    as.data.frame(cpm_mat) %>% rownames_to_column("ENSEMBL"),
    file = file.path(outdir, "CPM_matrix_all_samples.xlsx"),
    overwrite = TRUE
  )
}

############################################################
### 8. Mouse annotation
############################################################

all_ensembl <- rownames(count_mat)

gene_annot_all <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_ensembl,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

############################################################
### 9. QC model with all samples
### Untreated is used only as baseline/reference QC.
############################################################

dds_all <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_all,
  design = ~ Condition
)

keep_all <- rowSums(counts(dds_all) >= 10) >= 2
dds_all <- dds_all[keep_all, ]

dds_all <- DESeq(dds_all)
vsd_all <- vst(dds_all, blind = FALSE)

pcaData_all <- plotPCA(vsd_all, intgroup = c("Condition"), returnData = TRUE)
percentVar_all <- round(100 * attr(pcaData_all, "percentVar"))

p_pca_all <- ggplot(
  pcaData_all,
  aes(PC1, PC2, color = Condition, label = name)
) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
  ggtitle("PCA - Primary mouse microglia STAT1 siRNA, all samples") +
  theme_bw()

ggsave(
  file.path(outdir, "PCA_all_samples_including_untreated.png"),
  p_pca_all,
  width = 7,
  height = 5,
  dpi = 300
)

sample_dists_all <- dist(t(assay(vsd_all)))
sample_dist_mat_all <- as.matrix(sample_dists_all)

annotation_all <- as.data.frame(colData(dds_all)[, c("Condition"), drop = FALSE])
rownames(annotation_all) <- colnames(dds_all)

png(file.path(outdir, "Sample_distance_heatmap_all_samples.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sample_dist_mat_all,
  annotation_col = annotation_all,
  annotation_row = annotation_all,
  main = "Sample-to-sample distances: all samples"
)
dev.off()

############################################################
### 10. Main analysis: exclude untreated from DE contrasts
############################################################

main_samples <- meta_all %>%
  filter(Condition %in% c("LipoOnly", "siSTAT1_1", "siSTAT1_3")) %>%
  pull(Sample_name)

count_main <- count_mat[, main_samples, drop = FALSE]

meta_main <- meta_all %>%
  filter(Sample_name %in% main_samples) %>%
  arrange(match(Sample_name, main_samples)) %>%
  as.data.frame()

rownames(meta_main) <- meta_main$Sample_name

meta_main$Condition <- droplevels(meta_main$Condition)
meta_main$siRNA <- droplevels(meta_main$siRNA)
meta_main$Duplex <- droplevels(meta_main$Duplex)
meta_main$TransfectionContext <- droplevels(meta_main$TransfectionContext)

if (!is.null(cpm_mat)) {
  cpm_main <- cpm_mat[, main_samples, drop = FALSE]
} else {
  cpm_main <- NULL
}

print(colnames(count_main))
print(rownames(meta_main))

stopifnot(all(colnames(count_main) == rownames(meta_main)))

write.xlsx(
  meta_main,
  file = file.path(outdir, "Parsed_sample_metadata_main_analysis_excluding_untreated.xlsx"),
  overwrite = TRUE
)

dds <- DESeqDataSetFromMatrix(
  countData = count_main,
  colData = meta_main,
  design = ~ Condition
)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

gene_annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = rownames(dds),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

############################################################
### 11. Main DE contrasts
############################################################

make_res_table <- function(res_obj, annot_tbl) {
  as.data.frame(res_obj) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(annot_tbl, by = "ENSEMBL") %>%
    arrange(padj)
}

res_siSTAT1_1_vs_Lipo <- results(
  dds,
  contrast = c("Condition", "siSTAT1_1", "LipoOnly")
)

res_siSTAT1_3_vs_Lipo <- results(
  dds,
  contrast = c("Condition", "siSTAT1_3", "LipoOnly")
)

res_siSTAT1_3_vs_1 <- results(
  dds,
  contrast = c("Condition", "siSTAT1_3", "siSTAT1_1")
)

res_siSTAT1_1_vs_Lipo_df <- make_res_table(res_siSTAT1_1_vs_Lipo, gene_annot)
res_siSTAT1_3_vs_Lipo_df <- make_res_table(res_siSTAT1_3_vs_Lipo, gene_annot)
res_siSTAT1_3_vs_1_df    <- make_res_table(res_siSTAT1_3_vs_1, gene_annot)

write.xlsx(
  list(
    siSTAT1_1_vs_LipoOnly = res_siSTAT1_1_vs_Lipo_df,
    siSTAT1_3_vs_LipoOnly = res_siSTAT1_3_vs_Lipo_df,
    siSTAT1_3_vs_siSTAT1_1 = res_siSTAT1_3_vs_1_df
  ),
  file = file.path(outdir, "DE_results_STAT1_siRNA_vs_LipoOnly.xlsx"),
  overwrite = TRUE
)

############################################################
### 12. Optional pooled STAT1 siRNA analysis
############################################################

meta_pooled <- meta_main
meta_pooled$siRNA <- factor(
  ifelse(meta_pooled$Condition == "LipoOnly", "Control", "siSTAT1"),
  levels = c("Control", "siSTAT1")
)

dds_pooled <- DESeqDataSetFromMatrix(
  countData = count_main,
  colData = meta_pooled,
  design = ~ siRNA
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
    siSTAT1_pooled_vs_LipoOnly = res_siSTAT1_pooled_vs_Lipo_df
  ),
  file = file.path(outdir, "DE_results_STAT1_siRNA_pooled_vs_LipoOnly.xlsx"),
  overwrite = TRUE
)

contrast_results <- list(
  siSTAT1_1_vs_LipoOnly = res_siSTAT1_1_vs_Lipo_df,
  siSTAT1_3_vs_LipoOnly = res_siSTAT1_3_vs_Lipo_df,
  siSTAT1_3_vs_siSTAT1_1 = res_siSTAT1_3_vs_1_df,
  siSTAT1_pooled_vs_LipoOnly = res_siSTAT1_pooled_vs_Lipo_df
)

############################################################
### 13. Normalized counts and VST
############################################################

vsd <- vst(dds, blind = FALSE)

norm_counts <- counts(dds, normalized = TRUE)

write.xlsx(
  as.data.frame(norm_counts) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "Normalized_counts_main_analysis.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  as.data.frame(assay(vsd)) %>% rownames_to_column("ENSEMBL"),
  file = file.path(outdir, "VST_matrix_main_analysis.xlsx"),
  overwrite = TRUE
)

############################################################
### 14. PCA and sample distance for main analysis
############################################################

pcaData <- plotPCA(vsd, intgroup = c("Condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(
  pcaData,
  aes(PC1, PC2, color = Condition, label = name)
) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Primary mouse microglia STAT1 siRNA, untreated excluded") +
  theme_bw()

ggsave(
  file.path(outdir, "PCA_main_analysis_excluding_untreated.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

sample_dists <- dist(t(assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)

annotation_df <- as.data.frame(colData(dds)[, c("Condition"), drop = FALSE])
rownames(annotation_df) <- colnames(dds)

png(file.path(outdir, "Sample_distance_heatmap_main_analysis.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sample_dist_mat,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "Sample-to-sample distances: untreated excluded"
)
dev.off()

############################################################
### 15. Volcano plots
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
  siSTAT1_1_vs_LipoOnly = "Primary microglia: STAT1 siRNA #1 vs LipoOnly",
  siSTAT1_3_vs_LipoOnly = "Primary microglia: STAT1 siRNA #3 vs LipoOnly",
  siSTAT1_3_vs_siSTAT1_1 = "Primary microglia: STAT1 siRNA #3 vs #1",
  siSTAT1_pooled_vs_LipoOnly = "Primary microglia: pooled STAT1 siRNA vs LipoOnly"
)

for (nm in names(contrast_results)) {
  p <- plot_volcano(contrast_results[[nm]], volcano_titles[[nm]])

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
### 16. Top DEGs
############################################################

get_top_genes <- function(res_df, n = 20) {
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
  tops <- get_top_genes(contrast_results[[nm]], n = 20)
  top_deg_list[[paste0(nm, "_up")]] <- tops$up
  top_deg_list[[paste0(nm, "_down")]] <- tops$down
}

write.xlsx(
  top_deg_list,
  file = file.path(outdir, "Top_DEGs.xlsx"),
  overwrite = TRUE
)

############################################################
### 17. Focused primary microglia / STAT1 / lipid modules
############################################################

microglia_modules <- list(

  IFN_STAT1_ISG = c(
    "Stat1", "Stat2", "Irf1", "Irf7", "Irf9",
    "Isg15", "Ifit1", "Ifit2", "Ifit3",
    "Ifitm1", "Ifitm2", "Ifitm3",
    "Oas1a", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2",
    "Mx1", "Mx2", "Usp18", "Rsad2", "Zbp1", "Ddx58", "Ifih1"
  ),

  Inflammatory_cytokine_chemokine = c(
    "Tnf", "Il1b", "Il6", "Cxcl10", "Cxcl9", "Ccl2", "Ccl3", "Ccl4",
    "Nfkbia", "Nfkbiz", "Tnfaip3", "Ptgs2", "Nos2", "Il10", "Il12b"
  ),

  Antigen_presentation_MHC = c(
    "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1",
    "Cd74", "Ciita", "B2m", "Tap1", "Tap2", "Psmb8", "Psmb9"
  ),

  Microglia_identity = c(
    "Hexb", "Csf1r", "Cx3cr1", "Trem2", "Tyrobp", "Aif1",
    "P2ry12", "Tmem119", "Sall1", "Siglech", "Fcrls", "C1qa", "C1qb", "C1qc"
  ),

  Homeostatic_microglia = c(
    "P2ry12", "Tmem119", "Sall1", "Siglech", "Fcrls", "Cx3cr1",
    "Olfml3", "Gpr34", "P2ry13", "Mef2a", "Mef2c", "Sparc"
  ),

  DAM_lipid_associated = c(
    "Apoe", "Lpl", "Trem2", "Tyrobp", "Spp1", "Lgals3", "Cst7",
    "Ctsb", "Ctsd", "Ctsl", "Lyz2", "Itgax", "Clec7a", "Axl", "Cd68"
  ),

  Lipid_droplet_storage = c(
    "Plin1", "Plin2", "Plin3", "Plin4", "Plin5",
    "Dgat1", "Dgat2", "Soat1", "Soat2", "Lpin1", "Lpin2", "Fitm2"
  ),

  Lipid_uptake_processing = c(
    "Apoe", "Lpl", "Lipa", "Cd36", "Ldlr", "Lrp1", "Abca1", "Abcg1",
    "Nceh1", "Npc1", "Npc2", "Scarb1", "Fabp3", "Fabp4", "Fabp5"
  ),

  Lipolysis_FAO = c(
    "Pnpla2", "Abhd5", "Lipe", "Mgll",
    "Cpt1a", "Cpt2", "Acadm", "Acadvl", "Hadha", "Hadhb",
    "Acox1", "Ppara", "Ppargc1a"
  ),

  Lysosome_lipophagy = c(
    "Lipa", "Lamp1", "Lamp2", "Ctsb", "Ctsd", "Ctsl",
    "Rab7", "Rab7a", "Atg5", "Atg7", "Becn1", "Sqstm1", "Map1lc3b"
  ),

  OXPHOS_ETC = c(
    "Ndufa1", "Ndufa2", "Ndufa9", "Ndufb8", "Ndufs1", "Ndufs2",
    "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2",
    "Cox4i1", "Cox5a", "Cox6a1",
    "Atp5f1a", "Atp5f1b", "Atp5f1c", "Atp5mc1"
  ),

  Glycolysis = c(
    "Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkl", "Pfkp",
    "Aldoa", "Gapdh", "Pgk1", "Pgam1", "Eno1", "Pkm", "Ldha"
  ),

  Iron_redox_glutathione = c(
    "Fth1", "Ftl1", "Hmox1", "Slc40a1", "Tfrc", "Slc11a2",
    "Gpx4", "Gpx1", "Gclc", "Gclm", "Gsr", "Gss",
    "Txn1", "Txnrd1", "Prdx1", "Prdx2", "Nqo1", "Sod1", "Sod2"
  )
)

extract_module_results <- function(res_df, module_list, contrast_label) {
  bind_rows(lapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    res_df %>%
      filter(SYMBOL %in% genes) %>%
      mutate(Module = mod, Contrast = contrast_label)
  }))
}

module_result_list <- lapply(names(contrast_results), function(nm) {
  extract_module_results(contrast_results[[nm]], microglia_modules, nm)
})
names(module_result_list) <- names(contrast_results)

write.xlsx(
  module_result_list,
  file = file.path(outdir, "Focused_microglia_STAT1_lipid_module_gene_results.xlsx"),
  overwrite = TRUE
)

############################################################
### 18. Module-level summary
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
  summarize_modules(contrast_results[[nm]], microglia_modules, nm)
}))

write.xlsx(
  module_summary_all,
  file = file.path(outdir, "Focused_microglia_STAT1_lipid_module_summary.xlsx"),
  overwrite = TRUE
)

p_module_mean <- ggplot(module_summary_all, aes(x = Module, y = mean_log2FC, fill = Contrast)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Primary microglia STAT1 siRNA: focused module mean log2FC",
    x = "",
    y = "Mean log2FC"
  )

ggsave(
  file.path(outdir, "Focused_microglia_STAT1_lipid_module_summary_barplot.png"),
  p_module_mean,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
### 19. Focused heatmap
############################################################

all_module_genes <- unique(unlist(microglia_modules))

module_annot <- gene_annot %>%
  filter(SYMBOL %in% all_module_genes) %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(module_annot) > 0) {
  keep_genes <- intersect(module_annot$ENSEMBL, rownames(assay(vsd)))
  module_mat <- assay(vsd)[keep_genes, , drop = FALSE]
  rownames(module_mat) <- module_annot$SYMBOL[match(rownames(module_mat), module_annot$ENSEMBL)]

  ordered_samples <- meta_main %>%
    arrange(Condition, Replicate) %>%
    pull(Sample_name)

  ordered_samples <- intersect(ordered_samples, colnames(module_mat))
  module_mat <- module_mat[, ordered_samples, drop = FALSE]

  anno_module <- meta_main[, c("Sample_name", "Condition"), drop = FALSE] %>%
    as.data.frame()
  rownames(anno_module) <- anno_module$Sample_name
  anno_module$Sample_name <- NULL
  anno_module <- anno_module[ordered_samples, , drop = FALSE]

  png(file.path(outdir, "Focused_microglia_STAT1_lipid_module_heatmap.png"), width = 1800, height = 1600, res = 160)
  pheatmap(
    module_mat,
    scale = "row",
    annotation_col = anno_module,
    cluster_cols = FALSE,
    main = "Focused microglia / STAT1 / lipid / redox modules"
  )
  dev.off()
}

############################################################
### 20. Marker CPM summary
############################################################

marker_panel <- unique(c(
  microglia_modules$IFN_STAT1_ISG,
  microglia_modules$DAM_lipid_associated,
  microglia_modules$Lipid_droplet_storage,
  microglia_modules$Lipid_uptake_processing,
  microglia_modules$Iron_redox_glutathione,
  microglia_modules$Homeostatic_microglia
))

if (!is.null(cpm_mat)) {
  marker_annot <- gene_annot_all %>%
    filter(SYMBOL %in% marker_panel) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  keep_marker_genes <- intersect(marker_annot$ENSEMBL, rownames(cpm_mat))

  marker_cpm_long <- as.data.frame(cpm_mat[keep_marker_genes, , drop = FALSE]) %>%
    rownames_to_column("ENSEMBL") %>%
    pivot_longer(
      cols = -ENSEMBL,
      names_to = "Sample_name",
      values_to = "CPM"
    ) %>%
    left_join(marker_annot, by = "ENSEMBL") %>%
    left_join(meta_all %>% dplyr::select(Sample_name, Condition, Replicate), by = "Sample_name")

  marker_cpm_summary <- marker_cpm_long %>%
    group_by(SYMBOL, GENENAME, Condition) %>%
    summarise(
      mean_CPM = mean(CPM, na.rm = TRUE),
      sd_CPM = sd(CPM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = c(mean_CPM, sd_CPM)
    )

  write.xlsx(
    list(
      marker_CPM_long = marker_cpm_long,
      marker_CPM_summary = marker_cpm_summary
    ),
    file = file.path(outdir, "Focused_marker_CPM_summary_all_samples.xlsx"),
    overwrite = TRUE
  )
}

############################################################
### 21. GO enrichment
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
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

go_objects <- list()
go_export <- list()

for (nm in names(contrast_results)) {
  ego_up <- run_go(contrast_results[[nm]], "up")
  ego_down <- run_go(contrast_results[[nm]], "down")

  go_objects[[paste0(nm, "_up")]] <- ego_up
  go_objects[[paste0(nm, "_down")]] <- ego_down

  if (!is.null(ego_up)) go_export[[paste0("GO_", nm, "_up")]] <- as.data.frame(ego_up)
  if (!is.null(ego_down)) go_export[[paste0("GO_", nm, "_down")]] <- as.data.frame(ego_down)
}

make_safe_sheet_names <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- substr(x, 1, 31)
  
  # Ensure uniqueness after truncation
  if (any(duplicated(x))) {
    x <- make.unique(x, sep = "_")
    x <- substr(x, 1, 31)
  }
  
  return(x)
}

write_xlsx_safe <- function(x, file) {
  names(x) <- make_safe_sheet_names(names(x))
  openxlsx::write.xlsx(x, file = file, overwrite = TRUE)
}

if (length(go_export) > 0) {
  write_xlsx_safe(
    go_export,
    file.path(outdir, "GO_enrichment_results.xlsx")
  )
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
}

############################################################
### 22. Reactome GSEA
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
    organism = "mouse",
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
    }
  }
}

pattern_vec <- c(
  "interferon",
  "cytokine",
  "antigen",
  "cholesterol",
  "fatty acid",
  "lipid",
  "lysosome",
  "autophagy",
  "oxidative",
  "mitochond",
  "respiratory",
  "iron",
  "glutathione"
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
### 23. Final summary
############################################################

message("Primary mouse microglia STAT1 siRNA analysis completed successfully.")
message("Output directory: ", outdir)
message("All samples processed for QC: ", ncol(count_mat))
message("Main-analysis samples, untreated excluded: ", ncol(count_main))
message("Genes retained in main DESeq2 model: ", nrow(dds))
