############################################################
### Build SYMBOL-level VST matrix safely
############################################################

vsd_mat <- assay(vsd_group)

annot_for_vsd <- gene_annot %>%
  dplyr::filter(ENSEMBL %in% rownames(vsd_mat), !is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

vsd_df <- as.data.frame(vsd_mat, check.names = FALSE) %>%
  tibble::rownames_to_column("ENSEMBL") %>%
  dplyr::inner_join(
    annot_for_vsd[, c("ENSEMBL", "SYMBOL")],
    by = "ENSEMBL"
  )

# Identify sample columns explicitly
sample_cols <- setdiff(colnames(vsd_df), c("ENSEMBL", "SYMBOL"))

# Force sample columns to numeric
vsd_df[sample_cols] <- lapply(vsd_df[sample_cols], function(x) as.numeric(as.character(x)))

# If multiple Ensembl IDs map to same SYMBOL, keep the one with highest mean expression
vsd_df <- vsd_df %>%
  dplyr::mutate(mean_expr = rowMeans(dplyr::select(., dplyr::all_of(sample_cols)), na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(mean_expr)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# Build numeric matrix
vsd_symbol_mat <- as.matrix(vsd_df[, sample_cols, drop = FALSE])
rownames(vsd_symbol_mat) <- vsd_df$SYMBOL
storage.mode(vsd_symbol_mat) <- "numeric"

# sanity check
stopifnot(is.numeric(vsd_symbol_mat))
############################################################
### Define modules
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

############################################################
### Function: z-score based module score
############################################################

calc_module_score_z <- function(expr_mat, gene_list, min_genes = 2) {
  genes_present <- intersect(gene_list, rownames(expr_mat))
  
  if (length(genes_present) < min_genes) {
    return(list(
      score = rep(NA_real_, ncol(expr_mat)),
      genes_used = genes_present
    ))
  }
  
  submat <- expr_mat[genes_present, , drop = FALSE]
  
  # z-score each gene across samples
  zmat <- t(scale(t(submat)))
  
  # remove genes with zero variance / NA after scaling
  keep <- rowSums(is.na(zmat)) == 0
  zmat <- zmat[keep, , drop = FALSE]
  genes_used <- rownames(zmat)
  
  if (nrow(zmat) < min_genes) {
    return(list(
      score = rep(NA_real_, ncol(expr_mat)),
      genes_used = genes_used
    ))
  }
  
  module_score <- colMeans(zmat, na.rm = TRUE)
  
  return(list(
    score = module_score,
    genes_used = genes_used
  ))
}

############################################################
### Compute module scores
############################################################

module_score_list <- lapply(names(lipid_modules), function(mod) {
  out <- calc_module_score_z(vsd_symbol_mat, lipid_modules[[mod]], min_genes = 2)
  
  data.frame(
    Sample_name = colnames(vsd_symbol_mat),
    Module = mod,
    Score = out$score,
    Genes_used_n = length(out$genes_used),
    Genes_used = paste(out$genes_used, collapse = ", "),
    stringsAsFactors = FALSE
  )
})

module_scores_long <- dplyr::bind_rows(module_score_list) %>%
  dplyr::left_join(meta[, c("Sample_name", "Group", "Knockdown", "LPS")], by = "Sample_name")

write.xlsx(
  module_scores_long,
  file = file.path(outdir, "Module_scores_long.xlsx"),
  overwrite = TRUE
)

############################################################
### Plot module scores by group
############################################################

p_modules <- ggplot(module_scores_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Module scores across HMC3 STAT1/LPS conditions",
    x = "",
    y = "Module score (mean gene-wise z-score)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(outdir, "Module_scores_by_group.png"),
  p_modules,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(outdir, "Module_scores_by_group.eps"),
  p_modules,
  width = 12,
  height = 8,
  device = cairo_ps
)

#stats
############################################################
### Module-level statistical testing (2x2 design)
############################################################
library(ggpubr)
run_module_stats <- function(df) {
  # remove NA scores
  df <- df[!is.na(df$Score), ]
  
  if (nrow(df) < 4) return(NULL)
  
  fit <- lm(Score ~ Knockdown * LPS, data = df)
  
  anova_res <- as.data.frame(anova(fit))
  anova_res$Term <- rownames(anova_res)
  rownames(anova_res) <- NULL
  
  anova_res
}

module_stats_list <- lapply(split(module_scores_long, module_scores_long$Module), run_module_stats)
module_stats_list <- module_stats_list[!sapply(module_stats_list, is.null)]

# add module names
for (nm in names(module_stats_list)) {
  module_stats_list[[nm]]$Module <- nm
}

module_stats_df <- dplyr::bind_rows(module_stats_list)

write.xlsx(
  module_stats_df,
  file = file.path(outdir, "Module_ANOVA_results.xlsx"),
  overwrite = TRUE
)


#LD structure + Lipolysis + Lysosome/lipophagy
#STAT1 effect, no  LPS effect

#FAO
#STAT1 and LPS

#Esterification
#no main effects, interaction only

#IFN/STAT
#BOTH STAT1 AND LPS AFFECT IFN module

#nfkb/inflammation
#lps dominates, stat1 some contribution

#ER stress
#both factors, but #independent




############################################################
### Summarize module scores per group
############################################################

module_summary <- module_scores_long %>%
  dplyr::group_by(Module, Group) %>%
  dplyr::summarise(
    mean_score = mean(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    .groups = "drop"
  )

write.xlsx(
  module_summary,
  file = file.path(outdir, "Module_score_summary.xlsx"),
  overwrite = TRUE
)

############################################################
### Simple statistics: module score differences by group
############################################################

run_module_lm <- function(df) {
  if (all(is.na(df$Score))) return(NULL)
  
  fit <- lm(Score ~ Knockdown * LPS, data = df)
  out <- as.data.frame(anova(fit))
  out$Term <- rownames(out)
  rownames(out) <- NULL
  out
}

module_stats_list <- lapply(split(module_scores_long, module_scores_long$Module), run_module_lm)
module_stats_list <- module_stats_list[!sapply(module_stats_list, is.null)]

for (nm in names(module_stats_list)) {
  module_stats_list[[nm]]$Module <- nm
}

module_stats_df <- dplyr::bind_rows(module_stats_list)

write.xlsx(
  module_stats_df,
  file = file.path(outdir, "Module_score_ANOVA.xlsx"),
  overwrite = TRUE
)

############################################################
### Heatmap of module genes
############################################################

all_module_genes <- unique(unlist(lipid_modules))
genes_present <- intersect(all_module_genes, rownames(vsd_symbol_mat))

if (length(genes_present) > 1) {
  heat_mat <- vsd_symbol_mat[genes_present, , drop = FALSE]
  
  ordered_samples <- meta %>%
    dplyr::arrange(Group) %>%
    dplyr::pull(Sample_name)
  
  ordered_samples <- intersect(ordered_samples, colnames(heat_mat))
  heat_mat <- heat_mat[, ordered_samples, drop = FALSE]
  
  anno_df <- meta %>%
    dplyr::select(Sample_name, Group, Knockdown, LPS) %>%
    as.data.frame()
  
  rownames(anno_df) <- anno_df$Sample_name
  anno_df$Sample_name <- NULL
  anno_df <- anno_df[ordered_samples, , drop = FALSE]
  
  png(file.path(outdir, "Module_gene_heatmap.png"), width = 1600, height = 1400, res = 160)
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = anno_df,
    cluster_cols = FALSE,
    main = "Lipid / IFN / stress module genes"
  )
  dev.off()
  
  cairo_ps(file.path(outdir, "Module_gene_heatmap.eps"), width = 10, height = 9)
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = anno_df,
    cluster_cols = FALSE,
    main = "Lipid / IFN / stress module genes"
  )
  dev.off()
}
