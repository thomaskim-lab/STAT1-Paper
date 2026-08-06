# =============================================================================
# 03_promoter_distal_accessibility.R
# Summarizes pseudobulk DA into promoter/distal/locus metrics per target gene.
#
# Does not fill missing compartment values with zero.
#
# GenomeDK / Hypothalamus version.
# =============================================================================

# -----------------------------------------------------------------------------
# Robust source of 00_config.R
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- grep(file_arg, args, value = TRUE)

if (length(script_arg) > 0) {
  script_path <- normalizePath(sub(file_arg, "", script_arg[1]))
  script_dir <- dirname(script_path)
} else {
  script_dir <- getwd()
}

source(file.path(script_dir, "00_config.R"))

check_inputs(require_tools = FALSE)

message("============================================================")
message("03_promoter_distal_accessibility.R")
message("Started: ", Sys.time())
message("============================================================")

# =============================================================================
# Required inputs
# =============================================================================

compartment_map_rds <- file.path(TABLE_DIR, "peak_to_target_gene_compartments.rds")

check_file(compartment_map_rds, "peak_to_target_gene_compartments.rds")
check_file(DA_STAT1_CSV, "DA_global_STAT1_KO_vs_WT.csv")
check_file(DA_IRF1_CSV, "DA_global_IRF1_KO_vs_WT.csv")

message("Reading compartment map: ", compartment_map_rds)
compartment_map <- readRDS(compartment_map_rds)

required_comp_cols <- c("peak", "gene_name", "compartment")
missing_comp_cols <- setdiff(required_comp_cols, colnames(compartment_map))

if (length(missing_comp_cols) > 0) {
  stop(
    "Compartment map is missing required columns: ",
    paste(missing_comp_cols, collapse = ", ")
  )
}

compartment_map <- compartment_map %>%
  dplyr::distinct(
    peak,
    gene_name,
    compartment,
    .keep_all = TRUE
  )

message("Compartment map rows after distinct: ", nrow(compartment_map))
message("Unique target genes in compartment map: ", dplyr::n_distinct(compartment_map$gene_name))
message("Unique peaks in compartment map: ", dplyr::n_distinct(compartment_map$peak))

# =============================================================================
# Read DA tables
# =============================================================================

message("Reading DA tables...")

da_stat1 <- read_da_table(
  DA_STAT1_CSV,
  comparison_name = "STAT1_KO_vs_WT"
)

da_irf1 <- read_da_table(
  DA_IRF1_CSV,
  comparison_name = "IRF1_KO_vs_WT"
)

da_all <- dplyr::bind_rows(da_stat1, da_irf1)

required_da_cols <- c("peak", "logFC", "FDR", "comparison")
missing_da_cols <- setdiff(required_da_cols, colnames(da_all))

if (length(missing_da_cols) > 0) {
  stop(
    "DA tables are missing required columns: ",
    paste(missing_da_cols, collapse = ", ")
  )
}

da_all <- da_all %>%
  dplyr::mutate(
    logFC = as.numeric(logFC),
    FDR = as.numeric(FDR),
    PValue = as.numeric(PValue),
    logCPM = as.numeric(logCPM)
  )

message("DA rows total: ", nrow(da_all))
message("DA comparisons: ", paste(unique(da_all$comparison), collapse = ", "))

# =============================================================================
# Join peak DA to target-gene compartments
# =============================================================================

message("Joining DA peaks to target-gene compartments...")

# Keep only the compartment-map columns we want to avoid gene_name.x / gene_name.y
# after joining. The gene_name used here must come from 02_define_compartments.R.
compartment_map_for_join <- compartment_map %>%
  dplyr::select(
    peak,
    gene_name,
    compartment,
    dplyr::everything()
  ) %>%
  dplyr::distinct()

# Drop any gene_name-like columns from DA tables before joining.
da_all_for_join <- da_all %>%
  dplyr::select(
    -dplyr::any_of(c("gene_name", "gene", "nearest_gene", "symbol"))
  )

da_compartment <- da_all_for_join %>%
  dplyr::inner_join(
    compartment_map_for_join,
    by = "peak",
    relationship = "many-to-many"
  ) %>%
  dplyr::distinct(
    peak,
    gene_name,
    comparison,
    compartment,
    .keep_all = TRUE
  )

# Keep original downstream object name used by the rest of this script.
da_comp <- da_compartment

if (nrow(da_comp) == 0) {
  stop("No DA peaks overlapped target-gene compartments.")
}

da_comp_path <- file.path(TABLE_DIR, "DA_peaks_in_target_gene_compartments.csv")
readr::write_csv(da_comp, da_comp_path)

message("Wrote joined DA-compartment table: ", da_comp_path)
message("Rows in DA-compartment table: ", nrow(da_comp))
message("Genes represented: ", dplyr::n_distinct(da_comp$gene_name))
message("Peaks represented: ", dplyr::n_distinct(da_comp$peak))

# =============================================================================
# Helper functions
# =============================================================================

summarise_compartment <- function(df) {
  df %>%
    dplyr::group_by(gene_name, comparison, compartment) %>%
    dplyr::summarise(
      n_peaks = dplyr::n_distinct(peak),

      mean_logFC = safe_mean(logFC),
      median_logFC = safe_median(logFC),

      n_sig_lost_FDR10 = sum(sig_lost(logFC, FDR, FDR_CUTOFF), na.rm = TRUE),
      n_sig_gained_FDR10 = sum(sig_gained(logFC, FDR, FDR_CUTOFF), na.rm = TRUE),

      n_sig_lost_FDR05 = sum(sig_lost(logFC, FDR, STRICT_FDR_CUTOFF), na.rm = TRUE),
      n_sig_gained_FDR05 = sum(sig_gained(logFC, FDR, STRICT_FDR_CUTOFF), na.rm = TRUE),

      min_FDR = ifelse(all(is.na(FDR)), NA_real_, min(FDR, na.rm = TRUE)),
      min_PValue = ifelse(all(is.na(PValue)), NA_real_, min(PValue, na.rm = TRUE)),
      mean_logCPM = safe_mean(logCPM),

      .groups = "drop"
    )
}

ensure_column <- function(df, col, default = NA_real_) {
  if (!col %in% colnames(df)) {
    df[[col]] <- default
  }
  df
}

add_missing_wide_columns <- function(df) {
  expected_numeric_cols <- c(
    "mean_logFC_promoter_2kb",
    "mean_logFC_promoter_1kb",
    "mean_logFC_distal_TSS_100kb_noPromoter",
    "mean_logFC_locus_TSS_100kb",
    "mean_logFC_locus_gene_body_100kb",

    "median_logFC_promoter_2kb",
    "median_logFC_promoter_1kb",
    "median_logFC_distal_TSS_100kb_noPromoter",
    "median_logFC_locus_TSS_100kb",
    "median_logFC_locus_gene_body_100kb",

    "n_peaks_promoter_2kb",
    "n_peaks_promoter_1kb",
    "n_peaks_distal_TSS_100kb_noPromoter",
    "n_peaks_locus_TSS_100kb",
    "n_peaks_locus_gene_body_100kb",

    "n_sig_lost_FDR10_promoter_2kb",
    "n_sig_lost_FDR10_promoter_1kb",
    "n_sig_lost_FDR10_distal_TSS_100kb_noPromoter",
    "n_sig_lost_FDR10_locus_TSS_100kb",
    "n_sig_lost_FDR10_locus_gene_body_100kb",

    "n_sig_gained_FDR10_promoter_2kb",
    "n_sig_gained_FDR10_promoter_1kb",
    "n_sig_gained_FDR10_distal_TSS_100kb_noPromoter",
    "n_sig_gained_FDR10_locus_TSS_100kb",
    "n_sig_gained_FDR10_locus_gene_body_100kb",

    "n_sig_lost_FDR05_promoter_2kb",
    "n_sig_lost_FDR05_promoter_1kb",
    "n_sig_lost_FDR05_distal_TSS_100kb_noPromoter",
    "n_sig_lost_FDR05_locus_TSS_100kb",
    "n_sig_lost_FDR05_locus_gene_body_100kb",

    "n_sig_gained_FDR05_promoter_2kb",
    "n_sig_gained_FDR05_promoter_1kb",
    "n_sig_gained_FDR05_distal_TSS_100kb_noPromoter",
    "n_sig_gained_FDR05_locus_TSS_100kb",
    "n_sig_gained_FDR05_locus_gene_body_100kb",

    "min_FDR_promoter_2kb",
    "min_FDR_promoter_1kb",
    "min_FDR_distal_TSS_100kb_noPromoter",
    "min_FDR_locus_TSS_100kb",
    "min_FDR_locus_gene_body_100kb",

    "min_PValue_promoter_2kb",
    "min_PValue_promoter_1kb",
    "min_PValue_distal_TSS_100kb_noPromoter",
    "min_PValue_locus_TSS_100kb",
    "min_PValue_locus_gene_body_100kb",

    "mean_logCPM_promoter_2kb",
    "mean_logCPM_promoter_1kb",
    "mean_logCPM_distal_TSS_100kb_noPromoter",
    "mean_logCPM_locus_TSS_100kb",
    "mean_logCPM_locus_gene_body_100kb"
  )

  for (col in expected_numeric_cols) {
    df <- ensure_column(df, col, NA_real_)
  }

  df
}

# =============================================================================
# Summarize compartments
# =============================================================================

message("Summarizing DA by gene, comparison, and compartment...")

comp_summary_long <- summarise_compartment(da_comp)

comp_summary_long_path <- file.path(
  TABLE_DIR,
  "target_gene_compartment_DA_summary_long.csv"
)

readr::write_csv(comp_summary_long, comp_summary_long_path)

message("Wrote long compartment summary: ", comp_summary_long_path)

# =============================================================================
# Wide table for promoter vs distal/locus contrast
# =============================================================================

wide_logfc <- comp_summary_long %>%
  dplyr::select(
    gene_name,
    comparison,
    compartment,
    n_peaks,
    mean_logFC,
    median_logFC,
    n_sig_lost_FDR10,
    n_sig_gained_FDR10,
    n_sig_lost_FDR05,
    n_sig_gained_FDR05,
    min_FDR,
    min_PValue,
    mean_logCPM
  ) %>%
  tidyr::pivot_wider(
    names_from = compartment,
    values_from = c(
      n_peaks,
      mean_logFC,
      median_logFC,
      n_sig_lost_FDR10,
      n_sig_gained_FDR10,
      n_sig_lost_FDR05,
      n_sig_gained_FDR05,
      min_FDR,
      min_PValue,
      mean_logCPM
    )
  )

wide_logfc <- add_missing_wide_columns(wide_logfc)

wide_logfc <- wide_logfc %>%
  dplyr::mutate(
    redistribution_score = dplyr::if_else(
      !is.na(mean_logFC_promoter_2kb) &
        !is.na(mean_logFC_distal_TSS_100kb_noPromoter) &
        !is.na(n_peaks_promoter_2kb) &
        !is.na(n_peaks_distal_TSS_100kb_noPromoter) &
        n_peaks_promoter_2kb > 0 &
        n_peaks_distal_TSS_100kb_noPromoter > 0,
      mean_logFC_distal_TSS_100kb_noPromoter - mean_logFC_promoter_2kb,
      NA_real_
    ),

    redistribution_score_promoter1kb = dplyr::if_else(
      !is.na(mean_logFC_promoter_1kb) &
        !is.na(mean_logFC_distal_TSS_100kb_noPromoter) &
        !is.na(n_peaks_promoter_1kb) &
        !is.na(n_peaks_distal_TSS_100kb_noPromoter) &
        n_peaks_promoter_1kb > 0 &
        n_peaks_distal_TSS_100kb_noPromoter > 0,
      mean_logFC_distal_TSS_100kb_noPromoter - mean_logFC_promoter_1kb,
      NA_real_
    ),

    target_class = dplyr::case_when(
      gene_name %in% LIPID_GENES ~ "primary_lipid",
      gene_name %in% SECONDARY_LIPID_GENES ~ "secondary_lipid",
      TRUE ~ "other"
    ),

    is_lipid_target = gene_name %in% ALL_TARGET_GENES,

    interpretation_note = dplyr::case_when(
      comparison == "STAT1_KO_vs_WT" ~ "directional_STAT1_KO_n1",
      comparison == "IRF1_KO_vs_WT" ~ "small_n_IRF1_KO_vs_WT",
      TRUE ~ "interpret_with_caution"
    )
  ) %>%
  dplyr::arrange(comparison, target_class, gene_name)

wide_logfc_path <- file.path(
  TABLE_DIR,
  "target_gene_promoter_distal_redistribution_summary.csv"
)

readr::write_csv(wide_logfc, wide_logfc_path)

message("Wrote wide promoter/distal/locus summary: ", wide_logfc_path)

# =============================================================================
# Useful inspection table: most altered target-gene peaks
# =============================================================================

top_target_peak_table <- da_comp %>%
  dplyr::mutate(
    target_class = dplyr::case_when(
      gene_name %in% LIPID_GENES ~ "primary_lipid",
      gene_name %in% SECONDARY_LIPID_GENES ~ "secondary_lipid",
      TRUE ~ "other"
    ),
    direction = dplyr::case_when(
      sig_lost(logFC, FDR, FDR_CUTOFF) ~ "lost_FDR10",
      sig_gained(logFC, FDR, FDR_CUTOFF) ~ "gained_FDR10",
      TRUE ~ "not_FDR10"
    )
  ) %>%
  dplyr::arrange(comparison, FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::select(
    comparison,
    gene_name,
    target_class,
    compartment,
    peak,
    logFC,
    logCPM,
    PValue,
    FDR,
    direction,
    dplyr::everything()
  )

top_target_peak_path <- file.path(
  TABLE_DIR,
  "DA_peaks_in_target_gene_compartments_ranked.csv"
)

readr::write_csv(top_target_peak_table, top_target_peak_path)

message("Wrote ranked target-gene DA peak table: ", top_target_peak_path)

# =============================================================================
# Plots
# =============================================================================

message("Generating plots...")

plot_df <- wide_logfc %>%
  dplyr::filter(target_class %in% c("primary_lipid", "secondary_lipid")) %>%
  dplyr::filter(
    !is.na(mean_logFC_promoter_2kb),
    !is.na(mean_logFC_distal_TSS_100kb_noPromoter)
  )

if (nrow(plot_df) > 0) {
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = mean_logFC_promoter_2kb,
      y = mean_logFC_distal_TSS_100kb_noPromoter,
      label = gene_name
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ comparison) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Promoter versus distal accessibility at target lipid-handling loci",
      subtitle = "Distal compartment = TSS ±100 kb with all annotated 2 kb promoters excluded",
      x = "Promoter mean logFC (TSS ±2 kb)",
      y = "Distal mean logFC (TSS ±100 kb, promoters excluded)"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        max.overlaps = 50,
        size = 3
      )
  } else {
    warning("Package ggrepel not installed. Using geom_text instead.")
    p <- p +
      ggplot2::geom_text(
        size = 3,
        hjust = -0.05,
        vjust = 0.5,
        check_overlap = TRUE
      )
  }

  pdf_path <- file.path(PLOT_DIR, "target_genes_promoter_vs_distal_logFC.pdf")
  png_path <- file.path(PLOT_DIR, "target_genes_promoter_vs_distal_logFC.png")

  ggplot2::ggsave(pdf_path, p, width = 8, height = 5)
  ggplot2::ggsave(png_path, p, width = 8, height = 5, dpi = 300)

  message("Wrote plot: ", pdf_path)
  message("Wrote plot: ", png_path)
} else {
  warning(
    "No target genes had both promoter and distal mean logFC values. ",
    "Skipping promoter-vs-distal scatter plot."
  )
}

# -----------------------------------------------------------------------------
# Redistribution score plot
# -----------------------------------------------------------------------------

redistribution_df <- wide_logfc %>%
  dplyr::filter(target_class %in% c("primary_lipid", "secondary_lipid")) %>%
  dplyr::filter(!is.na(redistribution_score))

if (nrow(redistribution_df) > 0) {
  p_redist <- ggplot2::ggplot(
    redistribution_df,
    ggplot2::aes(
      x = stats::reorder(gene_name, redistribution_score),
      y = redistribution_score
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ comparison, scales = "free_y") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Target-gene distal-minus-promoter accessibility shift",
      x = NULL,
      y = "Redistribution score: distal mean logFC - promoter mean logFC"
    )

  redist_pdf <- file.path(PLOT_DIR, "target_genes_redistribution_score.pdf")
  redist_png <- file.path(PLOT_DIR, "target_genes_redistribution_score.png")

  ggplot2::ggsave(redist_pdf, p_redist, width = 8, height = 5)
  ggplot2::ggsave(redist_png, p_redist, width = 8, height = 5, dpi = 300)

  message("Wrote plot: ", redist_pdf)
  message("Wrote plot: ", redist_png)
} else {
  warning("No valid redistribution scores available. Skipping redistribution plot.")
}

# =============================================================================
# Final output checks
# =============================================================================

required_outputs <- c(
  da_comp_path,
  comp_summary_long_path,
  wide_logfc_path,
  top_target_peak_path
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]

if (length(missing_outputs) > 0) {
  stop(
    "The following expected outputs were not created:\n",
    paste(missing_outputs, collapse = "\n")
  )
}

message("============================================================")
message("03_promoter_distal_accessibility.R complete")
message("Finished: ", Sys.time())
message("============================================================")