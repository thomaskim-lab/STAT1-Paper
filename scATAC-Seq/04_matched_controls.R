# =============================================================================
# 04_matched_controls.R
# Builds matched-control genes and compares target lipid genes to controls.
#
# Matching variables available from scATAC alone:
#   - baseline WT promoter accessibility
#   - baseline WT total locus accessibility
#   - gene length
#   - number of promoter peaks
#   - number of locus peaks
#
# Optional RNA/GC can be added later.
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
message("04_matched_controls.R")
message("Started: ", Sys.time())
message("============================================================")

# =============================================================================
# Required inputs
# =============================================================================

peak_rds <- file.path(TABLE_DIR, "microglia_peakunion_peaks.granges.rds")
promoters_2kb_merged_rds <- file.path(TABLE_DIR, "all_promoters_2kb_merged.granges.rds")
genes_rds <- file.path(TABLE_DIR, "genes_mm10.granges.rds")

check_file(OBJECT_RDS, "microglia_all_final.rds")
check_file(DA_STAT1_CSV, "DA_global_STAT1_KO_vs_WT.csv")
check_file(DA_IRF1_CSV, "DA_global_IRF1_KO_vs_WT.csv")
check_file(peak_rds, "microglia_peakunion_peaks.granges.rds")
check_file(promoters_2kb_merged_rds, "all_promoters_2kb_merged.granges.rds")

# =============================================================================
# Load object and peak universe
# =============================================================================

message("Loading microglia object...")

obj <- load_microglia_object(verbose = TRUE)
DefaultAssay(obj) <- ASSAY_NAME

message("Reading peak universe: ", peak_rds)
peak_gr <- readRDS(peak_rds)
peak_gr <- clean_bed_gr(peak_gr)

if (!"peak_id" %in% colnames(mcols(peak_gr))) {
  peak_gr$peak_id <- gr_to_peak_string(peak_gr)
}

peak_names <- peak_gr$peak_id

message("Peak universe size: ", length(peak_gr))

# Get counts.
# slot = "counts" is retained for Seurat/Signac compatibility here.
counts <- Seurat::GetAssayData(
  obj,
  assay = ASSAY_NAME,
  slot = "counts"
)

if (nrow(counts) != length(peak_gr)) {
  stop(
    "Counts matrix row number does not match peak_gr length.\n",
    "nrow(counts): ", nrow(counts), "\n",
    "length(peak_gr): ", length(peak_gr)
  )
}

rownames(counts) <- peak_names

# =============================================================================
# Gene annotation and genome-wide promoter exclusion
# =============================================================================

if (file.exists(genes_rds)) {
  message("Reading gene annotation: ", genes_rds)
  genes <- readRDS(genes_rds)
  genes <- clean_bed_gr(genes)
} else {
  warning("genes_mm10.granges.rds not found. Rebuilding gene annotation from EnsDb.")
  genes <- make_gene_annotation()
  genes <- clean_bed_gr(genes)
}

message("Reading merged genome-wide promoters: ", promoters_2kb_merged_rds)
promoters_2kb_merged <- readRDS(promoters_2kb_merged_rds)
promoters_2kb_merged <- clean_bed_gr(promoters_2kb_merged)

# =============================================================================
# Construct promoter/locus/distal peak maps for all genes
# =============================================================================

message("Constructing all-gene promoter/locus/distal peak maps...")

all_gene_promoters <- GenomicRanges::promoters(
  genes,
  upstream = PROMOTER_UPSTREAM,
  downstream = PROMOTER_DOWNSTREAM
)

all_gene_promoters$gene_name <- genes$gene_name

if ("gene_id" %in% colnames(mcols(genes))) {
  all_gene_promoters$gene_id <- genes$gene_id
}

all_gene_promoters <- clean_bed_gr(all_gene_promoters)

all_gene_loci <- GenomicRanges::promoters(
  genes,
  upstream = LOCUS_WINDOW,
  downstream = LOCUS_WINDOW
)

all_gene_loci$gene_name <- genes$gene_name

if ("gene_id" %in% colnames(mcols(genes))) {
  all_gene_loci$gene_id <- genes$gene_id
}

all_gene_loci <- clean_bed_gr(all_gene_loci)

assign_gene_peaks <- function(comp_gr, label) {
  hits <- GenomicRanges::findOverlaps(
    peak_gr,
    comp_gr,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    return(data.frame())
  }

  out <- data.frame(
    peak = peak_names[queryHits(hits)],
    gene_name = comp_gr$gene_name[subjectHits(hits)],
    compartment = label,
    stringsAsFactors = FALSE
  )

  if ("gene_id" %in% colnames(mcols(comp_gr))) {
    out$gene_id <- comp_gr$gene_id[subjectHits(hits)]
  }

  dplyr::distinct(out)
}

promoter_map_all <- assign_gene_peaks(all_gene_promoters, "promoter")
locus_map_all <- assign_gene_peaks(all_gene_loci, "locus")

message("All-gene promoter peak records: ", nrow(promoter_map_all))
message("All-gene locus peak records: ", nrow(locus_map_all))

# Distal = locus-associated peaks that do not overlap any genome-wide 2 kb promoter.
peak_overlaps_any_promoter <- GenomicRanges::countOverlaps(
  peak_gr,
  promoters_2kb_merged,
  ignore.strand = TRUE
) > 0

nonpromoter_peak_names <- peak_names[!peak_overlaps_any_promoter]

all_gene_distal_map <- locus_map_all %>%
  dplyr::filter(peak %in% nonpromoter_peak_names) %>%
  dplyr::mutate(compartment = "distal") %>%
  dplyr::distinct()

all_gene_promoter_map <- promoter_map_all %>%
  dplyr::mutate(compartment = "promoter") %>%
  dplyr::distinct()

all_gene_locus_map <- locus_map_all %>%
  dplyr::mutate(compartment = "locus") %>%
  dplyr::distinct()

all_comp_map <- dplyr::bind_rows(
  all_gene_promoter_map,
  all_gene_distal_map,
  all_gene_locus_map
) %>%
  dplyr::distinct()

all_comp_map_path <- file.path(TABLE_DIR, "all_gene_peak_compartment_map.csv")
readr::write_csv(all_comp_map, all_comp_map_path)

message("Wrote all-gene peak compartment map: ", all_comp_map_path)
message("All-gene distal peak records: ", nrow(all_gene_distal_map))

# =============================================================================
# WT baseline accessibility by peak
# =============================================================================

message("Computing WT baseline accessibility by peak...")

if (!"genotype" %in% colnames(obj@meta.data)) {
  stop("obj@meta.data lacks required column: genotype")
}

wt_cells <- colnames(obj)[obj$genotype == "WT"]

if (length(wt_cells) == 0) {
  stop("No WT cells found in obj$genotype.")
}

message("WT cells used for baseline accessibility: ", length(wt_cells))

wt_counts <- counts[, wt_cells, drop = FALSE]
libsize <- Matrix::colSums(wt_counts)

if (any(libsize <= 0)) {
  stop("Some WT cells have zero total counts in the peakunion assay.")
}

# Counts per 10k per cell, then mean across WT cells.
cp10k <- Matrix::t(Matrix::t(wt_counts) / libsize) * 10000
peak_wt_mean <- Matrix::rowMeans(cp10k)

peak_metrics <- data.frame(
  peak = peak_names,
  wt_mean_accessibility = as.numeric(peak_wt_mean),
  peak_width = GenomicRanges::width(peak_gr),
  chr = as.character(GenomeInfoDb::seqnames(peak_gr)),
  start = GenomicRanges::start(peak_gr),
  end = GenomicRanges::end(peak_gr),
  stringsAsFactors = FALSE
)

peak_metrics_path <- file.path(TABLE_DIR, "peak_WT_baseline_accessibility_metrics.csv")
readr::write_csv(peak_metrics, peak_metrics_path)

message("Wrote peak baseline metrics: ", peak_metrics_path)

# =============================================================================
# Gene-level matching features
# =============================================================================

message("Building all-gene matching features...")

promoter_access <- promoter_map_all %>%
  dplyr::left_join(peak_metrics, by = "peak") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    baseline_promoter_accessibility = mean(wt_mean_accessibility, na.rm = TRUE),
    n_promoter_peaks = dplyr::n_distinct(peak),
    mean_promoter_peak_width = mean(peak_width, na.rm = TRUE),
    .groups = "drop"
  )

locus_access <- locus_map_all %>%
  dplyr::left_join(peak_metrics, by = "peak") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    baseline_locus_accessibility = mean(wt_mean_accessibility, na.rm = TRUE),
    n_locus_peaks = dplyr::n_distinct(peak),
    mean_locus_peak_width = mean(peak_width, na.rm = TRUE),
    .groups = "drop"
  )

distal_access <- all_gene_distal_map %>%
  dplyr::left_join(peak_metrics, by = "peak") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    baseline_distal_accessibility = mean(wt_mean_accessibility, na.rm = TRUE),
    n_distal_peaks = dplyr::n_distinct(peak),
    mean_distal_peak_width = mean(peak_width, na.rm = TRUE),
    .groups = "drop"
  )

gene_features <- data.frame(
  gene_name = genes$gene_name,
  chr = as.character(GenomeInfoDb::seqnames(genes)),
  gene_start = GenomicRanges::start(genes),
  gene_end = GenomicRanges::end(genes),
  gene_length = GenomicRanges::width(genes),
  strand = as.character(BiocGenerics::strand(genes)),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(promoter_access, by = "gene_name") %>%
  dplyr::left_join(locus_access, by = "gene_name") %>%
  dplyr::left_join(distal_access, by = "gene_name") %>%
  tidyr::replace_na(list(
    baseline_promoter_accessibility = 0,
    baseline_locus_accessibility = 0,
    baseline_distal_accessibility = 0,
    n_promoter_peaks = 0,
    n_locus_peaks = 0,
    n_distal_peaks = 0,
    mean_promoter_peak_width = 0,
    mean_locus_peak_width = 0,
    mean_distal_peak_width = 0
  )) %>%
  dplyr::mutate(
    is_target = gene_name %in% ALL_TARGET_GENES,
    is_primary_lipid = gene_name %in% LIPID_GENES,
    is_secondary_lipid = gene_name %in% SECONDARY_LIPID_GENES,
    target_class = dplyr::case_when(
      gene_name %in% LIPID_GENES ~ "primary_lipid",
      gene_name %in% SECONDARY_LIPID_GENES ~ "secondary_lipid",
      TRUE ~ "non_target"
    )
  )

gene_features_path <- file.path(TABLE_DIR, "all_gene_matching_features.csv")
readr::write_csv(gene_features, gene_features_path)

message("Wrote all-gene matching features: ", gene_features_path)

# =============================================================================
# Matched controls using nearest-neighbor matching in scaled feature space
# =============================================================================

message("Building matched-control genes...")

match_features <- c(
  "baseline_promoter_accessibility",
  "baseline_locus_accessibility",
  "baseline_distal_accessibility",
  "gene_length",
  "n_promoter_peaks",
  "n_locus_peaks",
  "n_distal_peaks"
)

candidate_pool <- gene_features %>%
  dplyr::filter(!gene_name %in% EXCLUDE_CONTROL_GENES) %>%
  dplyr::filter(!is_target) %>%
  dplyr::filter(n_locus_peaks > 0)

target_pool <- gene_features %>%
  dplyr::filter(gene_name %in% ALL_TARGET_GENES) %>%
  dplyr::filter(n_locus_peaks > 0)

if (nrow(target_pool) == 0) {
  stop("No target genes with n_locus_peaks > 0 were available for matching.")
}

if (nrow(candidate_pool) == 0) {
  stop("No candidate control genes with n_locus_peaks > 0 were available for matching.")
}

message("Target genes available for matching: ", nrow(target_pool))
message("Candidate control genes available for matching: ", nrow(candidate_pool))

transform_for_matching <- function(df) {
  df %>%
    dplyr::mutate(
      baseline_promoter_accessibility = log1p(baseline_promoter_accessibility),
      baseline_locus_accessibility = log1p(baseline_locus_accessibility),
      baseline_distal_accessibility = log1p(baseline_distal_accessibility),
      gene_length = log1p(gene_length),
      n_promoter_peaks = log1p(n_promoter_peaks),
      n_locus_peaks = log1p(n_locus_peaks),
      n_distal_peaks = log1p(n_distal_peaks)
    )
}

combined <- dplyr::bind_rows(
  transform_for_matching(target_pool) %>% dplyr::mutate(pool = "target"),
  transform_for_matching(candidate_pool) %>% dplyr::mutate(pool = "candidate")
)

feature_mat <- as.matrix(combined[, match_features])

# Drop zero-variance features before scaling.
feature_sd <- apply(feature_mat, 2, stats::sd, na.rm = TRUE)
keep_features <- names(feature_sd)[!is.na(feature_sd) & feature_sd > 0]

if (length(keep_features) < 2) {
  stop(
    "Fewer than two non-zero-variance matching features remain. ",
    "Cannot perform reliable nearest-neighbor matching."
  )
}

message("Matching features used: ", paste(keep_features, collapse = ", "))

scaled_mat <- scale(feature_mat[, keep_features, drop = FALSE])
rownames(scaled_mat) <- combined$gene_name

target_scaled <- scaled_mat[combined$pool == "target", , drop = FALSE]
candidate_scaled <- scaled_mat[combined$pool == "candidate", , drop = FALSE]

K_CONTROLS_PER_TARGET <- 20
set.seed(1234)

matches <- lapply(rownames(target_scaled), function(g) {
  diffs <- sweep(candidate_scaled, 2, target_scaled[g, ], FUN = "-")
  d <- sqrt(rowSums(diffs^2))
  ord <- order(d)

  n_take <- min(K_CONTROLS_PER_TARGET, length(ord))

  data.frame(
    target_gene = g,
    control_gene = rownames(candidate_scaled)[ord[seq_len(n_take)]],
    distance = d[ord[seq_len(n_take)]],
    rank = seq_len(n_take),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(
    target_pool %>% dplyr::select(target_gene = gene_name, target_class),
    by = "target_gene"
  )

matches_path <- file.path(TABLE_DIR, "matched_control_genes_per_target.csv")
readr::write_csv(matches, matches_path)

message("Wrote matched controls: ", matches_path)

# Matching diagnostics
matched_control_diagnostics <- matches %>%
  dplyr::group_by(target_gene, target_class) %>%
  dplyr::summarise(
    n_controls = dplyr::n(),
    mean_distance = mean(distance, na.rm = TRUE),
    median_distance = stats::median(distance, na.rm = TRUE),
    max_distance = max(distance, na.rm = TRUE),
    .groups = "drop"
  )

matched_control_diagnostics_path <- file.path(
  TABLE_DIR,
  "matched_control_diagnostics.csv"
)

readr::write_csv(matched_control_diagnostics, matched_control_diagnostics_path)

message("Wrote matching diagnostics: ", matched_control_diagnostics_path)

# =============================================================================
# DA summaries for all genes
# =============================================================================

message("Reading DA tables and summarizing all-gene promoter/distal/locus DA...")

da_stat1 <- read_da_table(DA_STAT1_CSV, "STAT1_KO_vs_WT")
da_irf1 <- read_da_table(DA_IRF1_CSV, "IRF1_KO_vs_WT")

da_all <- dplyr::bind_rows(da_stat1, da_irf1) %>%
  dplyr::mutate(
    logFC = as.numeric(logFC),
    FDR = as.numeric(FDR),
    PValue = as.numeric(PValue),
    logCPM = as.numeric(logCPM)
  )

# Keep all-gene compartment gene_name as the authoritative gene_name.
# The DA tables also contain gene_name from nearest-gene annotation, which would
# otherwise create gene_name.x / gene_name.y after joining.
all_comp_map_for_join <- all_comp_map %>%
  dplyr::select(
    peak,
    gene_name,
    compartment,
    dplyr::everything()
  ) %>%
  dplyr::distinct()

da_all_for_join <- da_all %>%
  dplyr::select(
    -dplyr::any_of(c("gene_name", "gene", "nearest_gene", "symbol"))
  )

da_all_comp <- da_all_for_join %>%
  dplyr::inner_join(
    all_comp_map_for_join,
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

if (nrow(da_all_comp) == 0) {
  stop("No DA peaks overlapped all-gene promoter/distal/locus compartments.")
}

da_all_comp_path <- file.path(TABLE_DIR, "all_gene_DA_peaks_in_compartments.csv")
readr::write_csv(da_all_comp, da_all_comp_path)

message("Wrote all-gene DA-compartment table: ", da_all_comp_path)

all_gene_da_summary_long <- da_all_comp %>%
  dplyr::group_by(gene_name, comparison, compartment) %>%
  dplyr::summarise(
    n_peaks = dplyr::n_distinct(peak),
    mean_logFC = safe_mean(logFC),
    median_logFC = safe_median(logFC),
    mean_logCPM = safe_mean(logCPM),

    n_DA_lost_FDR10 = sum(sig_lost(logFC, FDR, FDR_CUTOFF), na.rm = TRUE),
    n_DA_gained_FDR10 = sum(sig_gained(logFC, FDR, FDR_CUTOFF), na.rm = TRUE),

    n_DA_lost_FDR05 = sum(sig_lost(logFC, FDR, STRICT_FDR_CUTOFF), na.rm = TRUE),
    n_DA_gained_FDR05 = sum(sig_gained(logFC, FDR, STRICT_FDR_CUTOFF), na.rm = TRUE),

    min_FDR = ifelse(all(is.na(FDR)), NA_real_, min(FDR, na.rm = TRUE)),
    .groups = "drop"
  )

all_gene_da_summary_long_path <- file.path(
  TABLE_DIR,
  "all_gene_promoter_distal_locus_DA_summary_long.csv"
)

readr::write_csv(all_gene_da_summary_long, all_gene_da_summary_long_path)

all_gene_da_wide <- all_gene_da_summary_long %>%
  dplyr::select(
    gene_name,
    comparison,
    compartment,
    mean_logFC,
    median_logFC,
    mean_logCPM,
    n_peaks,
    n_DA_lost_FDR10,
    n_DA_gained_FDR10,
    n_DA_lost_FDR05,
    n_DA_gained_FDR05,
    min_FDR
  ) %>%
  tidyr::pivot_wider(
    names_from = compartment,
    values_from = c(
      mean_logFC,
      median_logFC,
      mean_logCPM,
      n_peaks,
      n_DA_lost_FDR10,
      n_DA_gained_FDR10,
      n_DA_lost_FDR05,
      n_DA_gained_FDR05,
      min_FDR
    )
  )

ensure_column <- function(df, col, default = NA_real_) {
  if (!col %in% colnames(df)) {
    df[[col]] <- default
  }
  df
}

expected_cols <- c(
  "mean_logFC_promoter",
  "mean_logFC_distal",
  "mean_logFC_locus",
  "median_logFC_promoter",
  "median_logFC_distal",
  "median_logFC_locus",
  "mean_logCPM_promoter",
  "mean_logCPM_distal",
  "mean_logCPM_locus",
  "n_peaks_promoter",
  "n_peaks_distal",
  "n_peaks_locus",
  "n_DA_lost_FDR10_promoter",
  "n_DA_lost_FDR10_distal",
  "n_DA_lost_FDR10_locus",
  "n_DA_gained_FDR10_promoter",
  "n_DA_gained_FDR10_distal",
  "n_DA_gained_FDR10_locus",
  "n_DA_lost_FDR05_promoter",
  "n_DA_lost_FDR05_distal",
  "n_DA_lost_FDR05_locus",
  "n_DA_gained_FDR05_promoter",
  "n_DA_gained_FDR05_distal",
  "n_DA_gained_FDR05_locus",
  "min_FDR_promoter",
  "min_FDR_distal",
  "min_FDR_locus"
)

for (col in expected_cols) {
  all_gene_da_wide <- ensure_column(all_gene_da_wide, col, NA_real_)
}

all_gene_da_wide <- all_gene_da_wide %>%
  dplyr::mutate(
    redistribution_score = dplyr::if_else(
      !is.na(mean_logFC_promoter) &
        !is.na(mean_logFC_distal) &
        !is.na(n_peaks_promoter) &
        !is.na(n_peaks_distal) &
        n_peaks_promoter > 0 &
        n_peaks_distal > 0,
      mean_logFC_distal - mean_logFC_promoter,
      NA_real_
    )
  ) %>%
  dplyr::left_join(
    gene_features %>%
      dplyr::select(
        gene_name,
        target_class,
        is_target,
        is_primary_lipid,
        is_secondary_lipid,
        baseline_promoter_accessibility,
        baseline_locus_accessibility,
        baseline_distal_accessibility,
        gene_length,
        n_promoter_peaks,
        n_locus_peaks,
        n_distal_peaks
      ),
    by = "gene_name"
  )

all_gene_da_wide_path <- file.path(
  TABLE_DIR,
  "all_gene_promoter_distal_locus_DA_summary.csv"
)

readr::write_csv(all_gene_da_wide, all_gene_da_wide_path)

message("Wrote all-gene DA summary: ", all_gene_da_wide_path)

# =============================================================================
# Matched comparison table
# =============================================================================

message("Building target-vs-matched-control DA value table...")

# all_gene_da_wide contains gene-level annotations such as target_class.
# For this matched-control table, target_class should come from the original
# target gene, not from the matched control gene. Remove these columns before
# joining to avoid target_class.x / target_class.y collisions.
da_wide_for_matching <- all_gene_da_wide %>%
  dplyr::select(
    -dplyr::any_of(c(
      "target_class",
      "is_target",
      "is_primary_lipid",
      "is_secondary_lipid"
    ))
  )

matched_values <- matches %>%
  dplyr::left_join(
    da_wide_for_matching,
    by = c("control_gene" = "gene_name"),
    relationship = "many-to-many"
  ) %>%
  dplyr::rename(matched_gene = control_gene) %>%
  dplyr::mutate(group = "matched_control")

target_values <- target_pool %>%
  dplyr::select(
    target_gene = gene_name,
    target_class
  ) %>%
  dplyr::left_join(
    da_wide_for_matching,
    by = c("target_gene" = "gene_name"),
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    matched_gene = target_gene,
    group = "target",
    distance = 0,
    rank = 0
  )

matched_test_input <- dplyr::bind_rows(
  target_values,
  matched_values
) %>%
  dplyr::arrange(target_class, target_gene, comparison, group, rank)

matched_test_input_path <- file.path(
  TABLE_DIR,
  "target_vs_matched_control_DA_values.csv"
)

readr::write_csv(matched_test_input, matched_test_input_path)

message("Wrote target-vs-matched-control values: ", matched_test_input_path)

# =============================================================================
# One-sided tests
# =============================================================================

message("Running one-sided matched-control tests...")

run_one_sided_tests <- function(df, metric, alternative, target_class_filter = NULL) {
  test_df <- df

  if (!is.null(target_class_filter)) {
    test_df <- test_df %>%
      dplyr::filter(target_class %in% target_class_filter)
  }

  if (!metric %in% colnames(test_df)) {
    warning("Metric not present, skipping: ", metric)
    return(data.frame())
  }

  test_df %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_target = sum(group == "target" & !is.na(.data[[metric]])),
      n_control = sum(group == "matched_control" & !is.na(.data[[metric]])),

      target_mean = mean(.data[[metric]][group == "target"], na.rm = TRUE),
      control_mean = mean(.data[[metric]][group == "matched_control"], na.rm = TRUE),

      target_median = stats::median(.data[[metric]][group == "target"], na.rm = TRUE),
      control_median = stats::median(.data[[metric]][group == "matched_control"], na.rm = TRUE),

      p_value = {
        x <- .data[[metric]][group == "target"]
        y <- .data[[metric]][group == "matched_control"]

        x <- x[!is.na(x)]
        y <- y[!is.na(y)]

        if (length(x) >= 3 && length(y) >= 3) {
          suppressWarnings(
            stats::wilcox.test(
              x,
              y,
              alternative = alternative,
              exact = FALSE
            )$p.value
          )
        } else {
          NA_real_
        }
      },

      metric = metric,
      alternative = alternative,
      target_class_tested = ifelse(
        is.null(target_class_filter),
        "all_targets",
        paste(target_class_filter, collapse = "+")
      ),

      .groups = "drop"
    )
}

test_specs <- dplyr::tribble(
  ~metric, ~alternative, ~interpretation,
  "mean_logFC_promoter", "less", "target promoter accessibility logFC is lower than matched controls",
  "mean_logFC_distal", "greater", "target distal accessibility logFC is higher than matched controls",
  "mean_logFC_locus", "less", "target locus accessibility logFC is lower than matched controls",
  "redistribution_score", "greater", "target distal-minus-promoter score is higher than matched controls",
  "n_DA_lost_FDR10_promoter", "greater", "target promoters have more FDR10 lost peaks than matched controls",
  "n_DA_gained_FDR10_distal", "greater", "target distal compartments have more FDR10 gained peaks than matched controls",
  "n_DA_lost_FDR05_promoter", "greater", "target promoters have more FDR05 lost peaks than matched controls",
  "n_DA_gained_FDR05_distal", "greater", "target distal compartments have more FDR05 gained peaks than matched controls"
)

tests_all <- lapply(seq_len(nrow(test_specs)), function(i) {
  run_one_sided_tests(
    matched_test_input,
    metric = test_specs$metric[i],
    alternative = test_specs$alternative[i],
    target_class_filter = c("primary_lipid", "secondary_lipid")
  ) %>%
    dplyr::mutate(interpretation = test_specs$interpretation[i])
}) %>%
  dplyr::bind_rows()

tests_primary <- lapply(seq_len(nrow(test_specs)), function(i) {
  run_one_sided_tests(
    matched_test_input,
    metric = test_specs$metric[i],
    alternative = test_specs$alternative[i],
    target_class_filter = "primary_lipid"
  ) %>%
    dplyr::mutate(interpretation = test_specs$interpretation[i])
}) %>%
  dplyr::bind_rows()

tests <- dplyr::bind_rows(tests_all, tests_primary) %>%
  dplyr::group_by(comparison, target_class_tested) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(comparison, target_class_tested, p_value)

tests_path <- file.path(TABLE_DIR, "matched_control_one_sided_tests.csv")
readr::write_csv(tests, tests_path)

message("Wrote matched-control one-sided tests: ", tests_path)

# =============================================================================
# Simple plots
# =============================================================================

message("Generating matched-control diagnostic plots...")

plot_df <- matched_test_input %>%
  dplyr::filter(
    target_class %in% c("primary_lipid", "secondary_lipid"),
    !is.na(comparison)
  )

if (nrow(plot_df) > 0) {
  p_redist <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = group,
      y = redistribution_score
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1) +
    ggplot2::facet_grid(target_class ~ comparison) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Target lipid-handling genes versus matched controls",
      subtitle = "Redistribution score = distal mean logFC - promoter mean logFC",
      x = NULL,
      y = "Redistribution score"
    )

  redist_pdf <- file.path(PLOT_DIR, "matched_controls_redistribution_score.pdf")
  redist_png <- file.path(PLOT_DIR, "matched_controls_redistribution_score.png")

  ggplot2::ggsave(redist_pdf, p_redist, width = 8, height = 5)
  ggplot2::ggsave(redist_png, p_redist, width = 8, height = 5, dpi = 300)

  message("Wrote plot: ", redist_pdf)
  message("Wrote plot: ", redist_png)

  p_prom_dist <- plot_df %>%
    dplyr::select(
      comparison,
      target_class,
      group,
      mean_logFC_promoter,
      mean_logFC_distal,
      mean_logFC_locus
    ) %>%
    tidyr::pivot_longer(
      cols = c(mean_logFC_promoter, mean_logFC_distal, mean_logFC_locus),
      names_to = "metric",
      values_to = "value"
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = group,
        y = value
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1) +
    ggplot2::facet_grid(metric ~ comparison, scales = "free_y") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Promoter, distal, and locus logFC in targets versus matched controls",
      x = NULL,
      y = "Mean logFC"
    )

  prom_dist_pdf <- file.path(PLOT_DIR, "matched_controls_promoter_distal_locus_logFC.pdf")
  prom_dist_png <- file.path(PLOT_DIR, "matched_controls_promoter_distal_locus_logFC.png")

  ggplot2::ggsave(prom_dist_pdf, p_prom_dist, width = 9, height = 7)
  ggplot2::ggsave(prom_dist_png, p_prom_dist, width = 9, height = 7, dpi = 300)

  message("Wrote plot: ", prom_dist_pdf)
  message("Wrote plot: ", prom_dist_png)
} else {
  warning("No matched-control plot data available.")
}

# =============================================================================
# Final output checks
# =============================================================================

required_outputs <- c(
  all_comp_map_path,
  peak_metrics_path,
  gene_features_path,
  matches_path,
  matched_control_diagnostics_path,
  da_all_comp_path,
  all_gene_da_summary_long_path,
  all_gene_da_wide_path,
  matched_test_input_path,
  tests_path
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]

if (length(missing_outputs) > 0) {
  stop(
    "The following expected outputs were not created:\n",
    paste(missing_outputs, collapse = "\n")
  )
}

message("============================================================")
message("04_matched_controls.R complete")
message("Finished: ", Sys.time())
message("============================================================")