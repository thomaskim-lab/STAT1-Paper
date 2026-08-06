# =============================================================================
# 02_define_compartments.R
# Defines target-gene promoter, distal, and whole-locus compartments.
#
# Distal peaks exclude all promoters genome-wide, not only target-gene promoters.
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
message("02_define_compartments.R")
message("Started: ", Sys.time())
message("============================================================")

# =============================================================================
# Required inputs from 01_export_peakunion_and_promoters.R
# =============================================================================

peak_rds <- file.path(TABLE_DIR, "microglia_peakunion_peaks.granges.rds")
promoters_2kb_merged_rds <- file.path(TABLE_DIR, "all_promoters_2kb_merged.granges.rds")
promoters_1kb_merged_rds <- file.path(TABLE_DIR, "all_promoters_1kb_merged.granges.rds")
genes_rds <- file.path(TABLE_DIR, "genes_mm10.granges.rds")

check_file(peak_rds, "microglia_peakunion_peaks.granges.rds")
check_file(promoters_2kb_merged_rds, "all_promoters_2kb_merged.granges.rds")
check_file(promoters_1kb_merged_rds, "all_promoters_1kb_merged.granges.rds")

message("Reading peak universe: ", peak_rds)
peak_gr <- readRDS(peak_rds)
peak_gr <- clean_bed_gr(peak_gr)

if (!"peak_id" %in% colnames(mcols(peak_gr))) {
  peak_gr$peak_id <- gr_to_peak_string(peak_gr)
}

message("Peak universe size: ", length(peak_gr))

message("Reading genome-wide merged promoters:")
message("  ", promoters_2kb_merged_rds)
message("  ", promoters_1kb_merged_rds)

promoters_2kb_merged <- readRDS(promoters_2kb_merged_rds)
promoters_1kb_merged <- readRDS(promoters_1kb_merged_rds)

promoters_2kb_merged <- clean_bed_gr(promoters_2kb_merged)
promoters_1kb_merged <- clean_bed_gr(promoters_1kb_merged)

if (file.exists(genes_rds)) {
  message("Reading gene annotation: ", genes_rds)
  genes <- readRDS(genes_rds)
  genes <- clean_bed_gr(genes)
} else {
  warning("genes_mm10.granges.rds not found. Rebuilding gene annotation from EnsDb.")
  genes <- make_gene_annotation()
  genes <- clean_bed_gr(genes)
}

# =============================================================================
# Select target genes
# =============================================================================

target_genes <- genes[genes$gene_name %in% ALL_TARGET_GENES]
target_genes <- sort(target_genes)

missing_genes <- setdiff(ALL_TARGET_GENES, target_genes$gene_name)

if (length(missing_genes) > 0) {
  warning(
    "Missing target genes in annotation: ",
    paste(missing_genes, collapse = ", ")
  )
}

if (length(target_genes) == 0) {
  stop("None of the requested target genes were found in the annotation.")
}

target_gene_table <- data.frame(
  gene_name = target_genes$gene_name,
  chr = as.character(seqnames(target_genes)),
  start = start(target_genes),
  end = end(target_genes),
  width = width(target_genes),
  strand = as.character(strand(target_genes)),
  stringsAsFactors = FALSE
)

if ("gene_id" %in% colnames(mcols(target_genes))) {
  target_gene_table$gene_id <- target_genes$gene_id
}

target_gene_table_path <- file.path(TABLE_DIR, "target_genes_found_in_annotation.csv")
readr::write_csv(target_gene_table, target_gene_table_path)

saveRDS(target_genes, file.path(TABLE_DIR, "target_genes_found_in_annotation.granges.rds"))

message("Target genes found: ", length(target_genes))
message("Wrote target gene table: ", target_gene_table_path)

# =============================================================================
# Helper functions
# =============================================================================

add_gene_metadata <- function(gr, genes_gr, compartment_label) {
  gr$gene_name <- genes_gr$gene_name
  gr$compartment <- compartment_label

  if ("gene_id" %in% colnames(mcols(genes_gr))) {
    gr$gene_id <- genes_gr$gene_id
  }

  gr
}

make_gene_body_locus <- function(genes_gr, flank = LOCUS_WINDOW) {
  out <- genes_gr

  start(out) <- pmax(1L, start(out) - flank)
  end(out) <- end(out) + flank

  out$gene_name <- genes_gr$gene_name
  out$compartment <- "locus_gene_body_100kb"

  if ("gene_id" %in% colnames(mcols(genes_gr))) {
    out$gene_id <- genes_gr$gene_id
  }

  clean_bed_gr(out)
}

make_distal_for_gene <- function(locus_gr, promoter_exclusion_gr) {
  # Distal interval = target gene TSS ±100 kb minus all genome-wide promoters.
  # This intentionally removes promoters for every annotated gene, not only the
  # target gene promoter.
  out <- GenomicRanges::setdiff(
    locus_gr,
    promoter_exclusion_gr,
    ignore.strand = TRUE
  )

  if (length(out) > 0) {
    out$gene_name <- locus_gr$gene_name
    out$compartment <- "distal_TSS_100kb_noPromoter"

    if ("gene_id" %in% colnames(mcols(locus_gr))) {
      out$gene_id <- locus_gr$gene_id
    }
  }

  out
}

assign_peaks_to_compartment <- function(peaks, compartments, label) {
  if (length(compartments) == 0) {
    return(data.frame())
  }

  hits <- findOverlaps(peaks, compartments, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(data.frame())
  }

  peak_hits <- peaks[queryHits(hits)]
  comp_hits <- compartments[subjectHits(hits)]

  overlap_gr <- suppressWarnings(
    pintersect(
      peak_hits,
      comp_hits,
      ignore.strand = TRUE
    )
  )

  out <- data.frame(
    peak = gr_to_peak_string(peak_hits),
    chr = as.character(seqnames(peak_hits)),
    start = start(peak_hits),
    end = end(peak_hits),
    width = width(peak_hits),
    gene_name = comp_hits$gene_name,
    compartment = label,
    overlap_width = width(overlap_gr),
    stringsAsFactors = FALSE
  )

  if ("gene_id" %in% colnames(mcols(comp_hits))) {
    out$gene_id <- comp_hits$gene_id
  }

  dplyr::distinct(out)
}

write_compartment <- function(gr, bed_name, rds_name) {
  gr <- clean_bed_gr(gr)

  bed_path <- file.path(BED_DIR, bed_name)
  rds_path <- file.path(TABLE_DIR, rds_name)

  rtracklayer::export(gr, bed_path, format = "BED")
  saveRDS(gr, rds_path)

  message("Wrote BED: ", bed_path, " [", length(gr), " ranges]")
  message("Wrote RDS: ", rds_path)

  invisible(list(bed = bed_path, rds = rds_path))
}

# =============================================================================
# Define target-gene compartments as genomic intervals
# =============================================================================

message("Defining target-gene compartments...")

# -----------------------------------------------------------------------------
# Promoter: TSS ±2 kb
# -----------------------------------------------------------------------------

target_promoters_2kb <- promoters(
  target_genes,
  upstream = PROMOTER_UPSTREAM,
  downstream = PROMOTER_DOWNSTREAM
)

target_promoters_2kb <- add_gene_metadata(
  target_promoters_2kb,
  target_genes,
  "promoter_2kb"
)

target_promoters_2kb <- clean_bed_gr(target_promoters_2kb)

# -----------------------------------------------------------------------------
# Sensitivity promoter: TSS ±1 kb
# -----------------------------------------------------------------------------

target_promoters_1kb <- promoters(
  target_genes,
  upstream = PROMOTER_SENS_UPSTREAM,
  downstream = PROMOTER_SENS_DOWNSTREAM
)

target_promoters_1kb <- add_gene_metadata(
  target_promoters_1kb,
  target_genes,
  "promoter_1kb"
)

target_promoters_1kb <- clean_bed_gr(target_promoters_1kb)

# -----------------------------------------------------------------------------
# Whole locus: TSS ±100 kb
# -----------------------------------------------------------------------------

target_loci_TSS_100kb <- promoters(
  target_genes,
  upstream = LOCUS_WINDOW,
  downstream = LOCUS_WINDOW
)

target_loci_TSS_100kb <- add_gene_metadata(
  target_loci_TSS_100kb,
  target_genes,
  "locus_TSS_100kb"
)

target_loci_TSS_100kb <- clean_bed_gr(target_loci_TSS_100kb)

# -----------------------------------------------------------------------------
# Alternative whole locus: gene body ±100 kb
# -----------------------------------------------------------------------------

target_loci_gene_body_100kb <- make_gene_body_locus(
  target_genes,
  flank = LOCUS_WINDOW
)

# -----------------------------------------------------------------------------
# Distal interval: TSS ±100 kb minus all genome-wide 2 kb promoters
# -----------------------------------------------------------------------------

target_distal_list <- GenomicRanges::GRangesList(
  lapply(seq_along(target_loci_TSS_100kb), function(i) {
    make_distal_for_gene(
      target_loci_TSS_100kb[i],
      promoters_2kb_merged
    )
  })
)

target_distal_100kb <- unlist(target_distal_list, use.names = FALSE)
target_distal_100kb <- clean_bed_gr(target_distal_100kb)

message("Target promoter 2 kb intervals: ", length(target_promoters_2kb))
message("Target promoter 1 kb intervals: ", length(target_promoters_1kb))
message("Target TSS ±100 kb loci: ", length(target_loci_TSS_100kb))
message("Target gene-body ±100 kb loci: ", length(target_loci_gene_body_100kb))
message("Target distal intervals after promoter exclusion: ", length(target_distal_100kb))

# =============================================================================
# Validate distal exclusion
# =============================================================================

distal_promoter_hits <- countOverlaps(
  target_distal_100kb,
  promoters_2kb_merged,
  ignore.strand = TRUE
)

if (any(distal_promoter_hits > 0)) {
  stop(
    "Distal intervals still overlap genome-wide 2 kb promoters. ",
    "Check promoter exclusion logic."
  )
}

message("Validated: distal intervals do not overlap genome-wide 2 kb promoters.")

# =============================================================================
# Export compartment BEDs and RDS files
# =============================================================================

write_compartment(
  target_promoters_2kb,
  "target_gene_promoters_2kb.bed",
  "target_gene_promoters_2kb.granges.rds"
)

write_compartment(
  target_promoters_1kb,
  "target_gene_promoters_1kb.bed",
  "target_gene_promoters_1kb.granges.rds"
)

write_compartment(
  target_loci_TSS_100kb,
  "target_gene_loci_TSS_100kb.bed",
  "target_gene_loci_TSS_100kb.granges.rds"
)

write_compartment(
  target_loci_gene_body_100kb,
  "target_gene_loci_gene_body_100kb.bed",
  "target_gene_loci_gene_body_100kb.granges.rds"
)

write_compartment(
  target_distal_100kb,
  "target_gene_distal_TSS_100kb_noPromoter.bed",
  "target_gene_distal_TSS_100kb_noPromoter.granges.rds"
)

# =============================================================================
# Assign peakunion peaks to compartments
# =============================================================================

message("Assigning peakunion peaks to target-gene compartments...")

peak_promoter_2kb <- assign_peaks_to_compartment(
  peak_gr,
  target_promoters_2kb,
  "promoter_2kb"
)

peak_promoter_1kb <- assign_peaks_to_compartment(
  peak_gr,
  target_promoters_1kb,
  "promoter_1kb"
)

peak_locus_TSS_100kb <- assign_peaks_to_compartment(
  peak_gr,
  target_loci_TSS_100kb,
  "locus_TSS_100kb"
)

peak_locus_gene_body_100kb <- assign_peaks_to_compartment(
  peak_gr,
  target_loci_gene_body_100kb,
  "locus_gene_body_100kb"
)

peak_distal_100kb <- assign_peaks_to_compartment(
  peak_gr,
  target_distal_100kb,
  "distal_TSS_100kb_noPromoter"
)

compartment_map <- dplyr::bind_rows(
  peak_promoter_2kb,
  peak_promoter_1kb,
  peak_locus_TSS_100kb,
  peak_locus_gene_body_100kb,
  peak_distal_100kb
)

compartment_map <- compartment_map %>%
  dplyr::arrange(gene_name, compartment, chr, start, end)

compartment_map_csv <- file.path(TABLE_DIR, "peak_to_target_gene_compartments.csv")
compartment_map_rds <- file.path(TABLE_DIR, "peak_to_target_gene_compartments.rds")

readr::write_csv(compartment_map, compartment_map_csv)
saveRDS(compartment_map, compartment_map_rds)

message("Wrote peak-to-compartment map CSV: ", compartment_map_csv)
message("Wrote peak-to-compartment map RDS: ", compartment_map_rds)

# =============================================================================
# Gene-level compartment peak count summaries
# =============================================================================

message("Writing gene-level peak-count summaries...")

compartment_peak_counts <- compartment_map %>%
  dplyr::group_by(gene_name, compartment) %>%
  dplyr::summarise(
    n_peak_gene_compartment_records = dplyr::n(),
    n_unique_peaks = dplyr::n_distinct(peak),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = compartment,
    values_from = c(n_peak_gene_compartment_records, n_unique_peaks),
    values_fill = 0
  )

target_gene_peak_count_path <- file.path(
  TABLE_DIR,
  "target_gene_compartment_peak_counts.csv"
)

readr::write_csv(compartment_peak_counts, target_gene_peak_count_path)

message("Wrote gene-level peak-count summary: ", target_gene_peak_count_path)

# Also write simple global compartment totals
global_compartment_counts <- compartment_map %>%
  dplyr::group_by(compartment) %>%
  dplyr::summarise(
    n_gene_peak_records = dplyr::n(),
    n_unique_peaks = dplyr::n_distinct(peak),
    n_genes = dplyr::n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  dplyr::arrange(compartment)

global_compartment_counts_path <- file.path(
  TABLE_DIR,
  "global_target_compartment_peak_counts.csv"
)

readr::write_csv(global_compartment_counts, global_compartment_counts_path)

message("Wrote global compartment count summary: ", global_compartment_counts_path)

print(global_compartment_counts)

# =============================================================================
# Final output checks
# =============================================================================

required_outputs <- c(
  file.path(BED_DIR, "target_gene_promoters_2kb.bed"),
  file.path(BED_DIR, "target_gene_promoters_1kb.bed"),
  file.path(BED_DIR, "target_gene_loci_TSS_100kb.bed"),
  file.path(BED_DIR, "target_gene_loci_gene_body_100kb.bed"),
  file.path(BED_DIR, "target_gene_distal_TSS_100kb_noPromoter.bed"),
  file.path(TABLE_DIR, "target_gene_promoters_2kb.granges.rds"),
  file.path(TABLE_DIR, "target_gene_promoters_1kb.granges.rds"),
  file.path(TABLE_DIR, "target_gene_loci_TSS_100kb.granges.rds"),
  file.path(TABLE_DIR, "target_gene_loci_gene_body_100kb.granges.rds"),
  file.path(TABLE_DIR, "target_gene_distal_TSS_100kb_noPromoter.granges.rds"),
  compartment_map_csv,
  compartment_map_rds,
  target_gene_peak_count_path,
  global_compartment_counts_path
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]

if (length(missing_outputs) > 0) {
  stop(
    "The following expected outputs were not created:\n",
    paste(missing_outputs, collapse = "\n")
  )
}

message("============================================================")
message("02_define_compartments.R complete")
message("Finished: ", Sys.time())
message("============================================================")