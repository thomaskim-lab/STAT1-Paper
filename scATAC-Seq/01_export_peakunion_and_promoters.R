# =============================================================================
# 01_export_peakunion_and_promoters.R
# Exports the peakunion universe and genome-wide promoter BED files.
# GenomeDK / Hypothalamus version.
#
# Run first.
# =============================================================================

# -----------------------------------------------------------------------------
# Robust source of 00_config.R
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", args[grep(file_arg, args)]))
script_dir <- dirname(script_path)

source(file.path(script_dir, "00_config.R"))

check_inputs(require_tools = FALSE)

message("============================================================")
message("01_export_peakunion_and_promoters.R")
message("Started: ", Sys.time())
message("============================================================")

# =============================================================================
# Load object
# =============================================================================

obj <- load_microglia_object(verbose = TRUE)

DefaultAssay(obj) <- ASSAY_NAME

message("Active assay: ", DefaultAssay(obj))
message("Assays: ", paste(Assays(obj), collapse = ", "))
message("Object dimensions: ", paste(dim(obj), collapse = " x "))

# Save basic sample summary for QC / manuscript tracking
sample_summary <- as.data.frame.matrix(table(obj$genotype, obj$mouse_id))
sample_summary_path <- file.path(TABLE_DIR, "sample_genotype_by_mouse_cell_counts.csv")
readr::write_csv(
  tibble::rownames_to_column(sample_summary, var = "genotype"),
  sample_summary_path
)
message("Wrote sample summary: ", sample_summary_path)

# =============================================================================
# Export peakunion peak universe
# =============================================================================

message("Exporting peakunion peak universe...")

peak_gr <- granges(obj[[ASSAY_NAME]])
peak_gr <- clean_bed_gr(peak_gr)

peak_gr$peak_id <- gr_to_peak_string(peak_gr)

peak_bed <- file.path(BED_DIR, "microglia_peakunion_peaks.bed")
peak_rds <- file.path(TABLE_DIR, "microglia_peakunion_peaks.granges.rds")

rtracklayer::export(peak_gr, peak_bed, format = "BED")
saveRDS(peak_gr, peak_rds)

message("Exported peak universe BED: ", peak_bed)
message("Exported peak universe RDS: ", peak_rds)
message("Number of peakunion peaks: ", length(peak_gr))

# Also save a simple peak metadata table
peak_table <- data.frame(
  peak = gr_to_peak_string(peak_gr),
  chr = as.character(seqnames(peak_gr)),
  start = start(peak_gr),
  end = end(peak_gr),
  width = width(peak_gr),
  stringsAsFactors = FALSE
)

peak_table_path <- file.path(TABLE_DIR, "microglia_peakunion_peaks_metadata.csv")
readr::write_csv(peak_table, peak_table_path)
message("Wrote peak metadata table: ", peak_table_path)

# =============================================================================
# Genome-wide promoters
# =============================================================================

message("Making genome-wide gene annotation and promoter BEDs...")

genes <- make_gene_annotation()
genes <- clean_bed_gr(genes)

genes$gene_id <- if ("gene_id" %in% colnames(mcols(genes))) genes$gene_id else NA_character_
genes$gene_name <- genes$gene_name

genes_rds <- file.path(TABLE_DIR, "genes_mm10.granges.rds")
genes_bed <- file.path(BED_DIR, "genes_mm10_gene_annotated.bed")

saveRDS(genes, genes_rds)
rtracklayer::export(genes, genes_bed, format = "BED")

message("Saved gene annotation RDS: ", genes_rds)
message("Saved gene annotation BED: ", genes_bed)
message("Number of genes: ", length(genes))

# -----------------------------------------------------------------------------
# TSS ±2 kb promoters
# -----------------------------------------------------------------------------

promoters_2kb <- promoters(
  genes,
  upstream = PROMOTER_UPSTREAM,
  downstream = PROMOTER_DOWNSTREAM
)

promoters_2kb$gene_name <- genes$gene_name

if ("gene_id" %in% colnames(mcols(genes))) {
  promoters_2kb$gene_id <- genes$gene_id
}

promoters_2kb <- clean_bed_gr(promoters_2kb)

# -----------------------------------------------------------------------------
# TSS ±1 kb sensitivity promoters
# -----------------------------------------------------------------------------

promoters_1kb <- promoters(
  genes,
  upstream = PROMOTER_SENS_UPSTREAM,
  downstream = PROMOTER_SENS_DOWNSTREAM
)

promoters_1kb$gene_name <- genes$gene_name

if ("gene_id" %in% colnames(mcols(genes))) {
  promoters_1kb$gene_id <- genes$gene_id
}

promoters_1kb <- clean_bed_gr(promoters_1kb)

# -----------------------------------------------------------------------------
# Merged promoters for genome-wide distal exclusion
# -----------------------------------------------------------------------------

promoters_2kb_merged <- reduce(promoters_2kb, ignore.strand = TRUE)
promoters_1kb_merged <- reduce(promoters_1kb, ignore.strand = TRUE)

promoters_2kb_merged <- clean_bed_gr(promoters_2kb_merged)
promoters_1kb_merged <- clean_bed_gr(promoters_1kb_merged)

# =============================================================================
# Export promoters
# =============================================================================

prom_2kb_annot_bed <- file.path(BED_DIR, "all_promoters_2kb_gene_annotated.bed")
prom_1kb_annot_bed <- file.path(BED_DIR, "all_promoters_1kb_gene_annotated.bed")
prom_2kb_merged_bed <- file.path(BED_DIR, "all_promoters_2kb_merged.bed")
prom_1kb_merged_bed <- file.path(BED_DIR, "all_promoters_1kb_merged.bed")

rtracklayer::export(promoters_2kb, prom_2kb_annot_bed, format = "BED")
rtracklayer::export(promoters_1kb, prom_1kb_annot_bed, format = "BED")
rtracklayer::export(promoters_2kb_merged, prom_2kb_merged_bed, format = "BED")
rtracklayer::export(promoters_1kb_merged, prom_1kb_merged_bed, format = "BED")

saveRDS(
  promoters_2kb,
  file.path(TABLE_DIR, "all_promoters_2kb_gene_annotated.granges.rds")
)

saveRDS(
  promoters_1kb,
  file.path(TABLE_DIR, "all_promoters_1kb_gene_annotated.granges.rds")
)

saveRDS(
  promoters_2kb_merged,
  file.path(TABLE_DIR, "all_promoters_2kb_merged.granges.rds")
)

saveRDS(
  promoters_1kb_merged,
  file.path(TABLE_DIR, "all_promoters_1kb_merged.granges.rds")
)

message("Exported promoter BED files:")
message("  ", prom_2kb_annot_bed)
message("  ", prom_1kb_annot_bed)
message("  ", prom_2kb_merged_bed)
message("  ", prom_1kb_merged_bed)

message("Number of annotated 2 kb promoters: ", length(promoters_2kb))
message("Number of merged 2 kb promoter intervals: ", length(promoters_2kb_merged))
message("Number of annotated 1 kb promoters: ", length(promoters_1kb))
message("Number of merged 1 kb promoter intervals: ", length(promoters_1kb_merged))

# =============================================================================
# Final checks
# =============================================================================

required_outputs <- c(
  peak_bed,
  peak_rds,
  genes_rds,
  genes_bed,
  prom_2kb_annot_bed,
  prom_1kb_annot_bed,
  prom_2kb_merged_bed,
  prom_1kb_merged_bed
)

missing_outputs <- required_outputs[!file.exists(required_outputs)]

if (length(missing_outputs) > 0) {
  stop(
    "The following expected outputs were not created:\n",
    paste(missing_outputs, collapse = "\n")
  )
}

message("============================================================")
message("01_export_peakunion_and_promoters.R complete")
message("Finished: ", Sys.time())
message("============================================================")