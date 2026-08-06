# =============================================================================
# 00_config.R
# Shared paths and parameters for scATAC Figure 2 promoter/distal/enhancer analysis
# GenomeDK / Hypothalamus project version
# =============================================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(EnsDb.Mmusculus.v79)
  library(ensembldb)
  library(rtracklayer)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# =============================================================================
# SERVER PATHS
# =============================================================================

PROJECT <- "/faststorage/project/Hypothalamus"
DATA_DIR <- file.path(PROJECT, "data")
SCRIPT_DIR <- file.path(PROJECT, "script")
TOOLS_DIR <- file.path(PROJECT, "tools")
UCSC_DIR <- file.path(TOOLS_DIR, "ucsc")

SCATAC_DIR <- file.path(
  DATA_DIR,
  "scATAC/merged_microglia_STAT1Project"
)

OBJECT_RDS <- file.path(
  SCATAC_DIR,
  "microglia_all_final.rds"
)

# Existing pseudobulk DA output from your edgeR script
PEAK_DA_DIR <- file.path(SCATAC_DIR, "peak")

DA_STAT1_CSV <- file.path(
  PEAK_DA_DIR,
  "DA_global_STAT1_KO_vs_WT.csv"
)

DA_IRF1_CSV <- file.path(
  PEAK_DA_DIR,
  "DA_global_IRF1_KO_vs_WT.csv"
)

DA_WITHIN_CLUSTER_CSV <- file.path(
  PEAK_DA_DIR,
  "DA_within_cluster_all.csv"
)

# New Figure 2 scATAC output directory
FIG2_DIR <- file.path(
  SCATAC_DIR,
  "Fig2_scATAC"
)

BED_DIR <- file.path(FIG2_DIR, "bed")
TABLE_DIR <- file.path(FIG2_DIR, "tables")
PLOT_DIR <- file.path(FIG2_DIR, "plots")
LOG_DIR <- file.path(FIG2_DIR, "logs")

COMPARTMENT_DIR <- file.path(FIG2_DIR, "compartments")
ENHANCER_DIR <- file.path(FIG2_DIR, "external_enhancers")
CIERNIA_DIR <- file.path(FIG2_DIR, "data/ciernia")

CIERNIA_HUB_DIR <- file.path(CIERNIA_DIR, "hub")
CIERNIA_BB_DIR <- file.path(CIERNIA_DIR, "bb")
CIERNIA_BED_DIR <- file.path(CIERNIA_DIR, "bed")

# Tool paths
BEDTOOLS <- Sys.which("bedtools")
BIGBEDTOBED <- file.path(UCSC_DIR, "bigBedToBed")

# Create output directories
for (d in c(
  FIG2_DIR,
  BED_DIR,
  TABLE_DIR,
  PLOT_DIR,
  LOG_DIR,
  COMPARTMENT_DIR,
  ENHANCER_DIR,
  CIERNIA_DIR,
  CIERNIA_HUB_DIR,
  CIERNIA_BB_DIR,
  CIERNIA_BED_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# INPUT CHECKS
# =============================================================================

check_file <- function(path, label = NULL) {
  if (is.null(label)) label <- basename(path)
  if (!file.exists(path)) {
    stop("Missing required file: ", label, "\nPath: ", path)
  }
  invisible(TRUE)
}

check_optional_file <- function(path, label = NULL) {
  if (is.null(label)) label <- basename(path)
  if (!file.exists(path)) {
    warning("Optional file not found: ", label, "\nPath: ", path)
    return(FALSE)
  }
  TRUE
}

check_executable <- function(path, label = NULL) {
  if (is.null(label)) label <- basename(path)
  if (is.na(path) || path == "" || !file.exists(path)) {
    stop("Missing executable: ", label, "\nPath: ", path)
  }
  if (file.access(path, mode = 1) != 0) {
    stop("File exists but is not executable: ", label, "\nPath: ", path)
  }
  invisible(TRUE)
}

check_inputs <- function(require_tools = FALSE) {
  check_file(OBJECT_RDS, "microglia_all_final.rds")
  check_file(DA_STAT1_CSV, "DA_global_STAT1_KO_vs_WT.csv")
  check_file(DA_IRF1_CSV, "DA_global_IRF1_KO_vs_WT.csv")
  check_optional_file(DA_WITHIN_CLUSTER_CSV, "DA_within_cluster_all.csv")

  if (require_tools) {
    if (BEDTOOLS == "") {
      stop("bedtools not found in PATH. Activate my_r_env first.")
    }
    check_executable(BEDTOOLS, "bedtools")
    check_executable(BIGBEDTOBED, "bigBedToBed")
  }

  invisible(TRUE)
}

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

GENOME <- "mm10"
PRIMARY_SEQLEVELS <- paste0("chr", c(1:19, "X", "Y"))
ASSAY_NAME <- "peakunion"

PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 2000

PROMOTER_SENS_UPSTREAM <- 1000
PROMOTER_SENS_DOWNSTREAM <- 1000

LOCUS_WINDOW <- 100000

FDR_CUTOFF <- 0.10
STRICT_FDR_CUTOFF <- 0.05

MIN_PEAKS_FOR_SUMMARY <- 1

# Main lipid-handling gene set
LIPID_GENES <- c(
  "Apoe", "Lpl", "Lipa", "Plin2", "Plin3", "Nceh1",
  "Sorl1", "Soat1", "Abca1", "Abcg1", "Trem2", "Tyrobp"
)

# Optional / secondary lipid genes
SECONDARY_LIPID_GENES <- c(
  "Acat1", "Acat2", "Cd36", "Fabp5", "Ctsb", "Ctsd", "Lamp1"
)

ALL_TARGET_GENES <- unique(c(LIPID_GENES, SECONDARY_LIPID_GENES))

# Controls: exclude these from candidate controls
EXCLUDE_CONTROL_GENES <- unique(c(ALL_TARGET_GENES))

# =============================================================================
# REPLICATE / INTERPRETATION NOTES
# =============================================================================

# Current scATAC mouse structure:
# TKK8 = STAT1_KO, 7516 cells
# TKK9 = IRF1_KO,  4818 cells
# TKL1 = IRF1_KO,  4344 cells
# TL41 = WT,       5811 cells
# TL53 = WT,        810 cells
#
# Interpretation:
# - STAT1_KO vs WT is directional/preliminary because STAT1_KO has n = 1 mouse.
# - IRF1_KO vs WT is more interpretable but still small, n = 2 vs n = 2.
# - STAT1_KO vs IRF1_KO should be treated as directional only.
#
# These scripts summarize peak-level DA into promoter/distal/locus accessibility
# architecture. They do not prove direct STAT1/IRF1 binding and do not assume
# same-cell multiome.

# =============================================================================
# HELPERS
# =============================================================================

load_microglia_object <- function(verbose = TRUE) {
  check_file(OBJECT_RDS, "microglia_all_final.rds")

  if (verbose) {
    message("Reading object: ", OBJECT_RDS)
  }

  obj <- readRDS(OBJECT_RDS)

  if (!ASSAY_NAME %in% Assays(obj)) {
    stop("Assay not found in object: ", ASSAY_NAME)
  }

  DefaultAssay(obj) <- ASSAY_NAME

  ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevels(ann) <- paste0("chr", seqlevels(ann))
  genome(ann) <- GENOME
  ann <- keepSeqlevels(ann, PRIMARY_SEQLEVELS, pruning.mode = "coarse")
  Annotation(obj) <- ann

  if (verbose) {
    message("Loaded object with dimensions:")
    print(dim(obj))
    message("Genotype table:")
    print(table(obj$genotype))
    message("Mouse table:")
    print(table(obj$mouse_id))
  }

  obj
}

make_gene_annotation <- function() {
  genes <- ensembldb::genes(
    EnsDb.Mmusculus.v79,
    columns = c("gene_id", "gene_name", "gene_biotype"),
    return.type = "GRanges"
  )

  if (!all(grepl("^chr", as.character(seqlevels(genes))))) {
    seqlevels(genes) <- paste0("chr", seqlevels(genes))
  }

  genome(genes) <- GENOME

  genes <- keepSeqlevels(
    genes,
    PRIMARY_SEQLEVELS,
    pruning.mode = "coarse"
  )

  genes <- genes[!is.na(genes$gene_name)]
  genes <- genes[genes$gene_name != ""]
  genes <- genes[!duplicated(genes$gene_name)]

  genes <- sortSeqlevels(genes)
  genes <- sort(genes)

  if (length(genes) == 0) {
    stop("make_gene_annotation() returned zero genes. Check EnsDb.Mmusculus.v79 / ensembldb.")
  }

  genes
}

clean_bed_gr <- function(gr) {
  gr <- keepSeqlevels(gr, PRIMARY_SEQLEVELS, pruning.mode = "coarse")
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  genome(gr) <- GENOME
  gr
}

peak_string_to_gr <- function(x) {
  Signac::StringToGRanges(x, sep = c("-", "-"))
}

gr_to_peak_string <- function(gr) {
  Signac::GRangesToString(gr, sep = c("-", "-"))
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

sig_lost <- function(logFC, FDR, fdr_cutoff = FDR_CUTOFF) {
  !is.na(FDR) & FDR < fdr_cutoff & !is.na(logFC) & logFC < 0
}

sig_gained <- function(logFC, FDR, fdr_cutoff = FDR_CUTOFF) {
  !is.na(FDR) & FDR < fdr_cutoff & !is.na(logFC) & logFC > 0
}

read_da_table <- function(path, comparison_name) {
  check_file(path, paste0("DA table: ", comparison_name))

  da <- read.csv(path, stringsAsFactors = FALSE)

  required_cols <- c("peak", "logFC", "FDR")
  missing_cols <- setdiff(required_cols, colnames(da))

  if (length(missing_cols) > 0) {
    stop(
      "DA table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nPath: ",
      path
    )
  }

  da$comparison <- comparison_name

  # Standardize optional columns if missing
  if (!"logCPM" %in% colnames(da)) da$logCPM <- NA_real_
  if (!"PValue" %in% colnames(da)) da$PValue <- NA_real_
  if (!"source" %in% colnames(da)) da$source <- "global"

  da
}

read_all_da_tables <- function() {
  stat1 <- read_da_table(
    DA_STAT1_CSV,
    comparison_name = "STAT1_KO_vs_WT"
  )

  irf1 <- read_da_table(
    DA_IRF1_CSV,
    comparison_name = "IRF1_KO_vs_WT"
  )

  bind_rows(stat1, irf1)
}

write_gr_bed <- function(gr, path) {
  gr <- clean_bed_gr(gr)
  rtracklayer::export(gr, path, format = "BED")
  message("Wrote BED: ", path, " [", length(gr), " ranges]")
  invisible(path)
}

read_bed_gr <- function(path) {
  check_file(path)
  gr <- rtracklayer::import(path, format = "BED")
  clean_bed_gr(gr)
}

run_cmd <- function(cmd) {
  message("Running command:\n", cmd)
  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ":\n", cmd)
  }
  invisible(TRUE)
}

save_csv <- function(x, path) {
  readr::write_csv(x, path)
  message("Wrote table: ", path, " [", nrow(x), " rows]")
  invisible(path)
}

# =============================================================================
# PRINT CONFIG SUMMARY WHEN SOURCED INTERACTIVELY
# =============================================================================

message("Loaded 00_config.R")
message("PROJECT: ", PROJECT)
message("SCATAC_DIR: ", SCATAC_DIR)
message("FIG2_DIR: ", FIG2_DIR)
message("OBJECT_RDS: ", OBJECT_RDS)
message("DA_STAT1_CSV: ", DA_STAT1_CSV)
message("DA_IRF1_CSV: ", DA_IRF1_CSV)
message("BEDTOOLS: ", ifelse(BEDTOOLS == "", "not found", BEDTOOLS))
message("BIGBEDTOBED: ", BIGBEDTOBED)