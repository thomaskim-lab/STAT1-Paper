############################################################
### FIGURE 6 (lipid arm) - Primary microglia STAT1 KD x Pioglitazone RNA-seq
### Plasmidsaurus order: YY2BYD
### Species: Mouse / Mus musculus
###
### Experimental design (2 x 2):
###   Ctrl - Pio   : PMG_Lipo_72h_R1/R2/R3
###   Ctrl + Pio   : PMG_Lipo_72h_Pioglitazone_R1/R2/R3
###   siSTAT1 - Pio: PMG_siSTAT1_3_72h_R1/R2/R3
###   siSTAT1 + Pio: PMG_siSTAT1_3_72h_Pioglitazone_R1/R2/R3
###
### This is the pioglitazone analogue of the BTKTYZ IFNg script, kept
### visually and structurally matched to Figure 4 / Figure 6.
###
### Biological questions (from the interpretation memo):
###   1) Does pioglitazone induce a CD36-dominant lipid-uptake / handling
###      program rather than classical de novo lipogenesis or adipogenesis?
###   2) What does siSTAT1 alone do to lipogenic vs cholesterol-handling arms?
###   3) Does STAT1 depletion modify the pioglitazone response
###      (STAT1 x Pioglitazone interaction)?
###
### CRITICAL CAVEAT BUILT INTO THE FIGURE:
###   The memo notes Stat1 mRNA is NOT reduced at 72 h. Panel 1 is a
###   knockdown-validation panel. If Stat1/ISGs are not suppressed, the
###   "siSTAT1" arm should be labelled "siSTAT1-transfected", not "STAT1 KD",
###   and STAT1-dependence claims must be withheld.
###
### Statistical design:
###   Model A (pairwise): ~ Replicate + Group   (falls back to ~ Group)
###   Model B (interaction): ~ Replicate + siRNA * Pio
###   Shrunken log2FC (ashr) for effect-size panels; raw Wald stat for GSEA.
############################################################

############################
### 0. Packages
############################

required_pkgs <- c(
  "DESeq2", "tidyverse", "pheatmap", "ggrepel", "openxlsx",
  "org.Mm.eg.db", "AnnotationDbi", "ashr", "fgsea", "limma",
  "patchwork", "svglite", "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall CRAN/Bioconductor packages before running this script."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(openxlsx)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ashr)
  library(fgsea)
  library(limma)
  library(patchwork)
  library(svglite)
})

set.seed(1234)

############################
### 1. Paths
############################

base_dir <- "C:/Users/au731006/Downloads/Plasmidsaurus"

matrix_file <- file.path(
  base_dir, "YY2BYD_results", "YY2BYD-expression-matrix.tsv"
)

sample_map_file <- file.path(
  base_dir, "YY2BYD_results", "YY2BYD.csv"
)

analysis_root <- file.path(
  base_dir, "YY2BYD_PrimaryMicroglia_STAT1_Pioglitazone_72h_analysis"
)

outdir   <- file.path(analysis_root, "Figure6_STAT1_KD_Pioglitazone_paper_ready")
figdir   <- file.path(outdir, "Figures")
tabledir <- file.path(outdir, "Tables")
objdir   <- file.path(outdir, "R_objects")

dir.create(figdir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tabledir, recursive = TRUE, showWarnings = FALSE)
dir.create(objdir,   recursive = TRUE, showWarnings = FALSE)

if (!file.exists(matrix_file)) stop("Expression matrix not found: ", matrix_file)

############################################################
### 2. Visual system - matched to Figure 4 / Figure 6
############################################################

COL_CONTROL  <- "#7F7F7F"   # Lipo only
COL_CTRL_PIO <- "#E08214"   # Lipo + pioglitazone
COL_KD       <- "#2166AC"   # siSTAT1
COL_KD_PIO   <- "#762A83"   # siSTAT1 + pioglitazone

COL_UP    <- "#B2182B"
COL_DOWN  <- "#2166AC"
COL_IFN   <- "#053061"      # STAT1/ISG accent
COL_LIPID <- "#E08214"
COL_REDOX <- "#762A83"
COL_BG    <- "grey82"

FDR_CUTOFF       <- 0.05
NOMINAL_P_CUTOFF <- 0.05
LFC_CUTOFF       <- 0.50   # volcano display threshold only
MIN_DIRECTIONAL_LFC <- 0.15 # descriptive rewiring classification only
COLOR_LIMIT         <- 1.50 # fixed effect-map colour range

GROUP_COLORS <- c(
  "Ctrl_NoPio" = COL_CONTROL,
  "Ctrl_Pio"   = COL_CTRL_PIO,
  "KD_NoPio"   = COL_KD,
  "KD_Pio"     = COL_KD_PIO
)

GROUP_LABELS <- c(
  "Ctrl_NoPio" = "Control",
  "Ctrl_Pio"   = "Control + Pio",
  "KD_NoPio"   = "siSTAT1",
  "KD_Pio"     = "siSTAT1 + Pio"
)

PIO_CONTRAST_COLORS <- c(
  "Pio effect in control" = COL_CTRL_PIO,
  "Pio effect in siSTAT1"  = COL_KD_PIO
)

HIGHLIGHT_COLORS <- c(
  "STAT1/ISG genes"      = COL_IFN,
  "Lipid-handling genes" = COL_LIPID,
  "Other genes"          = COL_BG
)

theme_fig <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(colour = "black"),
      axis.text  = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.spacing = grid::unit(2.5, "mm")
    )
}

save_plot <- function(p, stem, width, height) {
  svglite::svglite(
    filename = file.path(figdir, paste0(stem, ".svg")),
    width = width, height = height, bg = "transparent"
  )
  print(p); dev.off()
  
  pdf_file <- file.path(figdir, paste0(stem, ".pdf"))
  tryCatch(
    ggsave(pdf_file, p, width = width, height = height,
           device = grDevices::cairo_pdf, bg = "white"),
    error = function(e)
      ggsave(pdf_file, p, width = width, height = height,
             device = "pdf", bg = "white")
  )
  
  ggsave(file.path(figdir, paste0(stem, ".png")),
         p, width = width, height = height, dpi = 600, bg = "white")
}

make_safe_sheet_names <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  out <- character(length(x)); used <- character(0)
  for (i in seq_along(x)) {
    base <- substr(x[i], 1, 31); candidate <- base; counter <- 1
    while (tolower(candidate) %in% tolower(used)) {
      suffix <- paste0("_", counter)
      candidate <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
      counter <- counter + 1
    }
    out[i] <- candidate; used <- c(used, candidate)
  }
  out
}

write_xlsx_safe <- function(x, file) {
  names(x) <- make_safe_sheet_names(names(x))
  openxlsx::write.xlsx(x, file = file, overwrite = TRUE)
}

sig_mark <- function(padj, pvalue) {
  case_when(
    !is.na(padj)   & padj   < FDR_CUTOFF       ~ "*",
    !is.na(pvalue) & pvalue < NOMINAL_P_CUTOFF ~ "\u2020",
    TRUE ~ ""
  )
}

safe_neglog10 <- function(x, cap = 6) pmin(-log10(pmax(x, .Machine$double.xmin)), cap)

############################
### 3. Read expression matrix and sample map
############################

expr_df <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("gene_id", "gene_name", "gene_biotype") %in% colnames(expr_df)))

if (file.exists(sample_map_file)) {
  sample_map <- read.csv(sample_map_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  sample_map <- data.frame(
    run_id = paste0("YY2BYD_", 1:12),
    Sample_name = c(
      "PMG_Lipo_72h_R1", "PMG_Lipo_72h_R2", "PMG_Lipo_72h_R3",
      "PMG_Lipo_72h_Pioglitazone_R1", "PMG_Lipo_72h_Pioglitazone_R2",
      "PMG_Lipo_72h_Pioglitazone_R3",
      "PMG_siSTAT1_3_72h_R1", "PMG_siSTAT1_3_72h_R2", "PMG_siSTAT1_3_72h_R3",
      "PMG_siSTAT1_3_72h_Pioglitazone_R1", "PMG_siSTAT1_3_72h_Pioglitazone_R2",
      "PMG_siSTAT1_3_72h_Pioglitazone_R3"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(sample_map,
            file.path(tabledir, "Generated_YY2BYD_sample_map.csv"),
            row.names = FALSE)
}

if (!"run_id" %in% colnames(sample_map) && "Run_ID" %in% colnames(sample_map)) {
  sample_map <- dplyr::rename(sample_map, run_id = Run_ID)
}
stopifnot(all(c("run_id", "Sample_name") %in% colnames(sample_map)))

############################
### 4. Build count matrix
############################

count_cols <- grep("_count$", colnames(expr_df), value = TRUE)
if (length(count_cols) == 0) stop("No *_count columns found in expression matrix.")
strip_suffix <- function(x) sub("_(count|cpm)$", "", x)

count_mat <- as.matrix(expr_df[, count_cols, drop = FALSE])
rownames(count_mat) <- sub("\\..*$", "", expr_df$gene_id)
colnames(count_mat) <- strip_suffix(count_cols)
mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"
count_mat <- count_mat[!duplicated(rownames(count_mat)), , drop = FALSE]

############################
### 5. Match sequencing IDs to biological sample names
############################

matrix_ids <- colnames(count_mat)
if (all(matrix_ids %in% sample_map$run_id)) {
  sample_map <- sample_map[match(matrix_ids, sample_map$run_id), , drop = FALSE]
  colnames(count_mat) <- sample_map$Sample_name
} else if (all(matrix_ids %in% sample_map$Sample_name)) {
  sample_map <- sample_map[match(matrix_ids, sample_map$Sample_name), , drop = FALSE]
  colnames(count_mat) <- sample_map$Sample_name
} else {
  stop("Expression-matrix sample IDs do not match YY2BYD.csv.\nMatrix IDs: ",
       paste(matrix_ids, collapse = ", "))
}

############################
### 6. Metadata: explicit 2 x 2 design
############################

parse_metadata <- function(sample_names) {
  tibble(Sample_name = sample_names) %>%
    mutate(
      Group = case_when(
        grepl("^PMG_Lipo_72h_R[0-9]+$", Sample_name)                    ~ "Ctrl_NoPio",
        grepl("^PMG_Lipo_72h_Pioglitazone_R[0-9]+$", Sample_name)       ~ "Ctrl_Pio",
        grepl("^PMG_siSTAT1_3_72h_R[0-9]+$", Sample_name)               ~ "KD_NoPio",
        grepl("^PMG_siSTAT1_3_72h_Pioglitazone_R[0-9]+$", Sample_name)  ~ "KD_Pio",
        TRUE ~ NA_character_
      ),
      siRNA = case_when(
        Group %in% c("Ctrl_NoPio", "Ctrl_Pio") ~ "Control",
        Group %in% c("KD_NoPio", "KD_Pio")     ~ "STAT1_KD",
        TRUE ~ NA_character_
      ),
      Pio = case_when(
        Group %in% c("Ctrl_Pio", "KD_Pio")     ~ "Pio",
        Group %in% c("Ctrl_NoPio", "KD_NoPio") ~ "No_Pio",
        TRUE ~ NA_character_
      ),
      Replicate = sub(".*_R([0-9]+)$", "R\\1", Sample_name)
    ) %>%
    mutate(
      Group     = factor(Group, levels = c("Ctrl_NoPio", "Ctrl_Pio", "KD_NoPio", "KD_Pio")),
      siRNA     = factor(siRNA, levels = c("Control", "STAT1_KD")),
      Pio       = factor(Pio,   levels = c("No_Pio", "Pio")),
      Replicate = factor(Replicate, levels = sort(unique(Replicate)))
    )
}

meta <- parse_metadata(colnames(count_mat))
rownames(meta) <- meta$Sample_name

if (any(is.na(meta$Group)))
  stop("Unrecognized sample names:\n",
       paste(meta$Sample_name[is.na(meta$Group)], collapse = "\n"))

stopifnot(all(colnames(count_mat) == rownames(meta)))

# Is this a complete matched 2x2 (one sample per group per replicate block)?
block_table <- table(meta$Replicate, meta$Group)
print(block_table)
matched_blocks <- all(block_table == 1) && nrow(block_table) >= 2
if (matched_blocks) {
  message("Matched 2x2 block design detected -> using ~ Replicate + ... .")
  message("NOTE: this assumes R1/R2/R3 are PAIRED biological blocks. ",
          "If R1-R3 are independent replicates, set FORCE_NO_BLOCK <- TRUE below.")
} else {
  message("Design is NOT a complete matched 2x2 -> dropping Replicate term.")
}

FORCE_NO_BLOCK <- FALSE   # set TRUE if replicates are not paired blocks
use_block <- matched_blocks && !FORCE_NO_BLOCK

design_group       <- if (use_block) ~ Replicate + Group      else ~ Group
design_factorial   <- if (use_block) ~ Replicate + siRNA * Pio else ~ siRNA * Pio

write.xlsx(as.data.frame(meta),
           file.path(tabledir, "Figure6_Pio_sample_metadata.xlsx"),
           overwrite = TRUE)

############################
### 7. Gene annotation
############################

source_annot <- expr_df %>%
  transmute(ENSEMBL = sub("\\..*$", "", gene_id),
            source_gene_name = gene_name,
            gene_biotype = gene_biotype) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

orgdb_annot <- AnnotationDbi::select(
  org.Mm.eg.db, keys = unique(rownames(count_mat)),
  keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID", "GENENAME")
) %>% distinct(ENSEMBL, .keep_all = TRUE)

gene_annot <- tibble(ENSEMBL = rownames(count_mat)) %>%
  left_join(orgdb_annot, by = "ENSEMBL") %>%
  left_join(source_annot, by = "ENSEMBL") %>%
  mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "", source_gene_name, SYMBOL))

############################
### 8. Filtering
############################

keep <- rowSums(count_mat >= 10) >= 3
count_filt <- count_mat[keep, , drop = FALSE]
message("Genes before filtering: ", nrow(count_mat),
        " | retained: ", nrow(count_filt))

############################
### 9. Model A: pairwise contrasts (with ashr shrinkage)
############################

dds_group <- DESeqDataSetFromMatrix(count_filt, as.data.frame(meta), design_group)
dds_group <- DESeq(dds_group)

# Build a results table carrying raw Wald stat + shrunken log2FC.
make_res_table <- function(dds_obj, num, den, contrast_name) {
  raw <- results(dds_obj, contrast = c("Group", num, den), alpha = 0.05)
  shr <- lfcShrink(dds_obj, contrast = c("Group", num, den), res = raw, type = "ashr")
  
  raw_df <- as.data.frame(raw) %>%
    rownames_to_column("ENSEMBL") %>%
    transmute(ENSEMBL, baseMean,
              log2FC_raw = log2FoldChange, lfcSE_raw = lfcSE,
              stat, pvalue, padj)
  
  shr_df <- as.data.frame(shr) %>%
    rownames_to_column("ENSEMBL") %>%
    transmute(ENSEMBL, log2FoldChange, lfcSE)
  
  raw_df %>%
    left_join(shr_df, by = "ENSEMBL") %>%
    left_join(gene_annot, by = "ENSEMBL") %>%
    mutate(Contrast = contrast_name) %>%
    arrange(padj, pvalue)
}

contrast_tables <- list(
  Pio_in_Control     = make_res_table(dds_group, "Ctrl_Pio", "Ctrl_NoPio", "Pio_in_Control"),
  Pio_in_siSTAT1     = make_res_table(dds_group, "KD_Pio",   "KD_NoPio",   "Pio_in_siSTAT1"),
  siSTAT1_no_Pio     = make_res_table(dds_group, "KD_NoPio", "Ctrl_NoPio", "siSTAT1_no_Pio"),
  siSTAT1_with_Pio   = make_res_table(dds_group, "KD_Pio",   "Ctrl_Pio",   "siSTAT1_with_Pio"),
  combined_vs_ctrl   = make_res_table(dds_group, "KD_Pio",   "Ctrl_NoPio", "combined_vs_ctrl")
)

############################
### 10. Model B: STAT1 x Pioglitazone interaction
############################

dds_fac <- DESeqDataSetFromMatrix(count_filt, as.data.frame(meta), design_factorial)
dds_fac <- DESeq(dds_fac)
rn <- resultsNames(dds_fac); print(rn)

interaction_name <- rn[
  grepl("siRNA", rn, ignore.case = TRUE) &
    grepl("Pio",   rn, ignore.case = TRUE) &
    !grepl("_vs_", rn)
]
if (length(interaction_name) != 1)
  stop("Could not uniquely identify STAT1 x Pio interaction coefficient.\n",
       "resultsNames: ", paste(rn, collapse = ", "))
message("Interaction coefficient: ", interaction_name)

res_int_raw <- results(dds_fac, name = interaction_name, alpha = 0.05)
res_int_shr <- lfcShrink(dds_fac, coef = interaction_name, res = res_int_raw, type = "ashr")

interaction_table <- as.data.frame(res_int_raw) %>%
  rownames_to_column("ENSEMBL") %>%
  transmute(ENSEMBL, baseMean,
            log2FC_raw = log2FoldChange, lfcSE_raw = lfcSE, stat, pvalue, padj) %>%
  left_join(as.data.frame(res_int_shr) %>%
              rownames_to_column("ENSEMBL") %>%
              transmute(ENSEMBL, log2FoldChange, lfcSE), by = "ENSEMBL") %>%
  left_join(gene_annot, by = "ENSEMBL") %>%
  mutate(Contrast = "STAT1_x_Pio_interaction") %>%
  arrange(padj, pvalue)

contrast_tables$STAT1_x_Pio_interaction <- interaction_table

write_xlsx_safe(contrast_tables,
                file.path(tabledir, "Figure6_Pio_DESeq2_all_contrasts.xlsx"))

############################
### 11. VST + QC
############################

vsd <- vst(dds_group, blind = FALSE)
vst_mat <- assay(vsd)
saveRDS(dds_group, file.path(objdir, "dds_group.rds"))
saveRDS(dds_fac,   file.path(objdir, "dds_factorial.rds"))
saveRDS(vsd,       file.path(objdir, "vsd.rds"))

write.xlsx(as.data.frame(counts(dds_group, normalized = TRUE)) %>%
             rownames_to_column("ENSEMBL"),
           file.path(tabledir, "Figure6_Pio_normalized_counts.xlsx"), overwrite = TRUE)

pca_df <- plotPCA(vsd, intgroup = c("Group", "Replicate"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group, shape = Replicate)) +
  geom_point(size = 3.3, stroke = 0.8) +
  scale_colour_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = NULL) +
  xlab(paste0("PC1 (", percent_var[1], "%)")) +
  ylab(paste0("PC2 (", percent_var[2], "%)")) +
  labs(title = "Primary microglia: siSTAT1 x Pioglitazone (72 h)") +
  theme_fig(8) + theme(legend.position = "top")
save_plot(p_pca, "Fig6_Pio_QC_PCA", 5.0, 4.0)

sample_dist <- as.matrix(dist(t(vst_mat)))
anno <- meta %>% dplyr::select(Group, Replicate) %>% as.data.frame()
rownames(anno) <- rownames(meta)
distance_cols <- colorRampPalette(
  c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100)
annotation_colors_qc <- list(
  Group = GROUP_COLORS,
  Replicate = c("R1" = "#1B9E77", "R2" = "#D95F02", "R3" = "#7570B3"))

pdf(file.path(figdir, "Fig6_Pio_QC_sample_distance.pdf"), width = 6.2, height = 5.6)
pheatmap(sample_dist, annotation_col = anno, annotation_row = anno,
         annotation_colors = annotation_colors_qc, color = distance_cols,
         border_color = "white", fontsize = 7, main = "Sample-to-sample distance")
dev.off()

############################
### 12. Curated modules (ordered by the memo's pathway ranking)
############################

pio_modules <- list(
  Canonical_STAT1_ISG = c(
    "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15",
    "Ifit1", "Ifit2", "Ifit3", "Ifitm1", "Ifitm2", "Ifitm3",
    "Oas1a", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2",
    "Mx1", "Mx2", "Usp18", "Rsad2", "Zbp1", "Ddx58", "Ifih1",
    "Cxcl10", "Gbp2", "Socs1"
  ),
  Inflammatory_cytokine_chemokine = c(
    "Tnf", "Il1b", "Il6", "Cxcl10", "Cxcl9", "Ccl2", "Ccl3", "Ccl4",
    "Nfkbia", "Nfkbiz", "Tnfaip3", "Ptgs2", "Nos2"
  ),
  Antigen_presentation_MHC = c(
    "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "Cd74", "Ciita",
    "B2m", "Tap1", "Tap2", "Psmb8", "Psmb9"
  ),
  Homeostatic_microglia = c(
    "P2ry12", "Tmem119", "Sall1", "Siglech", "Fcrls", "Cx3cr1", "Olfml3",
    "Gpr34", "P2ry13", "Mef2a", "Mef2c"
  ),
  DAM_lipid_associated = c(
    "Apoe", "Lpl", "Trem2", "Tyrobp", "Spp1", "Lgals3", "Cst7", "Ctsb",
    "Ctsd", "Ctsl", "Lyz2", "Itgax", "Clec7a", "Axl", "Cd68"
  ),
  FA_uptake_transport = c(
    "Cd36", "Lpl", "Lipa", "Slc27a1", "Slc27a4", "Fabp3", "Fabp5",
    "Msr1", "Lrp1", "Ldlr", "Vldlr", "Apoe", "Nceh1"
  ),
  Lysosome_lipophagy = c(
    "Lipa", "Lamp1", "Lamp2", "Ctsb", "Ctsd", "Ctsl", "Rab7a", "Atg5",
    "Atg7", "Becn1", "Sqstm1", "Map1lc3b", "Npc1", "Npc2"
  ),
  Lipid_droplet_TAG = c(
    "Plin1", "Plin2", "Plin3", "Plin4", "Plin5", "Dgat1", "Dgat2",
    "Lpin1", "Lpin2", "Soat1", "Soat2", "Mogat1", "Fitm2", "Cidec"
  ),
  Cholesterol_transport_efflux = c(
    "Abca1", "Abcg1", "Scarb1", "Soat1", "Npc1", "Npc2", "Nceh1",
    "Ch25h", "Apoe", "Lipa"
  ),
  FA_synthesis_desaturation = c(
    "Fasn", "Acaca", "Acacb", "Scd1", "Scd2", "Elovl6", "Me1", "Srebf1",
    "Acly", "Lpcat3"
  ),
  FA_oxidation = c(
    "Pnpla2", "Abhd5", "Lipe", "Mgll", "Cpt1a", "Cpt2", "Acadm", "Acadvl",
    "Hadha", "Hadhb", "Acox1", "Ppara", "Ppargc1a"
  ),
  Cholesterol_biosynthesis = c(
    "Hmgcr", "Hmgcs1", "Idi1", "Fdps", "Fdft1", "Sqle", "Dhcr7", "Dhcr24",
    "Msmo1", "Srebf2"
  ),
  Sphingolipid_ceramide = c(
    "Sptlc1", "Sptlc2", "Cers5", "Cers6", "Ugcg", "Asah2", "Sgms1", "Smpd1"
  ),
  Eicosanoid = c(
    "Ptgs1", "Ptgs2", "Ptges", "Alox5", "Alox15", "Pla2g4a"
  ),
  PPAR_adipogenic = c(
    "Pparg", "Cd36", "Lpl", "Plin2", "Scd1", "Dgat1", "Acaca", "Fabp4",
    "Cpt1a", "Angptl4", "Pdk4", "Ucp2"
  ),
  OXPHOS_ETC = c(
    "Ndufa1", "Ndufa2", "Ndufa9", "Ndufb8", "Ndufs1", "Ndufs2", "Sdha",
    "Sdhb", "Uqcrc1", "Uqcrc2", "Cox4i1", "Cox5a", "Cox6a1", "Atp5f1a",
    "Atp5f1b", "Atp5f1c", "Atp5mc1"
  ),
  Glycolysis = c(
    "Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkl", "Pfkp", "Aldoa", "Gapdh",
    "Pgk1", "Pgam1", "Eno1", "Pkm", "Ldha"
  ),
  Iron_redox_glutathione = c(
    "Fth1", "Ftl1", "Hmox1", "Slc40a1", "Tfrc", "Slc11a2", "Gpx4", "Gpx1",
    "Gclc", "Gclm", "Gsr", "Gss", "Txn1", "Txnrd1", "Prdx1", "Prdx2",
    "Nqo1", "Sod1", "Sod2"
  )
)

module_labels <- c(
  Canonical_STAT1_ISG = "STAT1 / ISG (KD control)",
  Inflammatory_cytokine_chemokine = "Inflammatory cytokines",
  Antigen_presentation_MHC = "Antigen presentation / MHC",
  Homeostatic_microglia = "Homeostatic microglia",
  DAM_lipid_associated = "DAM / lipid-associated",
  FA_uptake_transport = "FA uptake / transport",
  Lysosome_lipophagy = "Lysosome / lipophagy",
  Lipid_droplet_TAG = "Lipid droplet / TAG",
  Cholesterol_transport_efflux = "Cholesterol transport / efflux",
  FA_synthesis_desaturation = "FA synthesis / desaturation",
  FA_oxidation = "Lipolysis / FA oxidation",
  Cholesterol_biosynthesis = "Cholesterol biosynthesis",
  Sphingolipid_ceramide = "Sphingolipid / ceramide",
  Eicosanoid = "Eicosanoid",
  PPAR_adipogenic = "PPAR / adipogenic lipid",
  OXPHOS_ETC = "OXPHOS / ETC",
  Glycolysis = "Glycolysis",
  Iron_redox_glutathione = "Iron / redox / glutathione"
)

# Priority genes used in effect maps, volcanoes, labels, and heatmaps.
isg_priority <- c(
  "Stat1", "Irf1", "Irf7", "Isg15", "Ifit1", "Ifit2", "Ifit3",
  "Oasl2", "Cxcl10", "Gbp2", "Usp18"
)

lipid_gene_groups <- list(
  "Disease/lipid-associated axis" = c(
    "Apoe", "Lpl", "Trem2", "Tyrobp", "Nceh1", "Sorl1"
  ),
  "Lipid-state regulation" = c(
    "Pparg", "Abhd5", "Abcg1", "Ch25h", "Soat1"
  ),
  "Lysosomal processing / efflux" = c(
    "Lipa", "Abca1", "Npc1", "Npc2", "Scarb1"
  ),
  "Uptake / droplet organisation" = c(
    "Cd36", "Msr1", "Lrp1", "Plin2", "Plin3", "Dgat1", "Lpin1"
  ),
  "Fatty-acid use / lipolysis" = c(
    "Cpt1a", "Pnpla2", "Lipe"
  ),
  "Synthesis / desaturation" = c(
    "Fasn", "Acaca", "Scd1", "Elovl6", "Lpcat3", "Srebf1"
  ),
  "Cholesterol biosynthesis" = c(
    "Hmgcr", "Hmgcs1", "Idi1", "Fdps", "Sqle", "Dhcr7"
  ),
  "Sphingolipid / eicosanoid" = c(
    "Sptlc2", "Ugcg", "Asah2", "Ptgs2"
  )
)

lipid_priority <- unique(unlist(lipid_gene_groups, use.names = FALSE))
lipid_group_map <- enframe(
  lipid_gene_groups,
  name = "GeneGroup",
  value = "SYMBOL"
) %>%
  unnest(SYMBOL) %>%
  distinct(SYMBOL, .keep_all = TRUE)

redox_priority <- c(
  "Fth1", "Ftl1", "Gpx4", "Gclc", "Gclm", "Gsr", "Txnrd1",
  "Slc40a1", "Hmox1", "Nqo1"
)

############################
### 13. Panel 1: KNOCKDOWN VALIDATION (must read first)
############################

extract_effects <- function(genes, contrasts) {
  bind_rows(lapply(contrasts, function(ct) {
    contrast_tables[[ct]] %>%
      filter(SYMBOL %in% genes, !is.na(SYMBOL)) %>%
      group_by(SYMBOL) %>%
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(SYMBOL, log2FoldChange, lfcSE, pvalue, padj, Contrast = ct)
  }))
}

kd_df <- extract_effects(isg_priority, c("siSTAT1_no_Pio", "siSTAT1_with_Pio")) %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = isg_priority),
    Contrast = recode(Contrast,
                      siSTAT1_no_Pio   = "siSTAT1 vs Ctrl",
                      siSTAT1_with_Pio = "siSTAT1+Pio vs Ctrl+Pio"),
    Contrast = factor(Contrast, levels = c("siSTAT1 vs Ctrl", "siSTAT1+Pio vs Ctrl+Pio")),
    mark = sig_mark(padj, pvalue)
  )

p_kd <- ggplot(kd_df, aes(SYMBOL, log2FoldChange, fill = Contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.76), width = 0.66) +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
                position = position_dodge(width = 0.76), width = 0.18,
                linewidth = 0.28, na.rm = TRUE) +
  geom_text(aes(label = mark,
                y = log2FoldChange +
                  if_else(log2FoldChange >= 0, 1, -1) * (replace_na(lfcSE, 0) + 0.10)),
            position = position_dodge(width = 0.76), size = 2.6, show.legend = FALSE) +
  scale_fill_manual(values = c("siSTAT1 vs Ctrl" = COL_KD,
                               "siSTAT1+Pio vs Ctrl+Pio" = COL_KD_PIO), name = NULL) +
  labs(title = "Knockdown validation: Stat1 + canonical ISGs",
       x = NULL, y = expression(log[2]~"fold change (shrunk)")) +
  theme_fig(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                       legend.position = "top")
save_plot(p_kd, "Fig6_Pio_KD_validation_STAT1_ISG", 5.0, 3.2)

############################
### 14. Formal module statistics: fgsea + cameraPR
############################

make_rank_vector <- function(res_df) {
  ranked <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(stat)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>% arrange(desc(stat))
  stats <- ranked$stat; names(stats) <- ranked$SYMBOL; stats
}

score_modules <- function(res_df, module_list, contrast_label) {
  rank_vec <- make_rank_vector(res_df)
  fg <- fgsea::fgseaMultilevel(module_list, rank_vec,
                               minSize = 5, maxSize = 500, eps = 0) %>%
    as_tibble() %>%
    transmute(Module = pathway, fgsea_NES = NES,
              fgsea_pvalue = pval, fgsea_FDR = padj, n_fgsea = size)
  
  camera_tbl <- tryCatch({
    idx <- limma::ids2indices(module_list, names(rank_vec))
    limma::cameraPR(rank_vec, idx, use.ranks = TRUE, sort = FALSE) %>%
      as.data.frame() %>% rownames_to_column("Module") %>% as_tibble() %>%
      transmute(Module, camera_Direction = Direction,
                camera_PValue = PValue, camera_FDR = FDR)
  }, error = function(e) tibble(Module = names(module_list),
                                camera_Direction = NA_character_,
                                camera_PValue = NA_real_, camera_FDR = NA_real_))
  
  full_join(fg, camera_tbl, by = "Module") %>%
    mutate(Contrast = contrast_label, ModuleLabel = unname(module_labels[Module]))
}

module_scores <- bind_rows(
  score_modules(contrast_tables$Pio_in_Control, pio_modules, "Pio in control"),
  score_modules(contrast_tables$Pio_in_siSTAT1, pio_modules, "Pio in siSTAT1"),
  score_modules(contrast_tables$siSTAT1_no_Pio, pio_modules, "siSTAT1 alone")
)

write_xlsx_safe(list(module_scores = module_scores),
                file.path(tabledir, "Figure6_Pio_module_scores.xlsx"))

module_order <- names(module_labels)
main_module_order <- c(
  "Canonical_STAT1_ISG",
  "DAM_lipid_associated",
  "FA_uptake_transport",
  "Lipid_droplet_TAG",
  "Lysosome_lipophagy",
  "FA_oxidation",
  "OXPHOS_ETC",
  "Iron_redox_glutathione"
)

make_module_dotplot <- function(df, selected_modules) {
  plot_df <- df %>%
    filter(Module %in% selected_modules) %>%
    mutate(
      ModuleLabel = factor(
        ModuleLabel,
        levels = rev(unname(module_labels[selected_modules]))
      ),
      Contrast = factor(
        Contrast,
        levels = c("Pio in control", "Pio in siSTAT1", "siSTAT1 alone")
      ),
      plot_FDR = pmax(replace_na(fgsea_FDR, 1), 1e-4),
      camera_outline = case_when(
        !is.na(camera_FDR) & camera_FDR < 0.05 ~ "camera FDR < 0.05",
        !is.na(camera_PValue) & camera_PValue < 0.05 ~ "camera P < 0.05",
        TRUE ~ "not camera-significant"
      )
    )

  ggplot(
    plot_df,
    aes(
      Contrast, ModuleLabel,
      fill = fgsea_NES,
      size = safe_neglog10(plot_FDR, 6),
      colour = camera_outline
    )
  ) +
    geom_point(shape = 21, stroke = 0.6) +
    scale_fill_gradient2(
      low = COL_DOWN, mid = "white", high = COL_UP,
      midpoint = 0, name = "fgsea NES"
    ) +
    scale_size_continuous(
      range = c(1.8, 6.5),
      name = expression(-log[10]~"FDR")
    ) +
    scale_colour_manual(
      values = c(
        "camera FDR < 0.05" = "black",
        "camera P < 0.05" = "grey25",
        "not camera-significant" = "grey70"
      ),
      name = "cameraPR"
    ) +
    labs(x = NULL, y = NULL) +
    theme_fig(8) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right"
    )
}

p_modules <- make_module_dotplot(module_scores, main_module_order)
save_plot(p_modules, "Fig6_Pio_module_dotplot", 6.2, 5.4)

p_modules_full <- make_module_dotplot(module_scores, module_order)
save_plot(p_modules_full, "Fig6_Pio_full_module_dotplot", 6.4, 8.5)

############################
### 15. Lipid-gene effect map across pairwise contrasts + interaction
############################

map_contrasts <- c(
  "Pio_in_Control", "Pio_in_siSTAT1",
  "siSTAT1_no_Pio", "STAT1_x_Pio_interaction"
)
map_labels <- c(
  Pio_in_Control = "Pio in ctrl",
  Pio_in_siSTAT1 = "Pio in siSTAT1",
  siSTAT1_no_Pio = "siSTAT1 alone",
  STAT1_x_Pio_interaction = "STAT1 x Pio"
)

lipid_map_df <- bind_rows(lapply(map_contrasts, function(ct) {
  contrast_tables[[ct]] %>%
    filter(SYMBOL %in% lipid_priority, !is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(SYMBOL, log2FoldChange, lfcSE, pvalue, padj, Contrast = ct)
})) %>%
  complete(SYMBOL = lipid_priority, Contrast = map_contrasts) %>%
  left_join(lipid_group_map, by = "SYMBOL") %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = rev(lipid_priority)),
    GeneGroup = factor(GeneGroup, levels = names(lipid_gene_groups)),
    Contrast = factor(map_labels[Contrast], levels = unname(map_labels[map_contrasts])),
    mark = sig_mark(padj, pvalue),
    cell_label = if_else(
      is.na(log2FoldChange), "\u00D7",
      paste0(sprintf("%.2f", log2FoldChange), mark)
    ),
    lower95 = log2FoldChange - 1.96 * lfcSE,
    upper95 = log2FoldChange + 1.96 * lfcSE
  )

p_lipid_map <- ggplot(
  lipid_map_df,
  aes(Contrast, SYMBOL, fill = log2FoldChange)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 2.0) +
  facet_grid(
    rows = vars(GeneGroup),
    scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_fill_gradient2(
    low = COL_DOWN, mid = "white", high = COL_UP,
    midpoint = 0, limits = c(-COLOR_LIMIT, COLOR_LIMIT),
    oob = scales::squish, na.value = "grey88",
    name = "Shrunken\nlog2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_fig(7.5) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing.y = grid::unit(1.2, "mm"),
    legend.position = "right"
  )

save_plot(p_lipid_map, "Fig6_Pio_lipid_gene_effect_map", 6.8, 10.0)

write.xlsx(
  lipid_map_df %>% mutate(
    SYMBOL = as.character(SYMBOL),
    GeneGroup = as.character(GeneGroup),
    Contrast = as.character(Contrast)
  ),
  file.path(tabledir, "Figure6_Pio_lipid_gene_effect_map.xlsx"),
  overwrite = TRUE
)

############################
### 16. Concordance: pioglitazone effect in control vs in siSTAT1
############################

isg_set   <- pio_modules$Canonical_STAT1_ISG
lipid_set <- unique(unlist(pio_modules[c(
  "DAM_lipid_associated", "FA_uptake_transport", "Lysosome_lipophagy",
  "Lipid_droplet_TAG", "Cholesterol_transport_efflux",
  "FA_synthesis_desaturation", "FA_oxidation",
  "Cholesterol_biosynthesis", "Sphingolipid_ceramide",
  "Eicosanoid", "PPAR_adipogenic"
)]))

concordance_df <- contrast_tables$Pio_in_Control %>%
  transmute(ENSEMBL, SYMBOL, ctrl_log2FC = log2FoldChange, ctrl_padj = padj) %>%
  inner_join(contrast_tables$Pio_in_siSTAT1 %>%
               transmute(ENSEMBL, kd_log2FC = log2FoldChange, kd_padj = padj),
             by = "ENSEMBL") %>%
  filter(!is.na(ctrl_log2FC), !is.na(kd_log2FC)) %>%
  mutate(
    Highlight = case_when(
      SYMBOL %in% isg_set   ~ "STAT1/ISG genes",
      SYMBOL %in% lipid_set ~ "Lipid-handling genes",
      TRUE ~ "Other genes"),
    label_gene = if_else(
      SYMBOL %in% c("Cd36", "Lpl", "Plin2", "Scd1", "Dgat1", "Cpt1a",
                    "Acaca", "Abca1", "Abcg1", "Pparg", "Stat1"),
      SYMBOL, NA_character_)
  )

rho <- suppressWarnings(cor(concordance_df$ctrl_log2FC, concordance_df$kd_log2FC,
                            method = "spearman", use = "complete.obs"))

p_conc <- ggplot(concordance_df, aes(ctrl_log2FC, kd_log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = "dotted") +
  geom_point(data = filter(concordance_df, Highlight == "Other genes"),
             colour = COL_BG, alpha = 0.45, size = 0.8) +
  geom_point(data = filter(concordance_df, Highlight != "Other genes"),
             aes(colour = Highlight), size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = HIGHLIGHT_COLORS,
                      breaks = c("STAT1/ISG genes", "Lipid-handling genes"), name = NULL) +
  ggrepel::geom_text_repel(aes(label = label_gene), na.rm = TRUE, size = 2.2,
                           max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2) +
  annotate("text", x = Inf, y = -Inf, label = paste0("Spearman rho = ", round(rho, 2)),
           hjust = 1.05, vjust = -0.7, size = 2.7) +
  labs(x = expression("Pio effect in control ("*log[2]*"FC)"),
       y = expression("Pio effect in siSTAT1 ("*log[2]*"FC)")) +
  theme_fig(8) + theme(legend.position = "top")
save_plot(p_conc, "Fig6_Pio_response_control_vs_siSTAT1", 4.8, 4.4)

############################
### 17. Descriptive VST module scores across the four groups
############################

vst_symbol <- as.data.frame(vst_mat) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(gene_annot %>% dplyr::select(ENSEMBL, SYMBOL), by = "ENSEMBL") %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(SYMBOL, .keep_all = TRUE)

module_score_list <- lapply(names(pio_modules), function(mod) {
  genes <- intersect(pio_modules[[mod]], vst_symbol$SYMBOL)
  if (length(genes) < 3) return(NULL)
  sub <- vst_symbol %>% filter(SYMBOL %in% genes) %>%
    dplyr::select(-ENSEMBL, -SYMBOL) %>% as.matrix()
  z <- t(scale(t(sub))); z[is.na(z)] <- 0
  tibble(Sample_name = colnames(z), Score = colMeans(z, na.rm = TRUE),
         Module = mod, n_genes = nrow(z))
})

module_vst_scores <- bind_rows(module_score_list) %>%
  left_join(meta %>% as.data.frame(), by = "Sample_name")

write.xlsx(module_vst_scores,
           file.path(tabledir, "Figure6_Pio_descriptive_VST_module_scores.xlsx"),
           overwrite = TRUE)

focus_modules <- c("Canonical_STAT1_ISG", "FA_uptake_transport",
                   "Lipid_droplet_TAG", "PPAR_adipogenic")

p_scores <- module_vst_scores %>%
  filter(Module %in% focus_modules) %>%
  mutate(Group = factor(Group, levels = c("Ctrl_NoPio", "Ctrl_Pio", "KD_NoPio", "KD_Pio")),
         Module = factor(Module, levels = focus_modules,
                         labels = module_labels[focus_modules])) %>%
  ggplot(aes(Group, Score, group = Replicate)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_line(colour = "grey72", linewidth = 0.45) +
  geom_point(aes(colour = Group), size = 2.4) +
  stat_summary(aes(group = Group, colour = Group), fun = mean,
               geom = "crossbar", width = 0.48, linewidth = 0.55) +
  scale_colour_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = NULL) +
  scale_x_discrete(labels = GROUP_LABELS) +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Relative module score") +
  theme_fig(8) + theme(axis.text.x = element_text(angle = 28, hjust = 1),
                       legend.position = "none")
save_plot(p_scores, "Fig6_Pio_descriptive_module_scores", 6.8, 5.3)

############################################################
### 18. Fig4-style descriptive lipid-response rewiring map
############################################################

# This classification is descriptive. It compares the Pio response in control
# with the Pio response in siSTAT1 and uses the interaction coefficient as the
# formal difference-of-differences estimate. It does not replace interaction P/FDR.
direction_code <- function(x, threshold = MIN_DIRECTIONAL_LFC) {
  case_when(
    is.na(x) ~ NA_integer_,
    x >= threshold ~ 1L,
    x <= -threshold ~ -1L,
    TRUE ~ 0L
  )
}

lipid_rewiring_wide <- lipid_map_df %>%
  mutate(
    SYMBOL = as.character(SYMBOL),
    GeneGroup = as.character(GeneGroup),
    Contrast_key = recode(
      as.character(Contrast),
      "Pio in ctrl" = "ctrl",
      "Pio in siSTAT1" = "kd",
      "siSTAT1 alone" = "kd_alone",
      "STAT1 x Pio" = "interaction"
    )
  ) %>%
  dplyr::select(
    SYMBOL, GeneGroup, Contrast_key,
    log2FoldChange, lfcSE, pvalue, padj
  ) %>%
  pivot_wider(
    names_from = Contrast_key,
    values_from = c(log2FoldChange, lfcSE, pvalue, padj)
  ) %>%
  mutate(
    dir_ctrl = direction_code(log2FoldChange_ctrl),
    dir_kd = direction_code(log2FoldChange_kd),
    dir_interaction = direction_code(log2FoldChange_interaction),
    delta_abs = abs(log2FoldChange_kd) - abs(log2FoldChange_ctrl),
    RewiringPattern = case_when(
      is.na(dir_ctrl) | is.na(dir_kd) | is.na(dir_interaction) ~
        "Not detected in all contrasts",
      dir_ctrl != 0L & dir_kd != 0L & dir_ctrl != dir_kd ~
        "Direction reversed in siSTAT1",
      dir_ctrl != 0L & dir_kd == 0L ~
        "Pio response lost/attenuated",
      dir_ctrl == 0L & dir_kd != 0L ~
        "Pio response gained in siSTAT1",
      dir_ctrl != 0L & dir_kd == dir_ctrl &
        delta_abs >= MIN_DIRECTIONAL_LFC & dir_interaction == dir_ctrl ~
        "Pio response potentiated",
      dir_ctrl != 0L & dir_kd == dir_ctrl &
        delta_abs <= -MIN_DIRECTIONAL_LFC & dir_interaction == -dir_ctrl ~
        "Pio response attenuated",
      dir_ctrl != 0L & dir_kd == dir_ctrl ~
        "Pio response preserved",
      TRUE ~ "Weak/mixed"
    ),
    Interaction_nominal =
      !is.na(pvalue_interaction) & pvalue_interaction < NOMINAL_P_CUTOFF,
    Interaction_FDR =
      !is.na(padj_interaction) & padj_interaction < FDR_CUTOFF
  )

rewiring_pattern_levels <- c(
  "Pio response preserved",
  "Pio response potentiated",
  "Pio response attenuated",
  "Direction reversed in siSTAT1",
  "Pio response gained in siSTAT1",
  "Pio response lost/attenuated",
  "Weak/mixed",
  "Not detected in all contrasts"
)

rewiring_pattern_colours <- c(
  "Pio response preserved" = "#4D9221",
  "Pio response potentiated" = "#D73027",
  "Pio response attenuated" = "#4575B4",
  "Direction reversed in siSTAT1" = "#C51B7D",
  "Pio response gained in siSTAT1" = "#E08214",
  "Pio response lost/attenuated" = "#762A83",
  "Weak/mixed" = "#BDBDBD",
  "Not detected in all contrasts" = "#F0F0F0"
)

lipid_rewiring_wide <- lipid_rewiring_wide %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = rev(lipid_priority)),
    GeneGroup = factor(GeneGroup, levels = names(lipid_gene_groups)),
    RewiringPattern = factor(RewiringPattern, levels = rewiring_pattern_levels)
  )

p_rewiring_pattern <- ggplot(
  lipid_rewiring_wide,
  aes(x = "Pattern", y = SYMBOL, fill = RewiringPattern)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(
    aes(label = if_else(Interaction_FDR, "*", if_else(Interaction_nominal, "\u2020", ""))),
    size = 2.4
  ) +
  facet_grid(rows = vars(GeneGroup), scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = rewiring_pattern_colours,
    drop = FALSE,
    name = "Descriptive Pio-response pattern"
  ) +
  labs(x = NULL, y = NULL) +
  theme_fig(7.5) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.y = grid::unit(1.2, "mm"),
    legend.position = "right"
  )

p_rewiring_combined <- p_lipid_map + p_rewiring_pattern +
  patchwork::plot_layout(widths = c(5.2, 1.4), guides = "collect") &
  theme(legend.position = "right")

save_plot(
  p_rewiring_combined,
  "Fig6_Pio_lipid_handling_rewiring_map",
  9.2, 10.0
)

pattern_counts <- lipid_rewiring_wide %>%
  count(RewiringPattern, name = "n_genes", .drop = FALSE)

p_pattern_counts <- ggplot(
  pattern_counts,
  aes(forcats::fct_reorder(RewiringPattern, n_genes), n_genes, fill = RewiringPattern)
) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = n_genes), hjust = -0.15, size = 2.6) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = rewiring_pattern_colours, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Genes") +
  theme_fig(8)

save_plot(
  p_pattern_counts,
  "Fig6_Pio_lipid_rewiring_pattern_counts",
  6.2, 4.2
)

############################################################
### 19. Fig4-style lipid forest plot with 95% confidence intervals
############################################################

forest_contrast_colours <- c(
  "Pio in ctrl" = COL_CTRL_PIO,
  "Pio in siSTAT1" = COL_KD_PIO,
  "siSTAT1 alone" = COL_KD,
  "STAT1 x Pio" = COL_REDOX
)

p_lipid_forest <- ggplot(
  lipid_map_df,
  aes(log2FoldChange, SYMBOL, colour = Contrast)
) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_errorbar(
    aes(xmin = lower95, xmax = upper95),
    orientation = "y", width = 0, linewidth = 0.35,
    position = position_dodge(width = 0.60), na.rm = TRUE
  ) +
  geom_point(
    size = 1.6,
    position = position_dodge(width = 0.60),
    na.rm = TRUE
  ) +
  facet_grid(
    rows = vars(GeneGroup),
    scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_colour_manual(values = forest_contrast_colours, name = NULL) +
  labs(x = "Shrunken log2 fold change with 95% CI", y = NULL) +
  theme_fig(7.5) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    legend.position = "top"
  )

save_plot(p_lipid_forest, "Fig6_Pio_lipid_handling_forest", 7.8, 11.0)

############################################################
### 20. Fig4-style highlighted volcanoes
############################################################

volcano_label_genes <- unique(c(
  isg_priority,
  lipid_priority,
  redox_priority
))

plot_volcano_curated <- function(res_df, title_text) {
  df <- res_df %>%
    mutate(
      y = safe_neglog10(pvalue, cap = 30),
      Direction = case_when(
        !is.na(padj) & padj < FDR_CUTOFF & log2FoldChange >= LFC_CUTOFF ~ "FDR up",
        !is.na(padj) & padj < FDR_CUTOFF & log2FoldChange <= -LFC_CUTOFF ~ "FDR down",
        TRUE ~ "Background"
      ),
      Priority = case_when(
        SYMBOL %in% isg_priority ~ "STAT1 / ISG",
        SYMBOL %in% lipid_priority ~ "Lipid handling",
        SYMBOL %in% redox_priority ~ "Redox / iron",
        TRUE ~ NA_character_
      )
    )

  background <- df %>% filter(Direction == "Background", is.na(Priority))
  sig_nonpriority <- df %>% filter(Direction != "Background", is.na(Priority))
  priority <- df %>% filter(!is.na(Priority))

  top_nonpriority <- sig_nonpriority %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(padj, pvalue) %>%
    slice_head(n = 5)

  label_df <- bind_rows(
    priority %>%
      filter(SYMBOL %in% volcano_label_genes) %>%
      arrange(padj, pvalue) %>%
      group_by(Priority) %>%
      slice_head(n = 10) %>%
      ungroup(),
    top_nonpriority
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  ggplot() +
    geom_point(
      data = background, aes(log2FoldChange, y),
      colour = COL_BG, size = 0.65, alpha = 0.55
    ) +
    geom_point(
      data = sig_nonpriority, aes(log2FoldChange, y, colour = Direction),
      size = 0.90, alpha = 0.75
    ) +
    geom_point(
      data = priority, aes(log2FoldChange, y, colour = Priority),
      size = 1.25, alpha = 0.90
    ) +
    geom_vline(
      xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
      linetype = "dashed", linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(NOMINAL_P_CUTOFF),
      linetype = "dashed", linewidth = 0.3
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(log2FoldChange, y, label = SYMBOL),
      size = 2.0, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.2,
      box.padding = 0.25
    ) +
    scale_colour_manual(
      values = c(
        "FDR up" = COL_UP,
        "FDR down" = COL_DOWN,
        "STAT1 / ISG" = COL_IFN,
        "Lipid handling" = COL_LIPID,
        "Redox / iron" = COL_REDOX
      ),
      name = NULL
    ) +
    labs(
      title = title_text,
      x = expression(log[2]~"fold change (shrunk)"),
      y = expression(-log[10]*"(nominal P)")
    ) +
    theme_fig(8) +
    theme(legend.position = "bottom")
}

volcano_specs <- tribble(
  ~contrast_name, ~title_text, ~file_stem,
  "Pio_in_Control", "Pioglitazone in control", "Fig6_Pio_volcano_Pio_in_control",
  "Pio_in_siSTAT1", "Pioglitazone in siSTAT1", "Fig6_Pio_volcano_Pio_in_siSTAT1",
  "siSTAT1_no_Pio", "siSTAT1 versus control, no pioglitazone", "Fig6_Pio_volcano_siSTAT1_alone",
  "STAT1_x_Pio_interaction", "STAT1 x pioglitazone interaction", "Fig6_Pio_volcano_interaction"
)

pwalk(volcano_specs, function(contrast_name, title_text, file_stem) {
  save_plot(
    plot_volcano_curated(contrast_tables[[contrast_name]], title_text),
    file_stem, 4.5, 4.0
  )
})

############################################################
### 21. Fig4-style sample-level priority-gene heatmap
############################################################

heatmap_groups <- list(
  "STAT1 / ISG" = isg_priority,
  "Lipid handling" = lipid_priority,
  "Redox / iron" = redox_priority
)
heatmap_gene_order <- unique(unlist(heatmap_groups, use.names = FALSE))
heatmap_group_map <- enframe(
  heatmap_groups,
  name = "GeneGroup",
  value = "SYMBOL"
) %>%
  unnest(SYMBOL) %>%
  distinct(SYMBOL, .keep_all = TRUE)

heatmap_annot <- gene_annot %>%
  filter(SYMBOL %in% heatmap_gene_order) %>%
  distinct(SYMBOL, .keep_all = TRUE)

heatmap_ids <- heatmap_annot$ENSEMBL[
  match(heatmap_gene_order, heatmap_annot$SYMBOL)
]
heatmap_ids <- heatmap_ids[!is.na(heatmap_ids)]
heatmap_ids <- heatmap_ids[heatmap_ids %in% rownames(vst_mat)]

hm <- vst_mat[heatmap_ids, , drop = FALSE]
rownames(hm) <- gene_annot$SYMBOL[match(rownames(hm), gene_annot$ENSEMBL)]
hm <- hm[!duplicated(rownames(hm)), , drop = FALSE]
hm <- hm[intersect(heatmap_gene_order, rownames(hm)), , drop = FALSE]

if (nrow(hm) < 2L) {
  stop("Fewer than two priority genes were found in the VST matrix; heatmap cannot be drawn.")
}

hm_z <- t(scale(t(hm)))
hm_z[!is.finite(hm_z)] <- 0
sample_order <- meta %>% arrange(Group, Replicate) %>% pull(Sample_name)
hm_z <- hm_z[, sample_order, drop = FALSE]

annotation_col_priority <- meta[sample_order, c("Group", "Replicate"), drop = FALSE]
annotation_row_priority <- data.frame(
  GeneGroup = heatmap_group_map$GeneGroup[
    match(rownames(hm_z), heatmap_group_map$SYMBOL)
  ],
  row.names = rownames(hm_z),
  check.names = FALSE
)

annotation_colors_priority <- list(
  Group = GROUP_COLORS,
  Replicate = c("R1" = "#1B9E77", "R2" = "#D95F02", "R3" = "#7570B3"),
  GeneGroup = c(
    "STAT1 / ISG" = COL_IFN,
    "Lipid handling" = COL_LIPID,
    "Redox / iron" = COL_REDOX
  )
)

priority_hm <- pheatmap(
  hm_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = as.data.frame(annotation_col_priority),
  annotation_row = as.data.frame(annotation_row_priority),
  annotation_colors = annotation_colors_priority,
  color = colorRampPalette(c(COL_DOWN, "white", COL_UP))(101),
  border_color = "white",
  fontsize = 6,
  fontsize_row = 5.5,
  fontsize_col = 6,
  angle_col = "45",
  main = "Priority genes: sample-level row z-scores",
  silent = TRUE
)

svglite::svglite(
  file.path(figdir, "Fig6_Pio_priority_gene_sample_heatmap.svg"),
  width = 7.2, height = 12.0, bg = "transparent"
)
grid::grid.newpage(); grid::grid.draw(priority_hm$gtable); dev.off()

pdf(file.path(figdir, "Fig6_Pio_priority_gene_sample_heatmap.pdf"),
    width = 7.2, height = 12.0)
grid::grid.newpage(); grid::grid.draw(priority_hm$gtable); dev.off()

png(file.path(figdir, "Fig6_Pio_priority_gene_sample_heatmap.png"),
    width = 7.2, height = 12.0, units = "in", res = 600)
grid::grid.newpage(); grid::grid.draw(priority_hm$gtable); dev.off()

############################################################
### 22. Fig4-style redox/iron effect map
############################################################

redox_map_df <- bind_rows(lapply(map_contrasts, function(ct) {
  contrast_tables[[ct]] %>%
    filter(SYMBOL %in% redox_priority, !is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(SYMBOL, log2FoldChange, pvalue, padj, Contrast = ct)
})) %>%
  complete(SYMBOL = redox_priority, Contrast = map_contrasts) %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = rev(redox_priority)),
    Contrast = factor(map_labels[Contrast], levels = unname(map_labels[map_contrasts])),
    mark = sig_mark(padj, pvalue),
    cell_label = if_else(
      is.na(log2FoldChange), "\u00D7",
      paste0(sprintf("%.2f", log2FoldChange), mark)
    )
  )

p_redox_map <- ggplot(
  redox_map_df,
  aes(Contrast, SYMBOL, fill = log2FoldChange)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 2.2) +
  scale_fill_gradient2(
    low = COL_DOWN, mid = "white", high = COL_UP,
    midpoint = 0, limits = c(-COLOR_LIMIT, COLOR_LIMIT),
    oob = scales::squish, na.value = "grey88",
    name = "Shrunken\nlog2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_fig(8) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

save_plot(p_redox_map, "Fig6_Pio_redox_gene_effect_map", 5.6, 4.4)

############################################################
### 23. Source data and output manifest for added Fig4-style panels
############################################################

rewiring_definitions <- tibble(
  Pattern = rewiring_pattern_levels,
  Definition = c(
    "Control and siSTAT1 show the same non-weak Pio-response direction without a clear magnitude shift.",
    "The same-direction Pio response is stronger in siSTAT1 and the interaction follows that direction.",
    "The same-direction Pio response is weaker in siSTAT1 and the interaction opposes the control response.",
    "Control and siSTAT1 Pio effects cross the descriptive direction threshold in opposite directions.",
    "The Pio effect is weak in control but non-weak in siSTAT1.",
    "The Pio effect is non-weak in control but weak in siSTAT1.",
    "Effects are weak or do not meet another descriptive category.",
    "At least one required contrast lacks an effect estimate."
  )
)

write_xlsx_safe(
  list(
    lipid_effect_long = lipid_map_df %>% mutate(
      SYMBOL = as.character(SYMBOL),
      GeneGroup = as.character(GeneGroup),
      Contrast = as.character(Contrast)
    ),
    lipid_rewiring = lipid_rewiring_wide %>% mutate(
      SYMBOL = as.character(SYMBOL),
      GeneGroup = as.character(GeneGroup),
      RewiringPattern = as.character(RewiringPattern)
    ),
    rewiring_definitions = rewiring_definitions,
    rewiring_pattern_counts = pattern_counts,
    redox_effects = redox_map_df %>% mutate(
      SYMBOL = as.character(SYMBOL),
      Contrast = as.character(Contrast)
    ),
    lipid_gene_definitions = lipid_group_map,
    priority_heatmap_matrix = as.data.frame(hm_z) %>% rownames_to_column("SYMBOL")
  ),
  file.path(tabledir, "Figure6_Pio_Fig4_style_source_data.xlsx")
)

manifest <- tribble(
  ~File_stem, ~Suggested_use,
  "Fig6_Pio_KD_validation_STAT1_ISG", "Validate Stat1/ISG suppression before interpreting STAT1 dependence",
  "Fig6_Pio_module_dotplot", "Main Fig4-style module summary",
  "Fig6_Pio_full_module_dotplot", "Full module-level analysis",
  "Fig6_Pio_lipid_gene_effect_map", "Grouped lipid-handling effects across pairwise and interaction contrasts",
  "Fig6_Pio_lipid_handling_rewiring_map", "Effect map plus descriptive Pio-response rewiring class",
  "Fig6_Pio_lipid_handling_forest", "Lipid-gene effects with 95% confidence intervals",
  "Fig6_Pio_lipid_rewiring_pattern_counts", "Counts of descriptive rewiring classes",
  "Fig6_Pio_volcano_Pio_in_control", "Curated volcano for Pio response in control",
  "Fig6_Pio_volcano_Pio_in_siSTAT1", "Curated volcano for Pio response in siSTAT1",
  "Fig6_Pio_volcano_siSTAT1_alone", "Curated volcano for siSTAT1 effect without Pio",
  "Fig6_Pio_volcano_interaction", "Curated volcano for STAT1 x Pio interaction",
  "Fig6_Pio_priority_gene_sample_heatmap", "Sample-level STAT1/ISG, lipid, and redox expression",
  "Fig6_Pio_redox_gene_effect_map", "Focused redox/iron effects"
)
write.csv(
  manifest,
  file.path(outdir, "Figure6_Pio_output_manifest.csv"),
  row.names = FALSE
)

############################
### 24. README / reproducibility
############################

capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

writeLines(c(
  "FIGURE 6 (lipid arm) - siSTAT1 x Pioglitazone RNA-seq (YY2BYD, 72 h)",
  "",
  "Interpretation framework:",
  "1. Panel 1 (KD validation) must be read FIRST. The interpretation memo notes",
  "   Stat1 mRNA is not reduced at 72 h. If Stat1/ISGs are not suppressed here,",
  "   relabel the arm 'siSTAT1-transfected' and withhold STAT1-dependence claims.",
  "2. Pioglitazone effect = CD36-dominant uptake/handling, not classical de novo",
  "   lipogenesis or full adipogenesis (Pparg, Adipoq, Plin1, Cidec stay low).",
  "3. STAT1 x Pio interaction terms are hypothesis-generating only (n=3; expect",
  "   limited FDR power). The rewiring classes are descriptive summaries and do",
  "   not replace the interaction coefficient, P value, or FDR.",
  "4. Fig4-style additions include broader microglial modules, grouped lipid maps,",
  "   a Pio-response rewiring map, 95% CI forest plot, curated volcanoes, a priority",
  "   gene heatmap, redox/iron map, and an output manifest.",
  "",
  "Statistical design:",
  paste0("Model A (pairwise): ", deparse(design_group)),
  paste0("Model B (interaction): ", deparse(design_factorial)),
  paste0("Replicate blocking used: ", use_block),
  paste0("Interaction coefficient: ", interaction_name),
  "",
  "Significance marks: * FDR < 0.05 ; dagger nominal P < 0.05 ; x = not estimable."
), con = file.path(outdir, "README_Figure6_Pio_analysis_logic.txt"))

message("============================================================")
message("Figure 6 siSTAT1 x Pioglitazone analysis complete.")
message("Output: ", outdir)
message("Genes retained: ", nrow(count_filt))
message("Replicate blocking used: ", use_block)
message("Interaction coefficient: ", interaction_name)
message("============================================================")
