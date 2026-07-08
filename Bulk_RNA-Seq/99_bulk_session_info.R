# 99_bulk_session_info.R
# Records R version and package versions for bulk RNA-seq analysis.

dir.create("../session_info", showWarnings = FALSE, recursive = TRUE)

# Core data handling
library(tidyverse)
library(readr)
library(readxl)
library(openxlsx)

# Differential expression
library(DESeq2)

# Annotation
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Enrichment and gene sets
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(fgsea)
library(enrichplot)

# Plotting
library(pheatmap)
library(ggrepel)
library(patchwork)
library(VennDiagram)
library(ggpubr)

# Optional but useful for cleaner session reporting
if (!requireNamespace("sessioninfo", quietly = TRUE)) {
  install.packages("sessioninfo")
}
library(sessioninfo)

# Base R session info
sink("../session_info/bulk_rnaseq_sessionInfo.txt")
cat("Bulk RNA-seq session information\n")
cat("Generated on: ", as.character(Sys.time()), "\n\n")
print(sessionInfo())
cat("\n\nDetailed sessioninfo::session_info()\n\n")
print(sessioninfo::session_info())
sink()