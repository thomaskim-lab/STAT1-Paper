# ============================================
# WT CORTEX: MERGE, ANNOTATE, SUBSET MICROGLIA
# Samples: TKH8 & TKH9 (WT Ctx Female 8 weeks)
# ============================================

library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(dplyr)

base_dir <- "/Users/au806094/Downloads/dataset/WT_Merge"
dir.create(file.path(base_dir, "Figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "Robj"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "DEG_output"), showWarnings = FALSE, recursive = TRUE)

# ============================================
# STEP 1: LOAD SAMPLES
# ============================================
print("Loading WT samples...")

load("/Users/au806094/Downloads/dataset/Stat1_Cortex_8W/TKH8_Preprocessing/Robj/TKH8_processed.Robj")
load("/Users/au806094/Downloads/dataset/Stat1_Cortex_8W/TKH9_Preprocessing/Robj/TKH9_processed.Robj")

print(paste("TKH8:", ncol(TKH8), "cells"))
print(paste("TKH9:", ncol(TKH9), "cells"))

# Add metadata
TKH8$Genotype <- "WT"
TKH9$Genotype <- "WT"
TKH8$SampleID <- "TKH8"
TKH9$SampleID <- "TKH9"

# ============================================
# STEP 2: MERGE
# ============================================
print("Merging WT samples...")

wt_merged <- merge(TKH8,
                   y = list(TKH9),
                   add.cell.ids = c("TKH8", "TKH9"),
                   project = "WT_Cortex")

print(paste("Total cells after merge:", ncol(wt_merged)))
table(wt_merged$SampleID)

rm(TKH8, TKH9); gc()

# ============================================
# STEP 3: QC VISUALIZATION
# ============================================
# Check what your mito column is called
print(grep("mt|mito", colnames(wt_merged@meta.data), value = TRUE, ignore.case = TRUE))
# If it prints "percent.mito" instead of "percent.mt", change all references below

pdf(file.path(base_dir, "Figures", "QC_by_sample.pdf"), width = 12, height = 6)
print(VlnPlot(wt_merged,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = "SampleID",
              pt.size = 0, ncol = 3))
dev.off()

# ============================================
# STEP 4: NORMALIZE + PCA
# ============================================
print("SCTransform...")
wt_merged <- SCTransform(wt_merged,
                         vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         verbose = FALSE)

print("Running PCA...")
wt_merged <- RunPCA(wt_merged, npcs = 50, verbose = FALSE)

pdf(file.path(base_dir, "Figures", "ElbowPlot.pdf"), width = 8, height = 6)
print(ElbowPlot(wt_merged, ndims = 50))
dev.off()

# ============================================
# STEP 5: CHECK BATCH EFFECT BEFORE DECIDING
# ============================================
# Look at this plot: if TKH8 and TKH9 overlap well, skip Harmony.
# If they separate, set use_harmony <- TRUE and re-run from here.

use_harmony <- TRUE  # <-- TOGGLE AFTER CHECKING THE PLOT

wt_merged@reductions$umap <- NULL
wt_merged <- RunUMAP(wt_merged, reduction = "pca", dims = 1:30, verbose = FALSE)

pdf(file.path(base_dir, "Figures", "UMAP_before_integration.pdf"), width = 12, height = 5)
p1 <- DimPlot(wt_merged, reduction = "umap", group.by = "SampleID", pt.size = 0.1) +
  ggtitle("By Sample - No Integration")
# Only plot Phase if cell cycle scoring was done in preprocessing
if ("Phase" %in% colnames(wt_merged@meta.data)) {
  p2 <- DimPlot(wt_merged, reduction = "umap", group.by = "Phase", pt.size = 0.1) +
    ggtitle("By Cell Cycle")
  print(plot_grid(p1, p2, ncol = 2))
} else {
  print(p1)
}
dev.off()

# ============================================
# STEP 6: HARMONY (OPTIONAL)
# ============================================
if (use_harmony) {
  print("Running Harmony on SampleID...")
  wt_merged@reductions$harmony <- NULL
  wt_merged <- RunHarmony(wt_merged,
                          group.by.vars = "SampleID",
                          reduction.use = "pca",
                          dims.use = 1:30,
                          reduction.save = "harmony")
  
  wt_merged@reductions$umap <- NULL
  wt_merged <- RunUMAP(wt_merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
  
  pdf(file.path(base_dir, "Figures", "UMAP_after_harmony.pdf"), width = 12, height = 5)
  p1 <- DimPlot(wt_merged, reduction = "umap", group.by = "SampleID", pt.size = 0.1) +
    ggtitle("By Sample - After Harmony")
  print(p1)
  dev.off()
  
  red_to_use <- "harmony"
} else {
  red_to_use <- "pca"
}

# ============================================
# STEP 7: CLUSTERING
# ============================================
print("Clustering...")

wt_merged <- FindNeighbors(wt_merged, reduction = red_to_use, dims = 1:30, verbose = FALSE)
wt_merged <- FindClusters(wt_merged, resolution = 0.3, verbose = FALSE)

# No UMAP re-run needed — already computed in Step 5/6

pdf(file.path(base_dir, "Figures", "UMAP_clusters.pdf"), width = 10, height = 8)
print(DimPlot(wt_merged, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend())
dev.off()

pdf(file.path(base_dir, "Figures", "UMAP_split_by_sample.pdf"), width = 12, height = 5)
print(DimPlot(wt_merged, reduction = "umap", split.by = "SampleID", pt.size = 0.1))
dev.off()

# ============================================
# STEP 8: FIND MARKERS (using RNA assay)
# ============================================
print("Finding cluster markers...")

# Standard approach: use log-normalized RNA for marker detection
# SCT is used for clustering/UMAP, RNA for DE
DefaultAssay(wt_merged) <- "RNA"
wt_merged <- JoinLayers(wt_merged)
wt_merged <- NormalizeData(wt_merged, verbose = FALSE)

all_markers <- FindAllMarkers(wt_merged,
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              only.pos = TRUE,
                              verbose = TRUE)

# Switch back
DefaultAssay(wt_merged) <- "SCT"

if (nrow(all_markers) > 0) {
  print(paste("Found", nrow(all_markers), "markers"))
  write.csv(all_markers, file.path(base_dir, "DEG_output", "cluster_markers.csv"))
  
  top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  write.csv(top10, file.path(base_dir, "DEG_output", "top10_markers.csv"))
  
  pdf(file.path(base_dir, "Figures", "Markers_heatmap.pdf"), width = 12, height = 10)
  print(DoHeatmap(wt_merged, features = top10$gene, size = 3))
  dev.off()
} else {
  warning("No markers found — check your clustering")
}

# Cortex markers quick check
cortex_markers <- c("Slc17a7", "Tbr1", "Gad1", "Gad2", "Pvalb", "Sst", "Vip",
                    "Aldh1l1", "Mbp", "Cx3cr1", "Pdgfra", "Pecam1")
available <- cortex_markers[cortex_markers %in% rownames(wt_merged)]

pdf(file.path(base_dir, "Figures", "Markers_featureplot.pdf"), width = 15, height = 12)
print(FeaturePlot(wt_merged, features = available, pt.size = 0.1, ncol = 4))
dev.off()

pdf(file.path(base_dir, "Figures", "Markers_dotplot.pdf"), width = 14, height = 6)
print(DotPlot(wt_merged, features = available) + RotatedAxis())
dev.off()

# ============================================
# STEP 9: ANNOTATE CLUSTERS
# ============================================
print(paste("Number of clusters:", length(unique(Idents(wt_merged)))))

celltype_map <- c(
  "0"  = "Astrocytes",
  "1"  = "Microglia",
  "2"  = "Endothelial",
  "3"  = "Microglia",
  "4"  = "Neurons",
  "5"  = "Oligodendrocytes",
  "6"  = "Pericytes",
  "7"  = "Astrocytes",
  "8"  = "OPCs",
  "9"  = "VSMC",
  "10" = "PVF",
  "11" = "VLMC",
  "12" = "Astrocytes",
  "13" = "RBC",
  "14" = "BAM",
  "15" = "Neurons",
  "16" = "Endothelial",
  "17" = "Lymphocytes",
  "18" = "BAM",
  "19" = "Ependymal",
  "20" = "NK",
  "21" = "Ependymal"
)

# Safety check
unmapped <- setdiff(as.character(unique(Idents(wt_merged))), names(celltype_map))
if (length(unmapped) > 0) stop(paste("Unmapped clusters:", paste(unmapped, collapse = ", ")))

# Assign annotation
current_clusters <- as.character(Idents(wt_merged))
wt_merged@meta.data$CellType <- celltype_map[current_clusters]

table(wt_merged$CellType)

pdf(file.path(base_dir, "Figures", "UMAP_annotated.pdf"), width = 10, height = 8)
print(DimPlot(wt_merged, group.by = "CellType", label = TRUE, pt.size = 0.1) + NoLegend())
dev.off()

# ============================================
# STEP 10: SAVE FULL OBJECT + MICROGLIA SUBSET
# ============================================
print("Saving objects...")

saveRDS(wt_merged, file.path(base_dir, "Robj", "WT_merged_annotated.rds"))

# Subset microglia — make sure the label matches exactly what you wrote above
wt_microglia <- subset(wt_merged, subset = CellType == "Microglia")
print(paste("WT Microglia:", ncol(wt_microglia), "cells"))

if (ncol(wt_microglia) == 0) warning("No microglia found — check your CellType labels")

saveRDS(wt_microglia, file.path(base_dir, "Robj", "WT_microglia.rds"))

print("========================================")
print("WT PIPELINE DONE")
print(paste("Total cells:", ncol(wt_merged)))
print(paste("Clusters:", length(unique(Idents(wt_merged)))))
print(paste("Microglia subset:", ncol(wt_microglia), "cells"))
print("========================================")