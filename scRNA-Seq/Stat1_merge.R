# ============================================
# STAT1 KO CORTEX: MERGE, ANNOTATE, SUBSET MICROGLIA
# Samples: TKG9, TKH1, TKH4 (STAT1 KO Ctx Female 8 weeks)
# ============================================

library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(dplyr)

base_dir <- "/Users/au806094/Downloads/scRNA_reprocess/Stat1_merge"
dir.create(file.path(base_dir, "Figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "Robj"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "DEG_output"), showWarnings = FALSE, recursive = TRUE)

# ============================================
# STEP 1: LOAD SAMPLES
# ============================================
print("Loading STAT1 KO samples...")

load("/Users/au806094/Downloads/dataset/Stat1_Cortex_8W/TKG9_Preprocessing/Robj/TKG9_processed.Robj") 
load("/Users/au806094/Downloads/dataset/Stat1_Cortex_8W/TKH1_Preprocessing/Robj/TKH1_processed.Robj")  
load("/Users/au806094/Downloads/dataset/Stat1_Cortex_8W/TKH4_Preprocessing/Robj/TKH4_processed.Robj")  

print(paste("TKG9:", ncol(TKG9), "cells"))
print(paste("TKH1:", ncol(TKH1), "cells"))
print(paste("TKH4:", ncol(TKH4), "cells"))

# Add metadata
TKG9$Genotype <- "Stat1_KO"
TKH1$Genotype <- "Stat1_KO"
TKH4$Genotype <- "Stat1_KO"
TKG9$SampleID <- "TKG9"
TKH1$SampleID <- "TKH1"
TKH4$SampleID <- "TKH4"

# ============================================
# STEP 2: MERGE
# ============================================
print("Merging STAT1 KO samples...")

ko_merged <- merge(TKG9,
                   y = list(TKH1, TKH4),
                   add.cell.ids = c("TKG9", "TKH1", "TKH4"),
                   project = "STAT1KO_Cortex")

print(paste("Total cells after merge:", ncol(ko_merged)))
table(ko_merged$SampleID)

rm(TKG9, TKH1, TKH4); gc()

# ============================================
# STEP 3: QC VISUALIZATION
# ============================================
print(grep("mt|mito", colnames(ko_merged@meta.data), value = TRUE, ignore.case = TRUE))

pdf(file.path(base_dir, "Figures", "QC_by_sample.pdf"), width = 14, height = 6)
print(VlnPlot(ko_merged,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = "SampleID",
              pt.size = 0, ncol = 3))
dev.off()

# ============================================
# STEP 4: NORMALIZE + PCA
# ============================================
print("SCTransform...")
ko_merged <- SCTransform(ko_merged,
                         vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         verbose = FALSE)

print("Running PCA...")
ko_merged <- RunPCA(ko_merged, npcs = 50, verbose = FALSE)

pdf(file.path(base_dir, "Figures", "ElbowPlot.pdf"), width = 8, height = 6)
print(ElbowPlot(ko_merged, ndims = 50))
dev.off()

# ============================================
# STEP 5: CHECK BATCH EFFECT BEFORE DECIDING
# ============================================
# With 3 replicates, more likely to see batch effects. Check the plot.

use_harmony <- TRUE  # <-- TOGGLE AFTER CHECKING THE PLOT

ko_merged@reductions$umap <- NULL
ko_merged <- RunUMAP(ko_merged, reduction = "pca", dims = 1:30, verbose = FALSE)

pdf(file.path(base_dir, "Figures", "UMAP_before_integration.pdf"), width = 14, height = 5)
p1 <- DimPlot(ko_merged, reduction = "umap", group.by = "SampleID", pt.size = 0.1) +
  ggtitle("By Sample - No Integration")
if ("Phase" %in% colnames(ko_merged@meta.data)) {
  p2 <- DimPlot(ko_merged, reduction = "umap", group.by = "Phase", pt.size = 0.1) +
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
  ko_merged@reductions$harmony <- NULL
  ko_merged <- RunHarmony(ko_merged,
                          group.by.vars = "SampleID",
                          reduction.use = "pca",
                          dims.use = 1:30,
                          reduction.save = "harmony")
  
  ko_merged@reductions$umap <- NULL
  ko_merged <- RunUMAP(ko_merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
  
  pdf(file.path(base_dir, "Figures", "UMAP_after_harmony.pdf"), width = 14, height = 5)
  p1 <- DimPlot(ko_merged, reduction = "umap", group.by = "SampleID", pt.size = 0.1) +
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

ko_merged <- FindNeighbors(ko_merged, reduction = red_to_use, dims = 1:30, verbose = FALSE)
ko_merged <- FindClusters(ko_merged, resolution = 0.2, verbose = FALSE)

pdf(file.path(base_dir, "Figures", "UMAP_clusters.pdf"), width = 10, height = 8)
print(DimPlot(ko_merged, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend())
dev.off()

pdf(file.path(base_dir, "Figures", "UMAP_split_by_sample.pdf"), width = 14, height = 5)
print(DimPlot(ko_merged, reduction = "umap", split.by = "SampleID", pt.size = 0.1))
dev.off()

# ============================================
# STEP 8: FIND MARKERS (using RNA assay)
# ============================================
print("Finding cluster markers...")

DefaultAssay(ko_merged) <- "RNA"
ko_merged <- JoinLayers(ko_merged)
ko_merged <- NormalizeData(ko_merged, verbose = FALSE)

all_markers <- FindAllMarkers(ko_merged,
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              only.pos = TRUE,
                              verbose = TRUE)

DefaultAssay(ko_merged) <- "SCT"

if (nrow(all_markers) > 0) {
  print(paste("Found", nrow(all_markers), "markers"))
  write.csv(all_markers, file.path(base_dir, "DEG_output", "cluster_markers.csv"))
  
  top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  write.csv(top10, file.path(base_dir, "DEG_output", "top10_markers.csv"))
  
  pdf(file.path(base_dir, "Figures", "Markers_heatmap.pdf"), width = 12, height = 10)
  print(DoHeatmap(ko_merged, features = top10$gene, size = 3))
  dev.off()
} else {
  warning("No markers found — check your clustering")
}

# Cortex markers quick check
cortex_markers <- c("Slc17a7", "Tbr1", "Gad1", "Gad2", "Pvalb", "Sst", "Vip",
                    "Aldh1l1", "Mbp", "Cx3cr1", "Pdgfra", "Pecam1")
available <- cortex_markers[cortex_markers %in% rownames(ko_merged)]

pdf(file.path(base_dir, "Figures", "Markers_featureplot.pdf"), width = 15, height = 12)
print(FeaturePlot(ko_merged, features = available, pt.size = 0.1, ncol = 4))
dev.off()

pdf(file.path(base_dir, "Figures", "Markers_dotplot.pdf"), width = 14, height = 6)
print(DotPlot(ko_merged, features = available) + RotatedAxis())
dev.off()

# ============================================
# STEP 9: ANNOTATE CLUSTERS
# ============================================
print(paste("Number of clusters:", length(unique(Idents(ko_merged)))))

# >> EDIT THIS MAP BASED ON YOUR MARKER RESULTS <<
celltype_map <- c(
  "0"  = "Astrocytes",
  "1"  = "Microglia",
  "2"  = "Endothelial",
  "3"  = "Oligodendrocytes",
  "4"  = "Neurons",
  "5"  = "VSMC",
  "6"  = "Astrocytes",
  "7"  = "OPCs",
  "8"  = "NPC",
  "9"  = "VLMC",
  "10" = "RBC",
  "11" = "Microglia",
  "12" = "BAM",
  "13" = "PVF",
  "14" = "Neurons",
  "15" = "Lymphocytes",
  "16" = "Astrocytes",
  "17" = "Ependymal",
  "18" = "Ependymal",
  "19" = "Pericytes"
)
# Add or remove lines to match your actual cluster count

# Safety check
unmapped <- setdiff(as.character(unique(Idents(ko_merged))), names(celltype_map))
if (length(unmapped) > 0) stop(paste("Unmapped clusters:", paste(unmapped, collapse = ", ")))

# Assign annotation
current_clusters <- as.character(Idents(ko_merged))
ko_merged@meta.data$CellType <- celltype_map[current_clusters]

table(ko_merged$CellType)

pdf(file.path(base_dir, "Figures", "UMAP_annotated.pdf"), width = 10, height = 8)
print(DimPlot(ko_merged, group.by = "CellType", label = TRUE, pt.size = 0.1) + NoLegend())
dev.off()

# ============================================
# STEP 10: SAVE FULL OBJECT + MICROGLIA SUBSET
# ============================================
print("Saving objects...")

saveRDS(ko_merged, file.path(base_dir, "Robj", "STAT1KO_merged_annotated.rds"))

ko_microglia <- subset(ko_merged, subset = CellType == "Microglia")
print(paste("STAT1 KO Microglia:", ncol(ko_microglia), "cells"))

if (ncol(ko_microglia) == 0) warning("No microglia found — check your CellType labels")

saveRDS(ko_microglia, file.path(base_dir, "Robj", "STAT1KO_microglia.rds"))

print("========================================")
print("STAT1 KO PIPELINE DONE")
print(paste("Total cells:", ncol(ko_merged)))
print(paste("Clusters:", length(unique(Idents(ko_merged)))))
print(paste("Microglia subset:", ncol(ko_microglia), "cells"))
print("========================================")