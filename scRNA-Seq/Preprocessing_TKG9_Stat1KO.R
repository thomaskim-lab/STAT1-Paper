# Single-cell RNA-seq Preprocessing Pipeline
# Sample: TKG9 (Stat1 KO, Female, 8 Week, Cortex)

# ============================================
# STEP 0: INSTALL AND LOAD PACKAGES
# ============================================

# Install required packages if not already installed
# Run these lines if this is your first time
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Seurat and dependencies
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

# Install other CRAN packages
packages_to_install <- c("cowplot", "dplyr", "Matrix", "RColorBrewer", "harmony")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
install.packages("hdf5r")


library(cowplot)
library(dplyr)
library(Matrix)
library(Seurat)
library(RColorBrewer)
library(harmony)
library(hdf5r)  # Required for reading .h5 files

# ============================================
# SETUP WORKING DIRECTORY AND CREATE FOLDERS
# ============================================

# Set base working directory
base_dir <- "/Users/au806094/Downloads/dataset"
setwd(base_dir)


dir.create(file.path(base_dir, "QC_output"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "DEG_output"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "Robj"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "Figures"), showWarnings = FALSE, recursive = TRUE)

# ============================================
# STEP 1: LOAD DATA
# ============================================
# Read10X_h5 reads the filtered feature barcode matrix from 10X Genomics
print("Loading TKG9 data...")
TKG9 <- Read10X_h5("/Users/au806094/Downloads/dataset/raw_from_genomedk/8WK/TKG9/Stat1KO_1_filtered_feature_bc_matrix.h5")

# Add unique prefix to avoid barcode clashes when merging later
# Each cell barcode gets "TKG9_" prefix
colnames(TKG9) = paste0("TKG9_", colnames(TKG9))

# ============================================
# STEP 2: CREATE SEURAT OBJECT WITH INITIAL FILTERING
# ============================================
# CreateSeuratObject creates the main data structure for analysis
# min.cells = 5: Keep genes expressed in at least 5 cells
# min.features = 1000: Keep cells with at least 1000 genes detected
print("Creating Seurat object with initial filtering...")
TKG9 <- CreateSeuratObject(counts = TKG9, 
                           project = "Stat1_KO",
                           min.cells = 5, 
                           min.features = 800)

# Print initial statistics
print(paste("Initial cells:", ncol(TKG9)))
print(paste("Initial genes:", nrow(TKG9)))

# Additional filtering based on UMI counts
# This removes low-quality cells with very few transcripts
#1000 if depth is higher
#UMI #1000/2000/3000
#Filter based on these units.  #not by mitochondria and ribosomal -> aging sample
TKG9 <- subset(TKG9, subset = nCount_RNA > 2000)
print(paste("Cells after UMI filter (>2000):", ncol(TKG9)))

# ============================================
# STEP 3: ADD METADATA
# ============================================
print("Adding metadata...")

# Check what identities currently exist
print("Current identity classes:")
print(table(Idents(TKG9)))
print("Original identity:")
print(table(TKG9$orig.ident))

# Method 2 (RECOMMENDED): Add metadata directly without RenameIdents
# This is cleaner and more straightforward
TKG9 <- AddMetaData(TKG9, "TKG9_Stat1KO", col.name = "Sample")
TKG9 <- AddMetaData(TKG9, "Stat1_KO", col.name = "Genotype")
TKG9 <- AddMetaData(TKG9, "8_Week", col.name = "Age")
TKG9 <- AddMetaData(TKG9, "Female", col.name = "Sex")
TKG9 <- AddMetaData(TKG9, "Cortex", col.name = "Region")
TKG9 <- AddMetaData(TKG9, "TKG9", col.name = "SampleID")

# Check metadata structure
print("Metadata structure:")
head(TKG9@meta.data)
print("Metadata columns available:")
print(colnames(TKG9@meta.data))
# ============================================
# STEP 4: CALCULATE QC METRICS
# ============================================
print("Calculating QC metrics...")
# Calculate percentage of mitochondrial genes (indicator of cell stress/death)
TKG9[["percent.mt"]] <- PercentageFeatureSet(TKG9, pattern = "^mt-")

# Calculate percentage of ribosomal proteins (can indicate cell type)
TKG9[["percent.RPS"]] <- PercentageFeatureSet(TKG9, pattern = "^Rps")
TKG9[["percent.RPL"]] <- PercentageFeatureSet(TKG9, pattern = "^Rpl")

# For sex check - X-linked genes
TKG9[["percent.sex"]] <- PercentageFeatureSet(TKG9, 
                                              features = c("Xist", "Tsix"))

# Print QC summary before filtering
print("QC metrics summary before filtering:")
print(summary(TKG9@meta.data[c("nCount_RNA", "nFeature_RNA", "percent.mt")]))

# ============================================
# STEP 5: FILTER BASED ON QC METRICS
# ============================================
print("Applying QC filters...")
# Remove cells with high mitochondrial content (likely dying cells)
TKG9 <- subset(TKG9, subset = percent.mt < 20)  # More stringent for brain tissue
print(paste("Cells after mitochondrial filter (<20%):", ncol(TKG9)))

# Remove cells with very high ribosomal content
TKG9 <- subset(TKG9, subset = percent.RPS < 25)
print(paste("Cells after ribosomal filter (<25%):", ncol(TKG9)))

# ============================================
# STEP 6: VISUALIZE QC METRICS
# ============================================
print("Creating QC plots...")
# Create violin plots to assess data quality
pdf(file.path(base_dir, "Figures", "TKG9_QC_violin.pdf"), width = 12, height = 8)
VlnPlot(TKG9, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                     "percent.RPS", "percent.RPL", "percent.sex"), 
        pt.size = 0,
        ncol = 3)
dev.off()

# Also save as PNG for easy viewing
png(file.path(base_dir, "Figures", "TKG9_QC_violin.png"), width = 1200, height = 800)
VlnPlot(TKG9, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                     "percent.RPS", "percent.RPL", "percent.sex"), 
        pt.size = 0,
        ncol = 3)
dev.off()

# Save QC metrics to file for records
QUALITY <- TKG9@meta.data
write.csv(QUALITY, file = file.path(base_dir, "QC_output", "TKG9_QC_metrics.csv"))
print(paste("QC metrics saved to:", file.path(base_dir, "QC_output", "TKG9_QC_metrics.csv")))

# ============================================
# STEP 7: NORMALIZATION WITH SCTRANSFORM
# ============================================
print("Running SCTransform normalization...")
# SCTransform performs normalization, variance stabilization, and scaling
# It also regresses out technical variables
TKG9 <- SCTransform(TKG9, 
                    vars.to.regress = c("nCount_RNA", "nFeature_RNA"),
                    verbose = FALSE)

# ============================================
# STEP 8: DIMENSIONALITY REDUCTION - PCA
# ============================================
print("Running PCA...")
# Run Principal Component Analysis
TKG9 <- RunPCA(TKG9, npcs = 50, verbose = FALSE)

# Visualize how much variance each PC explains
pdf(file.path(base_dir, "Figures", "TKG9_ElbowPlot.pdf"), width = 8, height = 6)
ElbowPlot(TKG9, ndims = 50, reduction = "pca")
dev.off()
# Look for the "elbow" where variance drops off

# ============================================
# STEP 9: CELL CYCLE SCORING
# ============================================
print("Performing cell cycle scoring...")
# Load cell cycle genes and convert from human to mouse
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Alternative: Use Seurat's built-in mouse cell cycle genes
# These are commonly used mouse S and G2/M phase markers
m.s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", 
               "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells",
               "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2",
               "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6",
               "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1",
               "Chaf1b", "Brip1", "E2f8")

m.g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", 
                 "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3",
                 "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", 
                 "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3",
                 "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", 
                 "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", 
                 "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3",
                 "Gas2l3", "Cbx5", "Cenpa")

TKG9 <- CellCycleScoring(TKG9, 
                         s.features = m.s.genes,
                         g2m.features = m.g2m.genes, 
                         set.ident = FALSE)

# ============================================
# STEP 10: RUN UMAP FOR VISUALIZATION
# ============================================
print("Running UMAP...")
# UMAP provides 2D visualization of cells
# Using first 20 PCs usually captures most variation
TKG9 <- RunUMAP(TKG9, 
                dims = 1:20,
                n.neighbors = 30L,  # Increased for larger datasets
                min.dist = 0.3,     # Standard for most analyses
                spread = 1)

# Visualize by cell cycle phase
pdf(file.path(base_dir, "Figures", "TKG9_UMAP_cellcycle.pdf"), width = 8, height = 6)
DimPlot(TKG9, reduction = "umap", group.by = "Phase", pt.size = 0.5)
dev.off()

# ============================================
# STEP 11: CLUSTERING
# ============================================
print("Finding clusters...")
# Find cell neighbors and clusters
TKG9 <- FindNeighbors(TKG9, dims = 1:20)
TKG9 <- FindClusters(TKG9, resolution = 0.5)  # Start with 0.5 for cortex

# Visualize clusters
pdf(file.path(base_dir, "Figures", "TKG9_UMAP_clusters.pdf"), width = 8, height = 6)
DimPlot(TKG9, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# ============================================
# STEP 12: CHECK CORTEX-SPECIFIC MARKERS
# ============================================
print("Checking cortical markers...")
# These are cortical layer and cell type markers
cortex_markers <- c(
  
  # Excitatory neurons
  "Slc17a7", # vGlut1
  "Satb2",   # Callosal neurons
  "Tbr1",    # Layer 6
  
  # Inhibitory neurons
  "Gad1", "Gad2",  # GABAergic
  "Pvalb",   # Parvalbumin
  "Sst",     # Somatostatin
  "Vip",     # VIP
  
  # Glia
  "Aldh1l1", # Astrocytes
  "Mbp",     # Oligodendrocytes  
  "Cx3cr1",  # Microglia
  "Pdgfra",  # OPC
  
  # Other
  "Pecam1"   # Endothelial
)

# Check which markers are present in the dataset
available_markers <- cortex_markers[cortex_markers %in% rownames(TKG9)]
print(paste("Available markers:", length(available_markers), "out of", length(cortex_markers)))

# Create feature plots for available markers
if(length(available_markers) > 0) {
  pdf(file.path(base_dir, "Figures", "TKG9_cortex_markers.pdf"), width = 12, height = 10)
  print(FeaturePlot(TKG9, 
                    features = available_markers[1:min(15, length(available_markers))],
                    cols = c("lightgrey", "red"),
                    pt.size = 0.1,
                    ncol = 3))
  dev.off()
}

# ============================================
# STEP 13: FIND INITIAL MARKERS
# ============================================
print("Finding cluster markers...")
# Identify genes that define each cluster
markers <- FindAllMarkers(TKG9, 
                          test.use = "wilcox",
                          logfc.threshold = 0.5,
                          min.pct = 0.25,
                          verbose = FALSE)

# Save markers
write.csv(markers, file = file.path(base_dir, "DEG_output", "TKG9_cluster_markers.csv"))
print(paste("Cluster markers saved to:", file.path(base_dir, "DEG_output", "TKG9_cluster_markers.csv")))

# Save top 10 markers per cluster
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = file.path(base_dir, "DEG_output", "TKG9_top10_markers.csv"))

# ============================================
# STEP 14: SAVE PROCESSED OBJECT
# ============================================
print("Saving processed Seurat object...")
save(TKG9, file = file.path(base_dir, "Robj", "TKG9_processed.Robj"))

# ============================================
# FINAL SUMMARY
# ============================================
print("========================================")
print("TKG9 preprocessing complete!")
print("========================================")
print(paste("Number of cells after QC:", ncol(TKG9)))
print(paste("Number of genes:", nrow(TKG9)))
print(paste("Number of clusters found:", length(unique(Idents(TKG9)))))
print("========================================")
print("Output files created:")
print(paste("- QC metrics:", file.path(base_dir, "QC_output", "TKG9_QC_metrics.csv")))
print(paste("- Cluster markers:", file.path(base_dir, "DEG_output", "TKG9_cluster_markers.csv")))
print(paste("- Seurat object:", file.path(base_dir, "Robj", "TKG9_processed.Robj")))
print(paste("- Figures saved in:", file.path(base_dir, "Figures")))
print("========================================")