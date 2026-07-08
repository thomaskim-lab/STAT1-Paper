**STAT1-Paper**

Analysis repository for the STAT1 paper, including bulk RNA-seq, single-cell RNA-seq, single-cell ATAC-seq, flow cytometry, and SABER-FISH analysis workflows.

**Repository structure**
- Bulk_RNA-Seq/: R scripts for Plasmidsaurus DGE RNA-seq and HMC3 bulk RNA-seq analyses.
- scRNA-Seq/: R scripts for 10x Genomics scRNA-seq import, quality control, integration, clustering, microglia subsetting, differential expression, and module scoring.
- scATAC-Seq/: R scripts for 10x Genomics scATAC-seq import, quality control, chromatin accessibility, promoter/distal accessibility, and motif analysis.
- metadata/: Sample annotation tables used by the analysis scripts.
- gene_sets/: Curated gene modules used for focused pathway and module-score analyses.
- results/: Processed analysis outputs, including differential-expression tables, enrichment results, and module-score summaries.
- session_info/: R session information files for reproducibility.

**Bulk RNA-seq datasets**
Dataset/order ID	Cell type	Species	Samples	Description
LPS	HMC3	Human	X	LPS time-course analysis
5SZCYD	HMC3	Human	6	STAT1 siRNA #1 and #3, 24 h
BN4GZZ	HMC3	Human	6	STAT1 siRNA #1 and #3, 24 h
R4NCSL	HMC3	Human	12	STAT1 siRNA #1 and #3 with or without LPS, 48 h
3BNGF5	Primary microglia	Mouse	8	Untreated, Lipofectamine-only, Stat1 siRNA #1 and #3, 48-72 h
BTKTYZ	Primary microglia	Mouse	12	Stat1 siRNA with or without IFNγ rescue


**Reproducibility**

Scripts were run in R. Package versions are recorded in the session_info/ directory. Where available, renv.lock records the package environment.

**Data availability**

Raw sequencing data will be deposited in [GEO/SRA accession]. Processed count matrices, sample metadata, differential-expression outputs, gene modules, and analysis scripts are provided in this repository.

**Notes**

Large raw sequencing files are not stored in this GitHub repository. Sensitive information, including login credentials, billing information, animal requisition forms, and internal order documents, has been excluded.
