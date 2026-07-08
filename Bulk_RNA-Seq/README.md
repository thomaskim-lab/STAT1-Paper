**Bulk RNA-seq analysis**

This folder contains R scripts used for bulk RNA-seq and Plasmidsaurus DGE RNA-seq analyses in the STAT1 paper.

**Datasets**
- LPS: HMC3 LPS time-course bulk RNA-seq.
- 5SZCYD and BN4GZZ: HMC3 STAT1 siRNA experiments at 24 h.
- R4NCSL: HMC3 STAT1 siRNA with or without LPS at 48 h.
- 3BNGF5: Primary mouse microglia Stat1 siRNA experiment.
- BTKTYZ: Primary mouse microglia Stat1 siRNA with or without IFNγ rescue.

**Suggested script order**
- 00_bulk_setup.R
- 01_hmc3_lps_analysis.R
- 02_hmc3_stat1_sirna_5SZCYD_BN4GZZ.R
- 03_hmc3_stat1_sirna_lps_R4NCSL.R
- 04_pmg_stat1_sirna_3BNGF5.R
- 05_bulk_module_scores.R
- 06_pmg_stat1_ifng_rescue_BTKTYZ.R
- 07_
- 99_bulk_session_info.R

**Outputs**

Bulk RNA-seq outputs include normalized expression matrices, DESeq2 differential-expression tables, PCA plots, sample-distance heatmaps, volcano plots, curated module summaries, GO enrichment, Reactome enrichment, and MSigDB/fgsea outputs where applicable.

**Reproducibility**

Run 99_bulk_session_info.R after installing and loading the required packages. This creates session_info/bulk_rnaseq_sessionInfo.txt.
