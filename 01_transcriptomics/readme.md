# Phase I: Transcriptomic Target Validation

This module executes the cross-dataset weighted meta-analysis of Alzheimer's disease brain tissue and cellular models to identify robust transcriptional signatures and validate primary molecular targets (MAPT and GSK3B).

## Datasets
* **GSE138260**: Agilent microarray, post-mortem temporal cortex (17 AD, 19 controls).
* **GSE118553**: Illumina RNA-seq, iPSC-derived neurons (familial AD mutations and isogenic controls, n = 207).

## Execution Pipeline
The analysis is modularized into 11 R scripts. Execute them sequentially:

### 1. Environment & Data Acquisition
* `00_install_packages.r`: Installs required CRAN and Bioconductor dependencies.
* `01_download_data.r`: Programmatically retrieves raw expression matrices and phenotype data from NCBI GEO.

### 2. Preprocessing & Meta-Analysis
* `02_preprocess_GSE138260.r`: Executes quality control, low-expression probe filtering, and limma-based differential expression.
* `03_preprocess_GSE118553.r`: Applies VST-normalization and DESeq2 differential expression.
* `04_meta_analysis.r`: Performs weighted fixed-effects meta-analysis using Stouffer's z-score method to isolate directionally consistent DEGs.

### 3. Visualization
* `05_volcano_v4.r`: Generates bilateral volcano plots for the 554 high-confidence DEGs.
* `06_scatter_v3.r`: Plots cross-dataset log₂ fold change consistency.
* `07a_go_dotplot.r`: Visualizes GO Biological Process enrichment.
* `07b_kegg_dotplot.r`: Visualizes KEGG Pathway enrichment.
* `08_lollipop_v3.r`: Generates error-bar profiles for the 15-gene targeted panel.
* `09_heatmap_v2.r`: Renders the z-scored expression heatmap for the top 40 reproducible DEGs.
