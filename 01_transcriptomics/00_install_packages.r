# ==============================================================================
# 00_install_packages.R
# AD Transcriptomic Meta-Analysis — Package Installation
#
# Run this ONCE before any other script.
# Installs all CRAN and Bioconductor dependencies.
# Safe to re-run: installed packages are skipped automatically.
# ==============================================================================

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  "BiocManager",   # Bioconductor installer
  "tidyverse",     # dplyr, ggplot2, readr, stringr, purrr, tidyr
  "ggrepel",       # non-overlapping labels in ggplot2
  "scales",        # axis formatting helpers
  "Cairo",         # high-quality PDF/PNG output (cairo_pdf device)
  "pheatmap",      # fallback heatmap (not used in final figures)
  "RColorBrewer",  # colour palettes
  "gridExtra",     # multi-panel layouts
  "cowplot",       # publication-quality plot composition
  "here"           # project-relative file paths
)

new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if (length(new_cran)) {
  message("Installing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org",
                   dependencies = TRUE, quiet = TRUE)
} else {
  message("All CRAN packages already installed.")
}

# ── Bioconductor packages ──────────────────────────────────────────────────────
bioc_pkgs <- c(
  # Data access
  "GEOquery",          # download datasets directly from NCBI GEO
  "Biobase",           # ExpressionSet class, pData(), exprs()

  # Microarray preprocessing
  "limma",             # differential expression for microarray and RNA-seq

  # RNA-seq preprocessing
  "DESeq2",            # differential expression for count data
  "edgeR",             # TMM normalisation, filterByExpr

  # Annotation
  "org.Hs.eg.db",      # human gene annotation (symbol ↔ Entrez ↔ Ensembl)
  "AnnotationDbi",     # annotation database interface
  "biomaRt",           # Ensembl ID conversion (used for GSE118553)

  # Functional enrichment
  "clusterProfiler",   # GO and KEGG enrichment
  "enrichplot",        # dotplot, barplot for enrichment results
  "DOSE",              # disease ontology (loaded by clusterProfiler)
  "ReactomePA",        # optional Reactome pathway enrichment

  # Heatmap
  "ComplexHeatmap",    # publication-quality annotated heatmaps
  "circlize",          # colorRamp2() for heatmap colour scales

  # Meta-analysis support
  "metafor"            # formal meta-analysis (used for heterogeneity checks)
)

new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[,"Package"]]
if (length(new_bioc)) {
  message("Installing Bioconductor packages: ", paste(new_bioc, collapse = ", "))
  BiocManager::install(new_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
} else {
  message("All Bioconductor packages already installed.")
}

# ── Version check ─────────────────────────────────────────────────────────────
message("\n── Installed versions ───────────────────────────────────────────────────")
key_pkgs <- c("limma","DESeq2","clusterProfiler","ComplexHeatmap",
              "GEOquery","org.Hs.eg.db","ggplot2","ggrepel","dplyr")
for (p in key_pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  message(sprintf("  %-22s %s", p, v))
}
message("─────────────────────────────────────────────────────────────────────────")
message("Installation complete. Proceed to 01_download_data.R")
