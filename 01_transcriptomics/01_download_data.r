# ==============================================================================
# 01_download_data.R
# AD Transcriptomic Meta-Analysis — Data Download from NCBI GEO
#
# Downloads and caches:
#   GSE138260 — Agilent microarray, post-mortem temporal cortex (n=36)
#   GSE118553 — RNA-seq, iPSC-derived neurons (n=207)
#
# Output: data/raw/  (relative to project root)
# Runtime: ~5–15 min depending on network speed
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(here)
})

# ── Directory setup ────────────────────────────────────────────────────────────
RAW_DIR <- here("data", "raw")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# GEOquery options: increase timeout for large supplementary files
options(timeout = 600)

# ── Helper: download and cache a GEO series ────────────────────────────────────
get_geo_series <- function(gse_id, destdir) {
  cache_file <- file.path(destdir, paste0(gse_id, "_eset.RData"))
  if (file.exists(cache_file)) {
    message(gse_id, ": loading from cache.")
    load(cache_file)
    return(eset)
  }
  message(gse_id, ": downloading from GEO...")
  gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = destdir)
  # getGEO may return a list when multiple platforms exist; take the first
  if (is.list(gse)) gse <- gse[[1]]
  eset <- gse
  save(eset, file = cache_file)
  message(gse_id, ": saved to ", cache_file)
  return(eset)
}

# ── Download GSE138260 (microarray) ───────────────────────────────────────────
message("\n── GSE138260: Agilent microarray, post-mortem temporal cortex ───────────")
eset_138260 <- get_geo_series("GSE138260", RAW_DIR)

# Quick sanity check
pd_138260 <- pData(eset_138260)
message("  Samples   : ", nrow(pd_138260))
message("  Disease col: ", names(pd_138260)[grep("disease", names(pd_138260),
                                                  ignore.case = TRUE)[1]])
print(table(pd_138260$`disease state:ch1`))

# ── Download GSE118553 (RNA-seq count matrix) ─────────────────────────────────
message("\n── GSE118553: RNA-seq, iPSC-derived neurons ─────────────────────────────")
eset_118553 <- get_geo_series("GSE118553", RAW_DIR)

pd_118553 <- pData(eset_118553)
message("  Samples: ", nrow(pd_118553))

# GSE118553 is stored as normalised expression in the Series Matrix;
# raw counts are in supplementary files.
# We use the Series Matrix VST-equivalent values for cross-platform comparison.
# Check expression range to confirm this is log-scale data:
ex_check <- exprs(eset_118553)[1:5, 1:3]
message("  Expression value range (first 5 genes): ",
        round(min(ex_check), 2), " to ", round(max(ex_check), 2))

# ── Save both raw ExpressionSets together ─────────────────────────────────────
save(eset_138260, eset_118553,
     file = file.path(RAW_DIR, "raw_esets.RData"))
message("\nBoth ExpressionSets saved to ", file.path(RAW_DIR, "raw_esets.RData"))
message("Proceed to 02_preprocess_GSE138260.R")
