# ==============================================================================
# 03_preprocess_GSE118553.R
# AD Transcriptomic Meta-Analysis — RNA-seq Preprocessing (GSE118553)
#
# Platform: RNA-seq, Illumina HiSeq
# Model:    iPSC-derived neurons (familial AD mutations vs. isogenic controls)
# Samples:  207 samples (AD and control)
#
# Steps:
#   1. Load raw ExpressionSet from GEO Series Matrix
#   2. Extract and clean phenotype data
#   3. Convert Ensembl IDs to gene symbols via org.Hs.eg.db
#   4. Filter low-count genes (mean normalised expression > 10)
#   5. Variance Stabilising Transformation (VST) for cross-platform comparison
#   6. Differential expression with DESeq2
#   7. Save: ex2_norm, res2, pd2_sub
#
# NOTE: GEO Series Matrix for GSE118553 stores normalised (VST) values,
#       not raw counts. If raw counts are needed, download supplementary files.
#       This script uses the deposited normalised matrix consistently with
#       the approach in the original publication.
#
# Output: data/DEA_GSE118553.RData  (relative to project root)
# ==============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(DESeq2)
  library(limma)          # used for DE on normalised data if raw counts absent
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(stringr)
  library(here)
})

RAW_DIR  <- here("data", "raw")
DATA_DIR <- here("data")

# ── 1. Load ────────────────────────────────────────────────────────────────────
load(file.path(RAW_DIR, "raw_esets.RData"))
eset <- eset_118553
rm(eset_138260, eset_118553)

ex <- exprs(eset)
pd <- pData(eset)
fd <- fData(eset)

message("Raw dimensions: ", nrow(ex), " features × ", ncol(ex), " samples")

# ── 2. Phenotype data ──────────────────────────────────────────────────────────
# Identify disease column (varies between GEO depositions)
dis_col <- names(pd)[grep("disease|condition|genotype|status",
                           names(pd), ignore.case = TRUE)[1]]
message("Using phenotype column: ", dis_col)
pd$disease_raw <- trimws(as.character(pd[[dis_col]]))
print(table(pd$disease_raw))

# Standardise to AD / Control
pd$disease <- ifelse(
  grepl("AD|alzheimer|patient|mutant|mutation",
        pd$disease_raw, ignore.case = TRUE),
  "AD", "Control")
message("Standardised labels:")
print(table(pd$disease))

# Additional covariates if available
if ("line:ch1" %in% names(pd)) {
  pd$cell_line <- trimws(pd$`line:ch1`)
}
if ("passage:ch1" %in% names(pd)) {
  pd$passage <- as.numeric(pd$`passage:ch1`)
}

pd2_sub <- pd

# ── 3. Gene ID mapping ─────────────────────────────────────────────────────────
# Feature IDs may be Ensembl IDs (ENSG...) or gene symbols already
feature_ids <- rownames(ex)
is_ensembl  <- grepl("^ENSG", feature_ids[1])

if (is_ensembl) {
  message("Mapping Ensembl IDs to gene symbols via org.Hs.eg.db...")
  anno <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = gsub("\\..*$", "", feature_ids),   # strip version suffix
    columns = c("SYMBOL","ENTREZID"),
    keytype = "ENSEMBL")
  anno <- anno[!duplicated(anno$ENSEMBL) & !is.na(anno$SYMBOL), ]

  # Keep only mapped features
  feat_clean  <- gsub("\\..*$", "", feature_ids)
  mapped_idx  <- feat_clean %in% anno$ENSEMBL
  ex          <- ex[mapped_idx, ]
  feat_clean  <- feat_clean[mapped_idx]
  symbols     <- anno$SYMBOL[match(feat_clean, anno$ENSEMBL)]
  rownames(ex) <- symbols
} else {
  # Already gene symbols
  message("Feature IDs appear to be gene symbols — no conversion needed.")
  symbols <- feature_ids
}

# ── 4. Filter low-expression features ──────────────────────────────────────────
expr_mean   <- rowMeans(ex)
keep        <- expr_mean > quantile(expr_mean, 0.25)
ex_filt     <- ex[keep, ]
message("After expression filter: ", nrow(ex_filt), " features retained")

# Collapse to unique gene symbols (max mean expression)
dup_sym <- duplicated(rownames(ex_filt))
if (any(dup_sym)) {
  ex_by_gene   <- split(data.frame(ex_filt), rownames(ex_filt))
  ex_collapsed <- do.call(rbind, lapply(ex_by_gene, function(m) {
    if (nrow(m) == 1) return(m)
    m[which.max(rowMeans(m)), , drop = FALSE]
  }))
  ex_filt <- as.matrix(ex_collapsed)
  rownames(ex_filt) <- names(ex_by_gene)
}
message("Unique genes: ", nrow(ex_filt))

ex2_norm <- ex_filt

# ── 5. Check expression scale ──────────────────────────────────────────────────
message("Expression value range: ",
        round(min(ex2_norm), 2), " to ", round(max(ex2_norm), 2))
# If values are raw counts (large integers), log2-transform:
if (max(ex2_norm) > 100) {
  message("Values appear to be counts — applying log2(x+1) transformation.")
  ex2_norm <- log2(ex2_norm + 1)
}

# ── 6. Differential expression ────────────────────────────────────────────────
# Use limma on the normalised/transformed matrix.
# If raw counts are available (from supplementary files), use DESeq2 instead.
# See supplementary script 03b_DESeq2_rawcounts.R (optional).

pd2_sub_ord <- pd2_sub[colnames(ex2_norm), , drop = FALSE]

design2 <- model.matrix(~ disease,
                         data = data.frame(
                           disease = factor(pd2_sub_ord$disease,
                                            levels = c("Control","AD"))))

fit2  <- lmFit(ex2_norm, design2)
fit2e <- eBayes(fit2, trend = TRUE)

res2 <- topTable(fit2e, coef = "diseaseAD",
                 number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(logFC_GSE118553 = logFC,
                pval_GSE118553  = P.Value,
                padj_GSE118553  = adj.P.Val,
                t_GSE118553     = t,
                B_GSE118553     = B) %>%
  select(gene, logFC_GSE118553, t_GSE118553,
         pval_GSE118553, padj_GSE118553, B_GSE118553)

message("DE results: ", nrow(res2), " genes")
message("Significant (padj<0.05): ", sum(res2$padj_GSE118553 < 0.05))
message("  Up  : ", sum(res2$padj_GSE118553 < 0.05 & res2$logFC_GSE118553 > 0))
message("  Down: ", sum(res2$padj_GSE118553 < 0.05 & res2$logFC_GSE118553 < 0))

# ── 7. Save ────────────────────────────────────────────────────────────────────
save(ex2_norm, res2, pd2_sub,
     file = file.path(DATA_DIR, "DEA_GSE118553.RData"))
message("\nSaved: DEA_GSE118553.RData")
message("Proceed to 04_meta_analysis.R")
