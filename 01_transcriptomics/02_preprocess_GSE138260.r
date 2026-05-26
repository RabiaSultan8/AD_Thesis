# ==============================================================================
# 02_preprocess_GSE138260.R
# AD Transcriptomic Meta-Analysis — Microarray Preprocessing (GSE138260)
#
# Platform: Agilent custom expression microarray (GPL27556), one-colour
# Tissue:   post-mortem temporal cortex
# Samples:  17 AD, 19 controls (n=36 total)
#
# Steps:
#   1. Load raw ExpressionSet
#   2. Extract and clean phenotype data
#   3. Filter low-expression probes (IQR > 1, detected in ≥1/3 samples)
#   4. Collapse probes to gene symbols (max mean expression per gene)
#   5. Quantile normalisation (already applied in deposited data — verify)
#   6. Differential expression with limma (AD vs. Control, covariate: age+sex)
#   7. Save: ex1_norm (genes × samples), res1 (limma results)
#
# Output: data/DEA_GSE138260.RData  (relative to project root)
# ==============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
  library(here)
})

RAW_DIR  <- here("data", "raw")
DATA_DIR <- here("data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load ────────────────────────────────────────────────────────────────────
load(file.path(RAW_DIR, "raw_esets.RData"))
eset <- eset_138260
rm(eset_138260, eset_118553)

ex  <- exprs(eset)           # probe × sample matrix
pd  <- pData(eset)
fd  <- fData(eset)           # probe annotation (includes gene symbol)

message("Raw dimensions: ", nrow(ex), " probes × ", ncol(ex), " samples")

# ── 2. Phenotype data ──────────────────────────────────────────────────────────
pd$disease <- trimws(pd$`disease state:ch1`)
pd$disease  <- ifelse(grepl("^AD$|alzheimer", pd$disease, ignore.case = TRUE),
                      "AD", "Control")
pd$age     <- as.numeric(pd$`age:ch1`)
pd$sex     <- trimws(pd$`gender:ch1`)

message("Sample composition:")
print(table(pd$disease))

# Confirm no missing covariates
stopifnot(!any(is.na(pd$disease)),
          !any(is.na(pd$age)),
          !any(is.na(pd$sex)))

# ── 3. Probe filtering ─────────────────────────────────────────────────────────
# (i) Remove probes detected above background in < 1/3 of samples
# Since deposited data is pre-filtered, we apply a quantile-based floor:
# keep probes with mean expression above the 20th percentile of all probes
expr_mean    <- rowMeans(ex)
thresh_mean  <- quantile(expr_mean, 0.20)
keep_mean    <- expr_mean > thresh_mean

# (ii) IQR filter: remove probes with near-zero variance
expr_iqr     <- apply(ex, 1, IQR)
keep_iqr     <- expr_iqr > quantile(expr_iqr, 0.25)

keep_probes  <- keep_mean & keep_iqr
ex_filt      <- ex[keep_probes, ]
fd_filt      <- fd[keep_probes, ]
message("After filtering: ", nrow(ex_filt), " probes retained")

# ── 4. Gene symbol annotation and probe collapse ───────────────────────────────
# Agilent GPL27556: gene symbol column is "Gene symbol" or "GENE_SYMBOL"
symbol_col <- intersect(c("Gene symbol","GENE_SYMBOL","gene_assignment","Symbol"),
                        colnames(fd_filt))[1]
if (is.na(symbol_col)) stop("Cannot locate gene symbol column in fData.")

fd_filt$symbol <- trimws(as.character(fd_filt[[symbol_col]]))

# Remove probes with missing or non-standard gene symbols
valid_probes <- fd_filt$symbol != "" &
                !is.na(fd_filt$symbol) &
                !grepl("^---$|^NA$", fd_filt$symbol)
ex_filt      <- ex_filt[valid_probes, ]
fd_filt      <- fd_filt[valid_probes, ]

# Collapse multiple probes per gene: keep probe with highest mean expression
rownames(ex_filt) <- fd_filt$symbol
ex_by_gene    <- split(data.frame(ex_filt), rownames(ex_filt))
ex_collapsed  <- do.call(rbind, lapply(ex_by_gene, function(m) {
  if (nrow(m) == 1) return(m)
  m[which.max(rowMeans(m)), , drop = FALSE]
}))
rownames(ex_collapsed) <- names(ex_by_gene)
ex1_norm <- as.matrix(ex_collapsed)
message("After gene collapsing: ", nrow(ex1_norm), " unique genes")

# ── 5. Verify quantile normalisation ───────────────────────────────────────────
# Check column means are approximately equal (within 1% CV)
col_means <- colMeans(ex1_norm)
cv_colmeans <- sd(col_means) / mean(col_means)
message("Column mean CV (should be <0.01 if quantile-normalised): ",
        round(cv_colmeans, 4))
if (cv_colmeans > 0.05) {
  message("  Warning: column means show >5% CV — applying quantile normalisation.")
  ex1_norm <- normalizeBetweenArrays(ex1_norm, method = "quantile")
}

# ── 6. Differential expression with limma ─────────────────────────────────────
# Model: ~ 0 + disease + age + sex
# (intercept-free contrast coding; age and sex as additive covariates)
pd_ord   <- pd[colnames(ex1_norm), , drop = FALSE]

design1  <- model.matrix(~ 0 + disease + age + sex,
                          data = data.frame(
                            disease = factor(pd_ord$disease,
                                             levels = c("Control","AD")),
                            age     = pd_ord$age,
                            sex     = factor(pd_ord$sex)))
colnames(design1) <- make.names(colnames(design1))
message("Design matrix columns: ", paste(colnames(design1), collapse = ", "))

# Empirical Bayes moderated t-statistics
fit1      <- lmFit(ex1_norm, design1)
contrast1 <- makeContrasts(AD_vs_Control = diseaseAD - diseaseControl,
                            levels = design1)
fit1c     <- contrasts.fit(fit1, contrast1)
fit1e     <- eBayes(fit1c, trend = TRUE)   # trend=TRUE for microarray

res1 <- topTable(fit1e, coef = "AD_vs_Control",
                 number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(logFC_GSE138260 = logFC,
                pval_GSE138260  = P.Value,
                padj_GSE138260  = adj.P.Val,
                t_GSE138260     = t,
                B_GSE138260     = B) %>%
  select(gene, logFC_GSE138260, t_GSE138260,
         pval_GSE138260, padj_GSE138260, B_GSE138260)

message("DE results: ", nrow(res1), " genes")
message("Significant (padj<0.05): ",
        sum(res1$padj_GSE138260 < 0.05))
message("  Up  : ", sum(res1$padj_GSE138260 < 0.05 & res1$logFC_GSE138260 > 0))
message("  Down: ", sum(res1$padj_GSE138260 < 0.05 & res1$logFC_GSE138260 < 0))

# ── 7. Save ────────────────────────────────────────────────────────────────────
save(ex1_norm, res1, pd_ord,
     file = file.path(DATA_DIR, "DEA_GSE138260.RData"))
message("\nSaved: DEA_GSE138260.RData")
message("Proceed to 03_preprocess_GSE118553.R")
