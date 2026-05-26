# ==============================================================================
# 04_meta_analysis.R
# AD Transcriptomic Meta-Analysis — Weighted Fixed-Effects Meta-Analysis
#
# Combines DE results from GSE138260 (n=36) and GSE118553 (n=207) using:
#   - Stouffer's z-score method (sample-size weighted p-value combination)
#   - Sample-size weighted mean log2 fold change
#   - Directional consistency filter (concordant sign in both datasets)
#
# Outputs:
#   meta_results — full gene-level table (all genes)
#   degs_meta    — high-confidence DEGs (padj<0.05, |FC|>=0.3, consistent)
#   DEA_results.RData — per-dataset DE tables (for lollipop CIs)
#   meta_results.RData — meta-analysis results (used by all plot scripts)
#   normalized_matrices.RData — eset1, ex1_norm, ex2_norm, pd2_sub
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)    # for formal heterogeneity check (I², Cochran Q)
  library(here)
})

DATA_DIR <- here("data")

# ── 1. Load DE results ─────────────────────────────────────────────────────────
load(file.path(DATA_DIR, "DEA_GSE138260.RData"))   # ex1_norm, res1, pd_ord
load(file.path(DATA_DIR, "DEA_GSE118553.RData"))   # ex2_norm, res2, pd2_sub

# Sample sizes
n1 <- ncol(ex1_norm)   # 36
n2 <- ncol(ex2_norm)   # 207
message("Sample sizes: n1 = ", n1, ", n2 = ", n2)

# ── 2. Merge on gene symbol ────────────────────────────────────────────────────
common_genes <- intersect(res1$gene, res2$gene)
message("Genes in both datasets: ", length(common_genes))

df <- inner_join(
  res1 %>% select(gene, logFC_GSE138260, pval_GSE138260, padj_GSE138260,
                  t_GSE138260),
  res2 %>% select(gene, logFC_GSE118553, pval_GSE118553, padj_GSE118553,
                  t_GSE118553),
  by = "gene")

# ── 3. Stouffer weighted z-score meta-p ───────────────────────────────────────
# z_i = qnorm(1 - p_i/2) * sign(FC_i)   (signed z-score)
# z_meta = (sqrt(n1)*z1 + sqrt(n2)*z2) / sqrt(n1 + n2)
# p_meta = 2 * pnorm(-|z_meta|)

df <- df %>%
  mutate(
    z1 = qnorm(1 - pval_GSE138260 / 2) * sign(logFC_GSE138260),
    z2 = qnorm(1 - pval_GSE118553  / 2) * sign(logFC_GSE118553),
    # Clamp extreme z-scores (numerical precision ceiling at ~37 for pnorm)
    z1 = pmax(pmin(z1, 37), -37),
    z2 = pmax(pmin(z2, 37), -37),
    z_meta = (sqrt(n1) * z1 + sqrt(n2) * z2) / sqrt(n1 + n2),
    p_meta = 2 * pnorm(-abs(z_meta)),
    p_meta = pmax(p_meta, 1e-300),    # floor to avoid log(0)
    p_meta_adj = p.adjust(p_meta, method = "BH")
  )

# ── 4. Weighted mean log2 fold change ─────────────────────────────────────────
df <- df %>%
  mutate(
    mean_logFC = (n1 * logFC_GSE138260 + n2 * logFC_GSE118553) / (n1 + n2),
    direction_consistent = sign(logFC_GSE138260) == sign(logFC_GSE118553)
  )

# ── 5. Classify DEGs ───────────────────────────────────────────────────────────
df <- df %>%
  mutate(
    sig_consistent = p_meta_adj < 0.05 & direction_consistent,
    is_DEG_up   = sig_consistent & mean_logFC >= 0.3,
    is_DEG_down = sig_consistent & mean_logFC <= -0.3
  )

meta_results <- df
degs_meta    <- df %>% filter(is_DEG_up | is_DEG_down) %>%
  arrange(p_meta_adj)

message("\n── Meta-analysis summary ────────────────────────────────────────────────")
message("Total genes tested  : ", nrow(meta_results))
message("BH-significant      : ", sum(meta_results$p_meta_adj < 0.05))
message("Directionally consistent & padj<0.05: ",
        sum(meta_results$sig_consistent))
message("DEGs (|FC|≥0.3)     : ", nrow(degs_meta))
message("  Upregulated       : ", sum(degs_meta$is_DEG_up))
message("  Downregulated     : ", sum(degs_meta$is_DEG_down))

# Top 10 up and down
message("\nTop 10 upregulated:")
print(degs_meta %>% filter(is_DEG_up) %>%
  select(gene, mean_logFC, p_meta_adj) %>% head(10))

message("\nTop 10 downregulated:")
print(degs_meta %>% filter(is_DEG_down) %>%
  select(gene, mean_logFC, p_meta_adj) %>% head(10))

# ── 6. Key target summary ─────────────────────────────────────────────────────
key_genes  <- c("MAPT","GSK3B","GRIN2B","YWHAB","APP","PSEN1","GFAP",
                "DLG4","SYNGR1","GAS7","RPH3A","APLP1","KLC1",
                "CDKN2D","ELMO1")
message("\n── Key target results ───────────────────────────────────────────────────")
print(meta_results %>%
  filter(gene %in% key_genes) %>%
  select(gene, logFC_GSE138260, logFC_GSE118553,
         mean_logFC, p_meta_adj, direction_consistent) %>%
  arrange(p_meta_adj))

# ── 7. Heterogeneity check (supplementary) ────────────────────────────────────
# Compute I² for the DEGs using a simple two-study RE model
# This is informational — reported in methods as a sensitivity check
message("\n── Heterogeneity (I²) for DEGs ─────────────────────────────────────────")
degs_sub <- meta_results %>% filter(is_DEG_up | is_DEG_down)
i2_values <- sapply(seq_len(nrow(degs_sub)), function(i) {
  yi <- c(degs_sub$logFC_GSE138260[i], degs_sub$logFC_GSE118553[i])
  # Approximate SE from t-statistic: SE ≈ logFC / t
  se1 <- abs(degs_sub$logFC_GSE138260[i] /
               (degs_sub$t_GSE138260[i] + 1e-8))
  se2 <- abs(degs_sub$logFC_GSE118553[i] /
               (degs_sub$t_GSE118553[i] + 1e-8))
  vi  <- c(se1^2, se2^2)
  tryCatch({
    m <- rma(yi = yi, vi = vi, method = "DL")
    as.numeric(m$I2)
  }, error = function(e) NA_real_)
})
message("Median I² across DEGs: ", round(median(i2_values, na.rm=TRUE), 1), "%")
message("Proportion I² > 50%  : ",
        round(mean(i2_values > 50, na.rm=TRUE) * 100, 1), "%")

# ── 8. Save all outputs ────────────────────────────────────────────────────────
# meta_results.RData — primary file read by all plot scripts
save(meta_results, degs_meta,
     res1, res2,
     file = file.path(DATA_DIR, "meta_results.RData"))

# DEA_results.RData — per-dataset tables (used by lollipop for CI estimation)
save(res1, res2,
     file = file.path(DATA_DIR, "DEA_results.RData"))

# normalized_matrices.RData — expression matrices + phenotype (used by heatmap)
eset1 <- tryCatch(expr = {
    tmp <- new.env()
    load(file.path(here("data", "raw"), "GSE138260_eset.RData"),
         envir = tmp)
    tmp$eset
  }, error = function(e) NULL)

save(eset1, ex1_norm, ex2_norm, pd2_sub,
     file = file.path(DATA_DIR, "normalized_matrices.RData"))

message("\nAll outputs saved to ", DATA_DIR)
message("Files created:")
message("  meta_results.RData        — meta-analysis results (→ plot scripts)")
message("  DEA_results.RData         — per-dataset DE tables (→ lollipop)")
message("  normalized_matrices.RData — expression matrices  (→ heatmap)")
message("\nAnalysis pipeline complete. Proceed to plot scripts (05–09).")
