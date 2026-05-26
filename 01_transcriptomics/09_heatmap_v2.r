# ──────────────────────────────────────────────────────────────
# FIGURE 2.3B  HEATMAP
#
# ISSUES FIXED:
#   1. COLOR SCALE MISMATCH (critical): colorRamp2 uses breakpoints ±2.5
#      but the heatmap_legend_param was NOT setting `at` values — so the
#      automatic legend showed ±4 (the actual data range), implying the
#      colour function spans ±4 when it does not. Fixed by adding
#      `at = c(-2.5,-1,0,1,2.5)` and `labels` with ≤/≥ clipping notation.
#   2. 5-stop colorRamp2 (blue→light-blue→white→pink→red) gives smoother
#      perceptual gradient than the original 3-stop version
#   3. Row gap: `row_split` added to separate Up-regulated and Down-regulated
#      gene blocks with a 2 mm gap — makes the direction structure visible
#   4. Row annotation: replaced anno_simple() with logical vector +
#      named colour list so a proper legend entry can be generated
#   5. pd2_sub disease capitalisation fixed: "control"→"Control" to match
#      pd1 and avoid legend having two differently-cased "Control" entries
#   6. gap = unit(4,"mm") in draw() separates the two heatmap panels
#   7. merge_legend = TRUE so Condition and Z-score legends are combined
#   8. Canvas: 8.5→9.0" wide, 7.5→8.0" tall to accommodate row titles
# ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
  library(here)
})

DATA_DIR <- here("data")
FIG_DIR  <- here("figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "normalized_matrices.RData"))
load(file.path(DATA_DIR, "meta_results.RData"))

pd1 <- eset1@phenoData@data
pd1$disease <- trimws(pd1$`disease state:ch1`)
pd1$disease <- ifelse(
  grepl("^AD$|Alzheimer|alzheimer", pd1$disease, ignore.case = TRUE),
  "AD", "Control")

top_up   <- degs_meta %>% filter(mean_logFC > 0) %>%
  arrange(p_meta_adj) %>% head(20) %>% pull(gene)
top_down <- degs_meta %>% filter(mean_logFC < 0) %>%
  arrange(p_meta_adj) %>% head(20) %>% pull(gene)
key_genes  <- c("MAPT","GSK3B","GRIN2B","YWHAB","GFAP")
show_genes <- unique(c(top_up, top_down, key_genes))

g1       <- intersect(show_genes, rownames(ex1_norm))
g2       <- intersect(show_genes, rownames(ex2_norm))
common_g <- intersect(g1, g2)

gene_order <- unique(c(
  intersect(top_up,   common_g),
  intersect(top_down, common_g),
  setdiff(intersect(key_genes, common_g),
          c(intersect(top_up, common_g), intersect(top_down, common_g)))))

z_scale <- function(mat) t(scale(t(mat)))
ex1_z <- z_scale(ex1_norm[gene_order, ])
ex2_z <- z_scale(ex2_norm[gene_order, ])

ord1  <- order(factor(pd1$disease, levels = c("Control", "AD")))
ex1_z <- ex1_z[, ord1]
dis1  <- pd1$disease[ord1]

# FIX: capitalise "Control" in pd2_sub to match pd1
pd2_sub$disease <- ifelse(
  grepl("^AD$|Alzheimer|alzheimer", pd2_sub$disease, ignore.case = TRUE),
  "AD", "Control")
ord2  <- order(factor(pd2_sub$disease, levels = c("Control", "AD")))
ex2_z <- ex2_z[, ord2]
dis2  <- pd2_sub$disease[ord2]

# FIX: 5-stop colour scale — smoother gradient; same ±2.5 range
col_fun <- colorRamp2(
  c(-2.5, -1.0, 0, 1.0, 2.5),
  c("#2471A3", "#AED6F1", "white", "#F1948A", "#B03A2E"))

# FIX: row_split creates a visible gap between Up and Down blocks
n_up_genes   <- length(intersect(top_up, gene_order))
row_split     <- factor(
  c(rep("Up", n_up_genes),
    rep("Down", length(gene_order) - n_up_genes)),
  levels = c("Up","Down"))

# FIX: row annotation — logical vector + named list for legend
row_ha <- rowAnnotation(
  Direction = gene_order %in% top_up,
  col = list(Direction = c("TRUE"  = "#B03A2E",
                           "FALSE" = "#1A5276")),
  annotation_width     = unit(2.5, "mm"),
  show_annotation_name = FALSE,
  show_legend          = FALSE)

# ── Shared condition colour map ────────────────────────────────
cond_col <- c("Control" = "#AED6F1", "AD" = "#E74C3C")

top1 <- HeatmapAnnotation(
  Condition = dis1,
  col = list(Condition = cond_col),
  annotation_height  = unit(3.5, "mm"),
  annotation_name_gp = gpar(fontsize = 0),   # shown in legend only
  show_legend        = TRUE,
  annotation_legend_param = list(
    title     = "Condition",
    title_gp  = gpar(fontsize = 9.5, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)))

top2 <- HeatmapAnnotation(
  Condition = dis2,
  col = list(Condition = cond_col),
  annotation_height    = unit(3.5, "mm"),
  show_legend          = FALSE,
  show_annotation_name = FALSE)

# ── Heatmap 1 (GSE138260) ─────────────────────────────────────
ht1 <- Heatmap(
  ex1_z, name = "Z-score",
  col               = col_fun,
  top_annotation    = top1,
  left_annotation   = row_ha,
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_names_side    = "left",
  row_names_gp      = gpar(fontsize = 9.5, fontface = "italic"),
  row_split         = row_split,    # FIX: block gap
  row_title_gp      = gpar(fontsize = 9, fontface = "bold"),
  row_gap           = unit(2, "mm"),
  cluster_columns   = FALSE,
  cluster_rows      = FALSE,
  column_title      = "GSE138260\n(post-mortem cortex)",
  column_title_gp   = gpar(fontsize = 10.5, fontface = "bold"),
  width             = unit(5, "cm"),
  border            = TRUE,
  # FIX: at= set to match colorRamp2 breakpoints; labels show clipping
  heatmap_legend_param = list(
    title         = "Z-score",
    title_gp      = gpar(fontsize = 9.5, fontface = "bold"),
    labels_gp     = gpar(fontsize = 9),
    at            = c(-2.5, -1, 0, 1, 2.5),
    labels        = c("\u2264-2.5", "-1", "0", "1", "\u22652.5"),
    grid_width    = unit(3.5, "mm"),
    legend_height = unit(4.0, "cm")))

# ── Heatmap 2 (GSE118553) ─────────────────────────────────────
ht2 <- Heatmap(
  ex2_z, name = "Z2",
  col               = col_fun,
  top_annotation    = top2,
  show_column_names = FALSE,
  show_row_names    = FALSE,
  row_split         = row_split,    # FIX: must match ht1 split
  row_title_gp      = gpar(fontsize = 0),   # suppress duplicate title
  row_gap           = unit(2, "mm"),
  cluster_columns   = FALSE,
  cluster_rows      = FALSE,
  column_title      = "GSE118553\n(iPSC neurons)",
  column_title_gp   = gpar(fontsize = 10.5, fontface = "bold"),
  width             = unit(9, "cm"),
  border            = TRUE,
  show_heatmap_legend = FALSE)

ht_list <- ht1 + ht2

draw_args <- list(
  ht_list,
  column_title    = "Top 40 DEGs \u2014 cross-dataset expression (z-scored)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  merge_legend    = TRUE,    # FIX: was FALSE — merges Condition + Z-score
  gap             = unit(4, "mm"),   # FIX: visual gap between panels
  padding         = unit(c(5, 5, 5, 5), "mm"))

pdf(file.path(FIG_DIR, "Fig2.3B_heatmap.pdf"),
    width = 9.0, height = 8.0)
do.call(draw, draw_args)
dev.off()

png(file.path(FIG_DIR, "Fig2.3B_heatmap.png"),
    width = 9.0, height = 8.0, units = "in", res = 600)
do.call(draw, draw_args)
dev.off()

cat("Heatmap saved.\n")