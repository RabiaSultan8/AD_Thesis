# ──────────────────────────────────────────────────────────────
# FIGURE 2.3A  LOLLIPOP PLOT
#
# ISSUES FIXED:
#   1. p-value labels: fontface = "plain" — p-values are NEVER italicised
#      in scientific typography; the original had fontface="italic" which
#      is wrong convention for numeric annotations
#   2. X-axis limits: replaced hard-coded padding with floor/ceiling on
#      actual CI range → prevents any CI whisker from being clipped
#   3. geom_segment alpha raised 0.55→0.68 — slightly bolder sticks
#   4. axis.line.y = element_blank() + axis.ticks.y = element_blank()
#      — removes the y-axis line; only x-axis line remains (cleaner look)
#   5. panel.grid.major.y added as hairline — gives subtle row bands
#      aiding legibility for 14 genes
#   6. legend.spacing.x added so top legend items don't crowd
#   7. bg = "white" on PNG; canvas 7.5→8.0" wide
# ──────────────────────────────────────────────────────────────

library(ggplot2); library(dplyr); library(ggrepel); library(here)

DATA_DIR <- here("data")
FIG_DIR  <- here("figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "meta_results.RData"))
load(file.path(DATA_DIR, "DEA_results.RData"))

target_genes <- c("MAPT","GSK3B","GRIN2B","YWHAB","DLG4",
                  "APP","PSEN1","GFAP","SYNGR1","GAS7",
                  "RPH3A","APLP1","KLC1","CDKN2D","ELMO1")

df <- meta_results %>%
  filter(gene %in% target_genes) %>%
  arrange(mean_logFC) %>%
  mutate(
    gene     = factor(gene, levels = gene),
    sig      = p_meta_adj < 0.05 & direction_consistent,
    category = case_when(
      gene %in% c("MAPT","GSK3B")          ~ "Primary AD target",
      gene %in% c("GRIN2B","YWHAB","DLG4") ~ "Synaptic hub",
      gene == "GFAP"                        ~ "Neuroinflammation",
      TRUE                                  ~ "Supporting DEG"),
    se1 = abs(logFC_GSE138260) /
            sqrt(-log10(pval_GSE138260 + 1e-10)),
    se2 = abs(logFC_GSE118553) /
            sqrt(-log10(pval_GSE118553 + 1e-10)),
    se_combined = sqrt((36*se1^2 + 207*se2^2) / (36+207)),
    ci_lo = mean_logFC - 1.96*se_combined,
    ci_hi = mean_logFC + 1.96*se_combined,
    padj_label = case_when(
      !sig               ~ "ns",
      p_meta_adj < 1e-10 ~ formatC(p_meta_adj, format = "e", digits = 1),
      TRUE               ~ formatC(p_meta_adj, format = "e", digits = 1)))

cat_cols <- c(
  "Primary AD target" = "#1E8449",
  "Synaptic hub"      = "#1A5276",
  "Neuroinflammation" = "#B03A2E",
  "Supporting DEG"    = "#7D6608")

# FIX: auto-compute x limits from actual CI range
x_min <- floor(min(df$ci_lo, na.rm = TRUE) * 10) / 10 - 0.15
x_max <- ceiling(max(df$ci_hi, na.rm = TRUE) * 10) / 10 + 0.15

p4 <- ggplot(df, aes(y = gene, x = mean_logFC, color = category)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.45) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted",
             color = "grey72", linewidth = 0.30) +
  geom_segment(aes(x = 0, xend = mean_logFC, y = gene, yend = gene),
               linewidth = 0.80, alpha = 0.68) +   # FIX: was 0.55
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y",
                width = 0.30, linewidth = 0.48, alpha = 0.80) +
  geom_point(aes(size = sig, shape = sig), stroke = 0.55) +
  scale_size_manual(values  = c("TRUE"=3.8,"FALSE"=2.4), guide = "none") +
  scale_shape_manual(values = c("TRUE"=19,"FALSE"=21),   guide = "none") +
  # p-value labels — left (negative FC)
  geom_label_repel(
    data          = filter(df, mean_logFC <= 0),
    aes(label = padj_label),
    direction     = "y",
    nudge_x       = x_min - 0.25 - filter(df, mean_logFC <= 0)$mean_logFC,
    xlim          = c(NA, x_min - 0.05),
    hjust         = 1,
    size          = 3.1,
    color         = "grey20",
    fontface      = "plain",      # FIX: was "italic" — p-values are not italicised
    label.size    = 0.10,
    label.padding = unit(0.10, "lines"),
    label.r       = unit(0.06, "lines"),
    segment.color = "grey55",
    segment.size  = 0.28,
    segment.linetype = 1,
    show.legend   = FALSE,
    seed          = 42) +
  # p-value labels — right (positive FC)
  geom_label_repel(
    data          = filter(df, mean_logFC > 0),
    aes(label = padj_label),
    direction     = "y",
    nudge_x       = x_max + 0.25 - filter(df, mean_logFC > 0)$mean_logFC,
    xlim          = c(x_max + 0.05, NA),
    hjust         = 0,
    size          = 3.1,
    color         = "grey20",
    fontface      = "plain",      # FIX: was "italic"
    label.size    = 0.10,
    label.padding = unit(0.10, "lines"),
    label.r       = unit(0.06, "lines"),
    segment.color = "grey55",
    segment.size  = 0.28,
    segment.linetype = 1,
    show.legend   = FALSE,
    seed          = 42) +
  scale_color_manual(values = cat_cols, name = NULL) +
  scale_x_continuous(
    limits = c(x_min - 0.45, x_max + 0.45),
    breaks = seq(-1.25, 1.25, 0.25)) +
  labs(
    x       = expression("Weighted mean "*log[2]*" fold change (AD vs. control)"),
    y       = NULL,
    caption = "Open circles: gene not meeting directional consistency threshold across both datasets.") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y         = element_text(size = 11, face = "italic", color = "grey10"),
    axis.text.x         = element_text(size = 10.5, color = "grey10"),
    axis.title.x        = element_text(size = 11),
    axis.line.y         = element_blank(),    # FIX: removes y-axis line
    axis.ticks.y        = element_blank(),    # FIX: removes y-axis ticks
    axis.ticks.length   = unit(2.5, "pt"),
    legend.position     = "top",
    legend.text         = element_text(size = 10),
    legend.key.size     = unit(0.48, "cm"),
    legend.spacing.x    = unit(0.30, "cm"),  # FIX: was missing
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.32),
    panel.grid.major.y  = element_line(color = "grey95", linewidth = 0.25),  # FIX
    plot.caption        = element_text(size = 8.5, color = "grey45",
                                       face = "italic", hjust = 0,
                                       margin = margin(t = 6)),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.margin         = margin(10, 28, 8, 28)) +
  guides(color = guide_legend(override.aes = list(size = 4.0, shape = 19)))

ggsave(file.path(FIG_DIR, "Fig2.3A_lollipop.pdf"),
       p4, width = 8.0, height = 7.0, dpi = 600, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2.3A_lollipop.png"),
       p4, width = 8.0, height = 7.0, dpi = 600, bg = "white")
cat("Lollipop saved.\n")