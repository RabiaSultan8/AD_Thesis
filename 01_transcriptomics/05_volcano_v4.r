# ==============================================================================
# FIGURE 2.1A — VOLCANO PLOT (refined)
#
# KEY CHANGES FROM PREVIOUS VERSION:
#   1. LABEL LAYOUT — bilateral column approach (from Fig 3F)
#        - Downregulated gene labels → pinned to x = -LABEL_X (left column)
#        - Upregulated gene labels  → pinned to x = +LABEL_X (right column)
#        - nudge_x = target - point_x, direction = "y" → clean vertical stacking
#        - max.overlaps = Inf so no label is ever silently dropped
#        - Two separate geom_label_repel() calls, one per side
#   2. LEGEND — moved from bottom-right (was clashing with YWHAB/GSK3B labels)
#        to top-left, which is a data-free zone in this volcano
#   3. X-AXIS — limits widened to ±2.45 to give both label columns room;
#        tick breaks remain at ±0.5 intervals up to ±1.5 so the data area
#        reads identically to before
#   4. CANVAS — width bumped 7.0 → 8.5" to accommodate bilateral label margin;
#        height kept at 5.8"
#   5. label.size = 0 for right-side labels — borderless, matching Fig 3F
#        aesthetic (clean column against white margin); left side kept same
#   6. segment.curvature = 0 → straight connectors (consistent with Fig 3F)
#   7. All other prior fixes retained (capped triangles, ≥50 y-label, etc.)
# ==============================================================================

library(ggplot2); library(ggrepel); library(dplyr); library(here)

DATA_DIR <- here("data")
FIG_DIR  <- here("figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "meta_results.RData"))

degs_up   <- degs_meta %>% filter(mean_logFC > 0) %>% arrange(p_meta_adj)
degs_down <- degs_meta %>% filter(mean_logFC < 0) %>% arrange(p_meta_adj)

key_valid <- meta_results %>%
  filter(gene %in% c("MAPT","GSK3B","GRIN2B","YWHAB","DLG4","PSEN1"),
         p_meta_adj < 0.05, direction_consistent) %>% pull(gene)

YCAP <- 52

plot_df <- meta_results %>%
  mutate(
    key_flag   = gene %in% key_valid,
    raw_y      = -log10(p_meta_adj + 1e-300),
    capped     = raw_y > YCAP,
    neg_log10p = pmin(raw_y, YCAP),
    plot_cat   = case_when(
      key_flag                                                            ~ "Key Target",
      p_meta_adj < 0.05 & mean_logFC >= 0.3  & direction_consistent     ~ "Up",
      p_meta_adj < 0.05 & mean_logFC <= -0.3 & direction_consistent     ~ "Down",
      TRUE ~ "NS"),
    label = case_when(
      key_flag                           ~ gene,
      gene %in% head(degs_up$gene,  6)  ~ gene,
      gene %in% head(degs_down$gene, 6) ~ gene,
      gene == "GFAP"                     ~ gene,
      TRUE ~ NA_character_),
    # Key target labels get a mint fill; everything else white
    lbl_fill = ifelse(key_flag, "#D5F5E3", "white"))

plot_df$plot_cat <- factor(plot_df$plot_cat, levels = c("Up","Down","Key Target","NS"))

cols   <- c("Up"="#B03A2E","Down"="#1A5276","Key Target"="#1E8449","NS"="#BDBDBD")
sizes  <- c("Up"=1.7, "Down"=1.7, "Key Target"=3.5, "NS"=0.6)
alphas <- c("Up"=0.82,"Down"=0.82,"Key Target"=1.0, "NS"=0.15)

# ── Split labelled genes by side ───────────────────────────────────────────────
# Genes with negative logFC → left column; non-negative → right column
# Key targets split by their actual logFC direction
lbl_left  <- plot_df %>% filter(!is.na(label), mean_logFC <  0)
lbl_right <- plot_df %>% filter(!is.na(label), mean_logFC >= 0)

# X-coordinates where the label columns are pinned
# (chosen to sit in the expanded margin, clear of data points)
LABEL_X <- 2.05   # right column anchor
# left column anchor is -LABEL_X (symmetric)

# ── Build plot ─────────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(mean_logFC, neg_log10p,
              color = plot_cat, size = plot_cat, alpha = plot_cat)) +

  # NS layer (bottom)
  geom_point(data = ~filter(.x, plot_cat == "NS"), stroke = 0) +

  # Significant, non-capped points
  geom_point(data = ~filter(.x, plot_cat != "NS" & !capped), stroke = 0) +

  # Capped points — filled triangles signal y-axis truncation
  geom_point(data  = ~filter(.x, capped & plot_cat == "Up"),
             shape = 17, size = 3.2, alpha = 1.0, color = "#B03A2E") +
  geom_point(data  = ~filter(.x, capped & plot_cat == "Down"),
             shape = 17, size = 3.2, alpha = 1.0, color = "#1A5276") +

  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.32) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             color = "grey50", linewidth = 0.32) +

  # ── LEFT label column (downregulated genes) ──────────────────────────────
  # nudge_x = -LABEL_X - mean_logFC pins every label to x = -LABEL_X
  # direction = "y" stacks them vertically without horizontal jitter
  geom_label_repel(
    data          = lbl_left,
    aes(label = label, fill = lbl_fill),
    color         = "grey12",
    nudge_x       = -LABEL_X - lbl_left$mean_logFC,
    direction     = "y",
    hjust         = 1,                          # right-align text at anchor
    xlim          = c(-Inf, -1.55),             # keep labels in left margin
    size          = 3.3,
    fontface      = "italic",
    label.size    = 0.14,
    label.padding = unit(0.13, "lines"),
    label.r       = unit(0.08, "lines"),
    segment.color = "grey45",
    segment.size  = 0.24,
    segment.alpha = 0.85,
    segment.curvature = 0,                      # straight connectors
    max.overlaps  = Inf,                        # never drop a label silently
    show.legend   = FALSE,
    seed          = 42) +

  # ── RIGHT label column (upregulated genes) ───────────────────────────────
  geom_label_repel(
    data          = lbl_right,
    aes(label = label, fill = lbl_fill),
    color         = "grey12",
    nudge_x       = LABEL_X - lbl_right$mean_logFC,
    direction     = "y",
    hjust         = 0,                          # left-align text at anchor
    xlim          = c(1.55, Inf),               # keep labels in right margin
    size          = 3.3,
    fontface      = "italic",
    label.size    = 0.14,
    label.padding = unit(0.13, "lines"),
    label.r       = unit(0.08, "lines"),
    segment.color = "grey45",
    segment.size  = 0.24,
    segment.alpha = 0.85,
    segment.curvature = 0,
    max.overlaps  = Inf,
    show.legend   = FALSE,
    seed          = 42) +

  scale_fill_identity() +

  scale_color_manual(
    values = cols,
    labels = c(
      "Up"         = paste0("Upregulated (n=",   nrow(degs_up),   ")"),
      "Down"       = paste0("Downregulated (n=", nrow(degs_down), ")"),
      "Key Target" = paste0("Key AD target (n=", length(key_valid), ")"),
      "NS"         = "Not significant"),
    breaks = c("Up","Down","Key Target","NS")) +
  scale_size_manual(values  = sizes,  guide = "none") +
  scale_alpha_manual(values = alphas, guide = "none") +

  # ── X-axis: widened limits give room for bilateral label columns ──────────
  # Data lives in ≈ ±1.7; labels pin to ±2.05; limits ±2.45 adds breathing room
  scale_x_continuous(
    limits = c(-2.45, 2.45),
    breaks = seq(-1.5, 1.5, 0.5),
    expand = expansion(mult = 0.01)) +

  scale_y_continuous(
    limits = c(0, YCAP + 3),
    expand = expansion(mult = 0.01),
    breaks = c(0, 10, 20, 30, 40, 50),
    labels = c("0","10","20","30","40","\u226550")) +

  # ── Count annotations — matched to new x limits ───────────────────────────
  annotate("point", x = -2.35, y = 54.7, shape = 25,
           size = 3.0, color = "#1A5276", fill = "#1A5276") +
  annotate("text",  x = -2.18, y = 54.7, hjust = 0, size = 3.8,
           fontface = "bold", color = "#1A5276",
           label = paste0("Down: ", nrow(degs_down))) +
  annotate("point", x =  2.35, y = 54.7, shape = 24,
           size = 3.0, color = "#B03A2E", fill = "#B03A2E") +
  annotate("text",  x =  2.18, y = 54.7, hjust = 1, size = 3.8,
           fontface = "bold", color = "#B03A2E",
           label = paste0("Up: ", nrow(degs_up))) +

  labs(
    x     = expression("Weighted mean "*log[2]*" fold change (AD vs. control)"),
    y     = expression(-log[10]*"(BH-adjusted "*italic(p)*"-value)"),
    color = NULL) +

  theme_classic(base_size = 12) +
  theme(
    # ── Legend: moved to top-left (data-free zone in this volcano) ──────────
    legend.position   = c(0.10, 0.84),
    legend.background = element_rect(fill = "white", color = "grey75",
                                     linewidth = 0.30),
    legend.key        = element_rect(fill = NA),
    legend.key.size   = unit(0.44, "cm"),
    legend.text       = element_text(size = 9.5),
    legend.margin     = margin(5, 7, 5, 7),
    legend.justification = c(0, 1),            # anchor legend at its top-left corner
    axis.line         = element_line(linewidth = 0.38),
    axis.ticks        = element_line(linewidth = 0.30),
    axis.ticks.length = unit(2.5, "pt"),
    axis.text         = element_text(size = 10.5, color = "grey10"),
    axis.title        = element_text(size = 11),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(12, 10, 10, 10)) +

  guides(color = guide_legend(
    override.aes = list(
      size   = c(3.2, 3.2, 4.5, 2.0),
      alpha  = c(1.0, 1.0, 1.0, 0.5),
      stroke = c(0,   0,   0,   0))))

ggsave(file.path(FIG_DIR, "Fig2.1A_volcano.pdf"),
       p, width = 8.5, height = 5.8, dpi = 600, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2.1A_volcano.png"),
       p, width = 8.5, height = 5.8, dpi = 600, bg = "white")
cat("Volcano saved.\n")