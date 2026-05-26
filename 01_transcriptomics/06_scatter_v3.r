# ──────────────────────────────────────────────────────────────
# FIGURE 2.1B  CONSISTENCY SCATTER
#
# ISSUES FIXED:
#   1. "Key AD target" legend entry now shows (n=X) count for consistency
#      with Up/Down entries — was missing
#   2. Correlation annotation: plain annotate("rect") + annotate("text")
#      replaced with annotate("label") using monospaced font for alignment
#      of the two r values — avoids manual box coordinate juggling
#   3. legend.key = element_rect(fill=NA) — same fix as volcano
#   4. bg = "white" on PNG save
#   5. Canvas: 5.5×5.5 → 6.0×6.0" (coord_fixed keeps square aspect)
# ──────────────────────────────────────────────────────────────

library(ggplot2); library(ggrepel); library(dplyr); library(here)

DATA_DIR <- here("data")
FIG_DIR  <- here("figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "meta_results.RData"))

key_valid <- meta_results %>%
  filter(gene %in% c("MAPT","GSK3B","GRIN2B","YWHAB","DLG4","PSEN1"),
         p_meta_adj < 0.05, direction_consistent) %>% pull(gene)

top5_up <- meta_results %>%
  filter(p_meta_adj < 0.05, mean_logFC >= 0.3, direction_consistent) %>%
  arrange(p_meta_adj) %>% head(5) %>% pull(gene)

top5_dn <- meta_results %>%
  filter(p_meta_adj < 0.05, mean_logFC <= -0.3, direction_consistent) %>%
  arrange(p_meta_adj) %>% head(5) %>% pull(gene)

scatter_df <- meta_results %>%
  mutate(
    cat = case_when(
      gene %in% key_valid                                            ~ "Key Target",
      p_meta_adj < 0.05 & mean_logFC >= 0.3  & direction_consistent ~ "Up",
      p_meta_adj < 0.05 & mean_logFC <= -0.3 & direction_consistent ~ "Down",
      TRUE ~ "NS"),
    label = case_when(
      gene %in% key_valid ~ gene,
      gene %in% top5_up   ~ gene,
      gene %in% top5_dn   ~ gene,
      TRUE ~ NA_character_))

scatter_df$cat <- factor(scatter_df$cat, levels = c("Up","Down","Key Target","NS"))

deg_only <- filter(scatter_df, cat != "NS")
r_all    <- round(cor(meta_results$logFC_GSE138260,
                      meta_results$logFC_GSE118553, method = "pearson"), 3)
r_degs   <- round(cor(deg_only$logFC_GSE138260,
                      deg_only$logFC_GSE118553, method = "pearson"), 3)

cols   <- c("Up"="#B03A2E","Down"="#1A5276","Key Target"="#1E8449","NS"="#BDBDBD")
ax_lim <- 1.65

# FIX: "Key AD target" now shows count
n_key <- sum(scatter_df$cat == "Key Target")
legend_breaks <- c("Up","Down","Key Target")
legend_labels <- c(
  paste0("Upregulated (n=",   sum(scatter_df$cat=="Up"),   ")"),
  paste0("Downregulated (n=", sum(scatter_df$cat=="Down"), ")"),
  paste0("Key AD target (n=", n_key, ")"))

# Monospaced correlation label for column alignment
r_label <- sprintf("r (all genes)  = %.3f\nr (DEGs only) = %.3f", r_all, r_degs)

p <- ggplot(scatter_df, aes(logFC_GSE138260, logFC_GSE118553)) +
  geom_point(data = ~filter(.x, cat == "NS"),
             color = "#C0C0C0", size = 0.55, alpha = 0.18, stroke = 0) +
  geom_point(data = ~filter(.x, cat != "NS"),
             aes(color = cat), size = 1.9, alpha = 0.85, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey35", linewidth = 0.32) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.22) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.22) +
  geom_smooth(data = ~filter(.x, cat != "NS"),
              method = "lm", se = TRUE,
              color = "grey25", linewidth = 0.48,
              fill = "grey75", alpha = 0.30) +
  geom_label_repel(
    data          = ~filter(.x, !is.na(label)),
    aes(label = label, color = cat),
    size          = 3.3,
    fontface      = "italic",
    fill          = "white",
    label.size    = 0.13,
    label.padding = unit(0.13, "lines"),
    label.r       = unit(0.08, "lines"),
    box.padding   = 0.52,
    point.padding = 0.22,
    segment.size  = 0.24,
    segment.alpha = 0.72,
    max.overlaps  = 40,
    show.legend   = FALSE,
    seed          = 42) +
  # FIX: clean annotate("label") with monospaced font — no manual rect needed
  annotate("label",
           x = -ax_lim + 0.05, y = ax_lim - 0.05,
           hjust = 0, vjust = 1,
           label  = r_label,
           size   = 3.3, family = "mono",
           color  = "grey12", fill = "white",
           label.size    = 0.20,
           label.padding = unit(0.28, "lines"),
           label.r       = unit(0.06, "lines")) +
  annotate("text",
           x = 0, y = -ax_lim + 0.08,
           hjust = 0.5, size = 2.9,
           color = "grey50", fontface = "italic",
           label = "Dashed line: y = x (perfect agreement)") +
  scale_color_manual(values = cols,
                     labels = legend_labels,
                     breaks = legend_breaks) +
  coord_fixed(ratio = 1,
              xlim = c(-ax_lim, ax_lim),
              ylim = c(-ax_lim, ax_lim)) +
  labs(
    x     = expression(log[2]*" FC \u2014 GSE138260 (post-mortem cortex, n=36)"),
    y     = expression(log[2]*" FC \u2014 GSE118553 (iPSC neurons, n=207)"),
    color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = c(0.81, 0.15),
    legend.background = element_rect(fill = "white", color = "grey75",
                                     linewidth = 0.30),
    legend.key        = element_rect(fill = NA),
    legend.key.size   = unit(0.44, "cm"),
    legend.text       = element_text(size = 9.5),
    legend.margin     = margin(5, 6, 5, 6),
    axis.line         = element_line(linewidth = 0.38),
    axis.ticks        = element_line(linewidth = 0.30),
    axis.ticks.length = unit(2.5, "pt"),
    axis.text         = element_text(size = 10.5, color = "grey10"),
    axis.title        = element_text(size = 11),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(12, 16, 10, 10)) +
  guides(color = guide_legend(
    override.aes = list(
      size   = c(3.2, 3.2, 4.5),
      alpha  = c(1.0, 1.0, 1.0),
      stroke = c(0,   0,   0))))

ggsave(file.path(FIG_DIR, "Fig2.1B_consistency.pdf"),
       p, width = 6.0, height = 6.0, dpi = 600, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2.1B_consistency.png"),
       p, width = 6.0, height = 6.0, dpi = 600, bg = "white")
cat("Consistency scatter saved.\n")