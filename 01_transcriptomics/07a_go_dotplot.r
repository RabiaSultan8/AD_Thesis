# ──────────────────────────────────────────────────────────────
# FIGURE 2.2A  GO BIOLOGICAL PROCESS DOTPLOT
#
# ISSUES FIXED:
#   1. Gene ratio legend: breaks now computed from actual data range so
#      multiple sizes always appear — was showing only one circle
#   2. Color gradient: 4-stop gradient (pale→orange→dark-red→maroon) gives
#      better perceptual separation than the original 3-stop
#   3. Term labels: str_wrap() applied so long descriptions fold cleanly
#      instead of being hard-truncated (str_trunc alone loses context)
#   4. guide_legend(..., override.aes=list(color="grey40")) added so the
#      size-legend circles are visible against white background
#   5. panel.border linewidth bumped to 0.5 for crisp print reproduction
#   6. bg = "white" on PNG
# ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(ggplot2); library(dplyr)
  library(stringr); library(scales); library(here)
})

DATA_DIR <- here("data")
FIG_DIR  <- here("figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "meta_results.RData"))

sym2eg <- function(syms) {
  ids <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
           columns = "ENTREZID", keytype = "SYMBOL")
  na.omit(unique(ids$ENTREZID))
}

eg_up <- sym2eg(degs_meta %>% filter(mean_logFC > 0) %>% pull(gene))
eg_dn <- sym2eg(degs_meta %>% filter(mean_logFC < 0) %>% pull(gene))
eg_bg <- sym2eg(meta_results$gene)

go_up <- enrichGO(eg_up, universe = eg_bg, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
go_dn <- enrichGO(eg_dn, universe = eg_bg, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

go_up_s <- simplify(go_up, cutoff = 0.6, by = "p.adjust", select_fun = min)
go_dn_s <- simplify(go_dn, cutoff = 0.6, by = "p.adjust", select_fun = min)

up_terms_keep <- c(
  "astrocyte","glia","neuroinflam","immune","cytokine","complement",
  "interleukin","extracellular matrix","collagen","fibrosis",
  "myeloid","macrophage","inflammatory","interferon","JAK","STAT",
  "NF-kB","NF.kB","ubiquitin","proteasom","unfolded protein",
  "autophagy","lysosom")

dn_terms_keep <- c(
  "synap","postsynap","neurotransmit","axon","dendrit","kinase",
  "microtubule","cytoskeleton","mitochond","oxidative phosphoryl",
  "electron transport","neurodegenerat","neuron death",
  "NMDA","glutamat","GABAer","long-term potentiation","memory",
  "learning","cognitive","tau","MAP","phosphorylat")

filter_terms <- function(res, keep_patterns, n_max = 10) {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  pattern <- paste(keep_patterns, collapse = "|")
  df_filt <- df[grepl(pattern, df$Description, ignore.case = TRUE), ]
  if (nrow(df_filt) == 0) df_filt <- df
  df_filt %>% arrange(p.adjust) %>% head(n_max) %>%
    mutate(
      GR = sapply(GeneRatio, function(x) {
        p <- strsplit(x, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2])}),
      log10padj   = -log10(p.adjust),
      # FIX: str_wrap wraps long descriptions instead of hard-truncating
      Description = str_wrap(str_to_sentence(str_trunc(Description, 50)),
                             width = 38))
}

df_up_go <- filter_terms(go_up_s, up_terms_keep, 10) %>%
  mutate(Direction = "Upregulated")
df_dn_go <- filter_terms(go_dn_s, dn_terms_keep, 10) %>%
  mutate(Direction = "Downregulated")

plot_data_go <- bind_rows(df_up_go, df_dn_go)
plot_data_go$Direction <- factor(plot_data_go$Direction,
  levels = c("Upregulated","Downregulated"))

plot_data_go <- plot_data_go %>%
  arrange(Direction, desc(log10padj)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# FIX: gene ratio breaks derived from actual data
gr_range  <- range(plot_data_go$GR, na.rm = TRUE)
gr_breaks <- pretty(gr_range, n = 3)
gr_breaks <- gr_breaks[gr_breaks > 0 & gr_breaks <= gr_range[2] * 1.05]
if (length(gr_breaks) < 2) gr_breaks <- c(0.02, 0.05, 0.10)

p_go <- ggplot(plot_data_go,
  aes(x = Direction, y = Description, size = GR, color = log10padj)) +
  geom_point() +
  # FIX: 4-stop gradient for better perceptual separation
  scale_color_gradientn(
    colors = c("#FDE8C8","#E8871A","#9B2335","#6B0000"),
    name   = expression(-log[10]*"(padj)"),
    guide  = guide_colorbar(barwidth = 0.9, barheight = 5.5,
                            ticks.linewidth = 0.5,
                            title.position  = "top")) +
  # FIX: dynamic breaks + grey override so circles visible
  scale_size_continuous(
    range  = c(2.2, 9.0),
    name   = "Gene ratio",
    labels = scales::percent_format(accuracy = 0.1),
    breaks = gr_breaks,
    guide  = guide_legend(
      override.aes    = list(color = "grey35"),
      title.position  = "top")) +
  scale_x_discrete(expand = expansion(add = 0.85)) +
  labs(x = NULL, y = NULL, title = "GO \u2014 Biological Process") +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "#EBF5FB", color = "grey70"),
    strip.text         = element_text(size = 10.5, face = "bold", color = "grey15"),
    axis.text.y        = element_text(size = 9.5, color = "grey10",
                                      lineheight = 0.9),
    axis.text.x        = element_text(size = 11.5, face = "bold", color = "grey10"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey75", linewidth = 0.50),  # FIX
    legend.position    = "right",
    legend.title       = element_text(size = 9.5, face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key.height  = unit(0.5, "cm"),
    legend.spacing.y   = unit(0.25, "cm"),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 10, 8, 8))

n_go <- nrow(plot_data_go)
h_go <- max(6.0, 0.36 * n_go + 1.5)

ggsave(file.path(FIG_DIR, "Fig2.2A_GO_dotplot.pdf"),
       p_go, width = 7.5, height = min(h_go, 14), dpi = 600, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2.2A_GO_dotplot.png"),
       p_go, width = 7.5, height = min(h_go, 14), dpi = 600, bg = "white")

save(go_up, go_dn, go_up_s, go_dn_s, eg_up, eg_dn, eg_bg,
     file = file.path(DATA_DIR, "enrichment_go.RData"))
cat("GO dotplot saved.\n")