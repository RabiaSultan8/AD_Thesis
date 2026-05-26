# ──────────────────────────────────────────────────────────────
# FIGURE 2.2B  KEGG PATHWAY DOTPLOT
#
# ISSUES FIXED (same rationale as GO, plus):
#   1. Gene ratio legend: same dynamic-break fix — was showing only 8.0%
#   2. exclude_terms comparison now uses tolower() on both sides —
#      the original str_to_sentence() comparison failed for mixed-case
#      KEGG term names (e.g. "Intestinal immune network for IGA production")
#   3. Color gradient consistent with GO panel
#   4. Width widened to 7.5" to match GO panel
#   5. bg = "white" on PNG
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
load(file.path(DATA_DIR, "enrichment_go.RData"))   # eg_up, eg_dn, eg_bg

kegg_up <- enrichKEGG(eg_up, universe = eg_bg, organism = "hsa",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05)
kegg_dn <- enrichKEGG(eg_dn, universe = eg_bg, organism = "hsa",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_dn <- setReadable(kegg_dn, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

up_terms_keep <- c(
  "astrocyte","glia","neuroinflam","immune","cytokine","complement",
  "interleukin","extracellular matrix","collagen","fibrosis",
  "myeloid","macrophage","inflammatory","interferon","JAK","STAT",
  "NF-kB","NF.kB","ubiquitin","proteasom","unfolded protein",
  "autophagy","lysosom","prion","Alzheimer","Parkinson")

dn_terms_keep <- c(
  "synap","postsynap","neurotransmit","axon","dendrit","kinase",
  "microtubule","cytoskeleton","mitochond","oxidative phosphoryl",
  "electron transport","neurodegenerat","neuron death",
  "NMDA","glutamat","GABAer","long-term potentiation","memory",
  "dopamin","serotoni","Huntington","Alzheimer","Parkinson")

# FIX: tolower() on both sides — original str_to_sentence() missed some terms
exclude_terms <- tolower(c(
  "Intestinal immune network for IGA production",
  "Inflammatory bowel disease",
  "Autoimmune thyroid disease",
  "Type I diabetes mellitus",
  "Graft-versus-host disease",
  "Allograft rejection"))

filter_terms <- function(res, keep_patterns, exclude, n_max = 8) {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df <- df[!tolower(df$Description) %in% exclude, ]   # FIX
  pattern <- paste(keep_patterns, collapse = "|")
  df_filt <- df[grepl(pattern, df$Description, ignore.case = TRUE), ]
  if (nrow(df_filt) == 0) df_filt <- df
  df_filt %>% arrange(p.adjust) %>% head(n_max) %>%
    mutate(
      GR = sapply(GeneRatio, function(x) {
        p <- strsplit(x, "/")[[1]]; as.numeric(p[1]) / as.numeric(p[2])}),
      log10padj   = -log10(p.adjust),
      Description = str_wrap(str_to_sentence(str_trunc(Description, 50)),
                             width = 38))
}

df_up_kegg <- filter_terms(kegg_up, up_terms_keep, exclude_terms, 8) %>%
  mutate(Direction = "Upregulated")
df_dn_kegg <- filter_terms(kegg_dn, dn_terms_keep, exclude_terms, 8) %>%
  mutate(Direction = "Downregulated")

plot_data_kegg <- bind_rows(df_up_kegg, df_dn_kegg)
plot_data_kegg$Direction <- factor(plot_data_kegg$Direction,
  levels = c("Upregulated","Downregulated"))

plot_data_kegg <- plot_data_kegg %>%
  arrange(Direction, desc(log10padj)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# FIX: dynamic breaks
gr_range  <- range(plot_data_kegg$GR, na.rm = TRUE)
gr_breaks <- pretty(gr_range, n = 3)
gr_breaks <- gr_breaks[gr_breaks > 0 & gr_breaks <= gr_range[2] * 1.05]
if (length(gr_breaks) < 2) gr_breaks <- c(0.04, 0.06, 0.08)

p_kegg <- ggplot(plot_data_kegg,
  aes(x = Direction, y = Description, size = GR, color = log10padj)) +
  geom_point() +
  scale_color_gradientn(
    colors = c("#FDE8C8","#E8871A","#9B2335","#6B0000"),
    name   = expression(-log[10]*"(padj)"),
    guide  = guide_colorbar(barwidth = 0.9, barheight = 5.5,
                            ticks.linewidth = 0.5,
                            title.position  = "top")) +
  scale_size_continuous(
    range  = c(2.5, 10.5),
    name   = "Gene ratio",
    labels = scales::percent_format(accuracy = 0.1),
    breaks = gr_breaks,
    guide  = guide_legend(
      override.aes   = list(color = "grey35"),
      title.position = "top")) +
  scale_x_discrete(expand = expansion(add = 0.85)) +
  labs(x = NULL, y = NULL, title = "KEGG Pathway") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y        = element_text(size = 9.5, color = "grey10",
                                      lineheight = 0.9),
    axis.text.x        = element_text(size = 11.5, face = "bold", color = "grey10"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey75", linewidth = 0.50),
    legend.position    = "right",
    legend.title       = element_text(size = 9.5, face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key.height  = unit(0.5, "cm"),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 10, 8, 8))

n_kegg <- nrow(plot_data_kegg)
h_kegg <- max(4.5, 0.36 * n_kegg + 1.5)

ggsave(file.path(FIG_DIR, "Fig2.2B_KEGG_dotplot.pdf"),
       p_kegg, width = 7.5, height = min(h_kegg, 12), dpi = 600, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Fig2.2B_KEGG_dotplot.png"),
       p_kegg, width = 7.5, height = min(h_kegg, 12), dpi = 600, bg = "white")

save(kegg_up, kegg_dn, file = file.path(DATA_DIR, "enrichment_kegg.RData"))
cat("KEGG dotplot saved.\n")