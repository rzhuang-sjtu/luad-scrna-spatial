#!/usr/bin/env Rscript
# Figure 5 v2 — Neutrophil heterogeneity in LUAD
# strict-aligned with Fig 1/2/3/4 templates per user 2026-04-27 spec.
# Drop showtext, use sans. No labs(title=...). No composite.
# Output dir: ~/luad/results/fig5_panels/  (overwrite v1)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(survival)
  library(survminer)
  library(ggrastr)
  library(showtext)
  library(sysfonts)
})

DAT <- "${PROJECT_ROOT}/results/fig5_plot_data"
OUT <- "${PROJECT_ROOT}/results/fig5_panels"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── Font: Arial via showtext (fall back to sans if missing) ──
.find_arial <- function() {
  cands <- c("~/.local/share/fonts/arial.ttf",
             "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
             "/mnt/c/Windows/Fonts/arial.ttf")
  for (p in cands) {
    pp <- path.expand(p); if (file.exists(pp)) return(pp)
  }
  return(NA_character_)
}
.arial_reg <- .find_arial()
.font_ok <- FALSE
if (!is.na(.arial_reg)) {
  .arial_dir <- dirname(.arial_reg)
  .arial_bold <- file.path(.arial_dir, "arialbd.ttf")
  .arial_italic <- file.path(.arial_dir, "ariali.ttf")
  if (!file.exists(.arial_bold))   .arial_bold <- .arial_reg
  if (!file.exists(.arial_italic)) .arial_italic <- .arial_reg
  .font_ok <- tryCatch({
    sysfonts::font_add("Arial",
      regular = .arial_reg, bold = .arial_bold, italic = .arial_italic)
    showtext::showtext_auto(); showtext::showtext_opts(dpi = 300)
    TRUE
  }, error = function(e) FALSE)
}
my_font <- if (.font_ok) "Arial" else "sans"
cat("[font] Arial registered:", .font_ok, "| using:", my_font,
    "| path:", if (!is.na(.arial_reg)) .arial_reg else "(none)", "\n")

theme_pub <- function(base_size = 8) {
  theme_classic(base_family = my_font, base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "pt"),
      plot.title = element_blank(),
      legend.text = element_text(size = rel(0.9)),
      legend.title = element_text(size = rel(0.9)),
      legend.key.height = unit(8, "pt"),
      legend.key.width = unit(8, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(0.95))
    )
}

neu_colors <- c(
  "Neu_Inflammatory"   = "#E64B35",
  "Neu_Angiogenic"     = "#F39B7F",
  "Neu_Metastatic"     = "#3C5488",
  "Neu_ECM_remodeling" = "#4DBBD5",
  "Neu_OSM_priming"    = "#00A087",
  "Neu_OSM_low"        = "#8491B4",
  "Neu_IFN_response"   = "#91D1C2",
  "Neu_unclassified"   = "#D9D9D9"
)
neu_labels <- c(
  "Neu_Inflammatory"   = "Inflammatory",
  "Neu_Angiogenic"     = "Angiogenic",
  "Neu_Metastatic"     = "Metastatic",
  "Neu_ECM_remodeling" = "ECM remod.",
  "Neu_OSM_priming"    = "OSM priming",
  "Neu_OSM_low"        = "OSM low",
  "Neu_IFN_response"   = "IFN response",
  "Neu_unclassified"   = "Unclassified"
)
neu_order <- c("Neu_Inflammatory", "Neu_Angiogenic", "Neu_Metastatic",
               "Neu_IFN_response", "Neu_OSM_priming", "Neu_OSM_low",
               "Neu_ECM_remodeling", "Neu_unclassified")
neu_order_no_un <- neu_order[neu_order != "Neu_unclassified"]

mp_colors <- c("MP1" = "#E64B35", "MP2" = "#4DBBD5",
               "MP3" = "#00A087", "MP4" = "#3C5488")

tissue_labels <- c(
  "Normal_Lung"        = "Normal",
  "Adjacent_Normal"    = "Adjacent",
  "Normal_LN"          = "Normal LN",
  "Precancerous"       = "Precanc.",
  "Primary_Tumor"      = "Tumor",
  "LN_Metastasis"      = "LN met",
  "Distant_Metastasis" = "Dist. met",
  "Brain_Metastasis"   = "Brain met",
  "Pleural_Effusion"   = "PE"
)
tissue_order <- names(tissue_labels)

cat("[init] OUT =", OUT, "\n")

# ── UMAP corner-arrow helper ──
umap_arrow_axes <- function(data, x_col, y_col,
                            frac = 0.16, label_x = "UMAP 1", label_y = "UMAP 2",
                            text_size = 2.4, line_size = 0.45, arrow_mm = 1.6,
                            inset_frac = 0.02) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr); y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr);          y1 <- y0 + frac * diff(yr)
  arr <- grid::arrow(length = grid::unit(arrow_mm, "mm"), ends = "last", type = "closed")
  fam <- if (exists("my_font")) my_font else "sans"
  list(
    annotate("segment", x=x0, xend=x1, y=y0, yend=y0, arrow=arr, linewidth=line_size, color="black"),
    annotate("segment", x=x0, xend=x0, y=y0, yend=y1, arrow=arr, linewidth=line_size, color="black"),
    annotate("text", x=(x0+x1)/2, y=y0, label=label_x, vjust=2.4, size=text_size, family=fam),
    annotate("text", x=x0, y=(y0+y1)/2, label=label_y, angle=90, vjust=-1.4, size=text_size, family=fam),
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          panel.grid=element_blank())
  )
}

# 5A — Neutrophil UMAP
cat("[5A] UMAP\n")
df_a <- fread(file.path(DAT, "fig5a_umap_metadata.csv.gz"))
df_a$neu_subtype <- factor(df_a$neu_subtype, levels = neu_order)
set.seed(42)
df_a <- df_a[sample(nrow(df_a)), ]

p5a <- ggplot(as.data.frame(df_a), aes(UMAP1, UMAP2,
              color = factor(neu_subtype, levels = neu_order))) +
  ggrastr::rasterise(geom_point(size = 0.45, stroke = 0, alpha = 0.55,
                                shape = 16), dpi = 600) +
  scale_color_manual(values = neu_colors, labels = neu_labels, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL) +
  theme_pub(8) +
  theme(plot.background = element_blank(), panel.background = element_blank()) +
  umap_arrow_axes(as.data.frame(df_a), "UMAP1", "UMAP2")

ggsave(file.path(OUT, "fig5a_umap.pdf"), p5a, width = 100, height = 85, units = "mm")
ggsave(file.path(OUT, "fig5a_umap.png"), p5a, width = 100, height = 85, units = "mm", dpi = 300)

# 5B — Tissue stacked bar
cat("[5B] tissue stacked bar\n")
df_b <- fread(file.path(DAT, "fig5b_tissue_proportion.csv"))
df_b$tissue_type <- factor(df_b$tissue_type,
                           levels = intersect(tissue_order, unique(df_b$tissue_type)))

p5b <- ggplot(as.data.frame(df_b),
              aes(x = tissue_type, y = proportion_pct,
                  fill = factor(neu_subtype, levels = rev(neu_order)))) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = neu_colors, labels = neu_labels,
                    name = NULL, breaks = neu_order) +
  scale_x_discrete(labels = tissue_labels) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

ggsave(file.path(OUT, "fig5b_tissue_bar.pdf"), p5b, width = 110, height = 80, units = "mm")
ggsave(file.path(OUT, "fig5b_tissue_bar.png"), p5b, width = 110, height = 80, units = "mm", dpi = 300)

# 5C — Canonical marker dot plot
cat("[5C] canonical marker dot\n")
df_c <- fread(file.path(DAT, "fig5c_canonical_markers.csv"))
df_c$neu_subtype <- factor(df_c$neu_subtype, levels = rev(neu_order))

# normalize pct_expressing to 0-100
if (max(df_c$pct_expressing, na.rm = TRUE) <= 1.5) {
  df_c$pct_expressing <- df_c$pct_expressing * 100
}

# per-gene z-score
df_c <- df_c %>%
  group_by(gene) %>%
  mutate(mean_z = (mean_expression - mean(mean_expression)) /
         max(sd(mean_expression), 1e-6)) %>%
  ungroup()

# gene order: by canonical_for in neu_order
df_c$canonical_for <- factor(df_c$canonical_for, levels = neu_order)
gene_levels <- df_c %>%
  group_by(gene) %>% slice_head(n = 1) %>%
  arrange(canonical_for, gene) %>% pull(gene) %>% unique()
df_c$gene <- factor(df_c$gene, levels = gene_levels)

p5c <- ggplot(as.data.frame(df_c), aes(gene, neu_subtype)) +
  geom_point(aes(size = pct_expressing, fill = mean_z),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(
    name = "% expressed",
    range = c(0.3, 5.5),
    breaks = c(25, 50, 75),
    limits = c(0, 100)
  ) +
  scale_fill_gradient2(
    low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish,
    name = "Scaled expr.",
    guide = guide_colorbar(
      barwidth = unit(3, "mm"), barheight = unit(15, "mm"),
      frame.colour = "black", frame.linewidth = 0.3
    )
  ) +
  scale_y_discrete(labels = neu_labels) +
  labs(x = NULL, y = NULL) +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank()
  )

ggsave(file.path(OUT, "fig5c_canonical_dotplot.pdf"), p5c,
       width = 200, height = 85, units = "mm")
ggsave(file.path(OUT, "fig5c_canonical_dotplot.png"), p5c,
       width = 200, height = 85, units = "mm", dpi = 300)

# 5D — EMT ligand dot plot (faceted by gene_family)
cat("[5D] EMT ligand dot\n")
df_d <- fread(file.path(DAT, "fig5d_emt_ligand_dotplot.csv"))
df_d$neu_subtype <- factor(df_d$neu_subtype, levels = rev(neu_order))

if (max(df_d$pct_expressing, na.rm = TRUE) <= 1.5) {
  df_d$pct_expressing <- df_d$pct_expressing * 100
}

df_d <- df_d %>%
  group_by(gene) %>%
  mutate(mean_z = (mean_expression - mean(mean_expression)) /
         max(sd(mean_expression), 1e-6)) %>%
  ungroup()

emt_gene_order <- c(
  "TGFB1",
  "TNF","IL1A","IL1B","IL6","OSM",
  "CXCL1","CXCL2","CXCL8","CCL2","CCL3","CCL4","CCL5",
  "MMP9",
  "VEGFA","VEGFB","FN1","SPP1","SERPINE1","PLAU","PLAUR",
  "AREG","EREG","WNT5A","PDGFB"
)
df_d$gene <- factor(df_d$gene, levels = intersect(emt_gene_order, unique(df_d$gene)))
df_d$gene_family <- factor(df_d$gene_family,
                           levels = c("TGFb","TNF/IL","Chemokine","MMP",
                                      "Angio_EMT","EGF_other"))

p5d <- ggplot(as.data.frame(df_d), aes(gene, neu_subtype)) +
  geom_point(aes(size = pct_expressing, fill = mean_z),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(
    name = "% expressed",
    range = c(0.3, 5.5),
    breaks = c(25, 50, 75),
    limits = c(0, 100)
  ) +
  scale_fill_gradient2(
    low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish,
    name = "Scaled expr.",
    guide = guide_colorbar(
      barwidth = unit(3, "mm"), barheight = unit(15, "mm"),
      frame.colour = "black", frame.linewidth = 0.3
    )
  ) +
  scale_y_discrete(labels = neu_labels) +
  facet_grid(. ~ gene_family, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank(),
    panel.spacing.x = unit(4, "pt")
  )

ggsave(file.path(OUT, "fig5d_emt_dotplot.pdf"), p5d,
       width = 200, height = 85, units = "mm")
ggsave(file.path(OUT, "fig5d_emt_dotplot.png"), p5d,
       width = 200, height = 85, units = "mm", dpi = 300)

# 5E — EMT-promoting ligand specificity heatmap
#       per-gene z-score of mean expression across Neu subtypes,
#       grouped by ligand pathway (TGFb / TNF-IL / Chemokine / MMP /
#       Angio-EMT / Other). Replaces the sparse-logFC matrix.
cat("[5E] EMT ligand specificity heatmap (z-scored mean expression)\n")

mean_raw <- as.data.frame(fread(file.path(DAT, "fig5e_ligand_mean_raw.csv")))
rename_map <- as.data.frame(fread(file.path(DAT, "fig5e_rename_map.csv")))

# Long form, apply rename map, drop unclassified
mean_long <- pivot_longer(mean_raw, -scanvi_predicted,
                          names_to = "gene", values_to = "mean_expr")
mean_long <- mean_long %>%
  rename(scanvi_label = scanvi_predicted) %>%
  inner_join(rename_map, by = "scanvi_label") %>%
  filter(neu_subtype != "Neu_unclassified")

# Drop genes that are essentially silent across all subtypes (max < 0.05
# log1p mean), and genes that don't vary (sd < 0.02)
gene_stats <- mean_long %>%
  group_by(gene) %>%
  summarize(max_e = max(mean_expr, na.rm = TRUE),
            sd_e  = sd(mean_expr, na.rm = TRUE), .groups = "drop")
keep_g <- gene_stats$gene[gene_stats$max_e > 0.05 & gene_stats$sd_e > 0.02]
mean_long <- mean_long %>% filter(gene %in% keep_g)

# Per-gene z-score across Neu subtypes
mean_long <- mean_long %>%
  group_by(gene) %>%
  mutate(z_expr = (mean_expr - mean(mean_expr)) /
                  pmax(sd(mean_expr), 1e-6)) %>%
  ungroup()

# Wide z-matrix (gene × subtype)
m_e <- mean_long %>%
  select(gene, neu_subtype, z_expr) %>%
  pivot_wider(names_from = neu_subtype, values_from = z_expr) %>%
  as.data.frame()
rownames(m_e) <- m_e$gene; m_e$gene <- NULL
present_subtypes <- intersect(neu_order_no_un, colnames(m_e))
m_e <- as.matrix(m_e[, present_subtypes, drop = FALSE])
m_e[!is.finite(m_e)] <- 0
colnames(m_e) <- neu_labels[colnames(m_e)]

# Pathway family annotation per gene
pathway_map <- c(
  "TGFB1" = "TGFb", "TGFB2" = "TGFb", "TGFB3" = "TGFb",
  "TNF" = "TNF/IL", "IL6" = "TNF/IL", "IL1A" = "TNF/IL",
  "IL1B" = "TNF/IL", "OSM" = "TNF/IL",
  "CXCL1" = "Chemokine", "CXCL2" = "Chemokine", "CXCL5" = "Chemokine",
  "CXCL8" = "Chemokine", "CCL2" = "Chemokine", "CCL3" = "Chemokine",
  "CCL4" = "Chemokine", "CCL5" = "Chemokine",
  "MMP9" = "MMP", "MMP2" = "MMP", "MMP8" = "MMP", "MMP25" = "MMP",
  "VEGFA" = "Angio/EMT", "VEGFB" = "Angio/EMT", "FN1" = "Angio/EMT",
  "SPP1" = "Angio/EMT", "SERPINE1" = "Angio/EMT",
  "PLAU" = "Angio/EMT", "PLAUR" = "Angio/EMT",
  "HGF" = "Other", "AREG" = "Other", "EREG" = "Other",
  "WNT5A" = "Other", "PDGFB" = "Other", "FGF2" = "Other"
)
fam_levels <- c("TGFb", "TNF/IL", "Chemokine", "MMP", "Angio/EMT", "Other")
fam_pal <- c("TGFb"      = "#7E57C2", "TNF/IL"   = "#E64B35",
             "Chemokine" = "#F39B7F", "MMP"      = "#00A087",
             "Angio/EMT" = "#3C5488", "Other"    = "grey75")
gene_fam <- factor(unname(pathway_map[rownames(m_e)]), levels = fam_levels)

# Order rows by family, then by max z-score within family
ord <- order(gene_fam, -apply(m_e, 1, max))
m_e <- m_e[ord, , drop = FALSE]
gene_fam <- gene_fam[ord]

ht_pal <- colorRamp2(c(-2, 0, 2), c("#3C5488", "#F7F7F7", "#E64B35"))

ha_col_e <- HeatmapAnnotation(
  Subtype = colnames(m_e),
  col = list(Subtype = setNames(neu_colors[present_subtypes],
                                neu_labels[present_subtypes])),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(2.5, "mm")
)
ha_row_e <- rowAnnotation(
  Pathway = gene_fam,
  col = list(Pathway = fam_pal),
  show_annotation_name = FALSE, show_legend = TRUE,
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(
    Pathway = list(
      title_gp = gpar(fontsize = 7, fontfamily = my_font, fontface = "bold"),
      labels_gp = gpar(fontsize = 6, fontfamily = my_font),
      grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")
    )
  )
)

ht_e <- Heatmap(
  m_e, name = "z-score", col = ht_pal,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = ha_col_e, left_annotation = ha_row_e,
  row_split = gene_fam, row_title = NULL, row_gap = unit(1, "mm"),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6.5, fontfamily = my_font, fontface = "italic"),
  column_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_rot = 45,
  border = TRUE, border_gp = gpar(col = "black", lwd = 0.8),
  rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = my_font, fontface = "bold"),
    labels_gp = gpar(fontsize = 6, fontfamily = my_font),
    legend_height = unit(22, "mm"), grid_width = unit(2.5, "mm")
  ),
  width  = unit(ncol(m_e) * 7, "mm"),
  height = unit(nrow(m_e) * 5, "mm")
)
.canvas_w <- ncol(m_e) * 7 + 55      # left annotation + legend padding
.canvas_h <- nrow(m_e) * 5 + 25
pdf(file.path(OUT, "fig5e_emt_logfc_heatmap.pdf"),
    width = .canvas_w/25.4, height = .canvas_h/25.4)
ComplexHeatmap::draw(ht_e, heatmap_legend_side = "right",
                     annotation_legend_side = "right",
                     merge_legend = TRUE,
                     padding = unit(c(2, 3, 2, 3), "mm"))
dev.off()
png(file.path(OUT, "fig5e_emt_logfc_heatmap.png"),
    width = .canvas_w, height = .canvas_h, units = "mm", res = 300)
ComplexHeatmap::draw(ht_e, heatmap_legend_side = "right",
                     annotation_legend_side = "right",
                     merge_legend = TRUE,
                     padding = unit(c(2, 3, 2, 3), "mm"))
dev.off()

# 5F — Neu × MP correlation heatmap (single-cell, n=58)
cat("[5F] Neu x MP correlation heatmap\n")
df_f <- fread(file.path(DAT, "fig5f_neu_mp_correlation.csv"))
m_f <- df_f %>%
  pivot_wider(id_cols = neu_subtype, names_from = MP, values_from = spearman_rho) %>%
  as.data.frame()
rownames(m_f) <- m_f$neu_subtype; m_f$neu_subtype <- NULL
m_f <- m_f[neu_order_no_un, , drop = FALSE]
m_f <- as.matrix(m_f)

p_f <- df_f %>%
  pivot_wider(id_cols = neu_subtype, names_from = MP, values_from = pvalue) %>%
  as.data.frame()
rownames(p_f) <- p_f$neu_subtype; p_f$neu_subtype <- NULL
p_f <- as.matrix(p_f[neu_order_no_un, , drop = FALSE])

# rename row labels
rownames(m_f) <- neu_labels[rownames(m_f)]
rownames(p_f) <- rownames(m_f)

cor_pal <- colorRamp2(c(-0.6, 0, 0.6), c("#3C5488", "#F7F7F7", "#E64B35"))
ht_f <- Heatmap(
  m_f, name = "rho", col = cor_pal,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_gp = gpar(fontsize = 8, fontfamily = my_font),
  column_names_rot = 0, column_names_centered = TRUE,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    rho <- m_f[i, j]; pv <- p_f[i, j]
    star <- ifelse(pv < 0.001, "***",
            ifelse(pv < 0.01, "**",
            ifelse(pv < 0.05, "*", "")))
    txt_col <- ifelse(abs(rho) > 0.35, "white", "black")
    # number on lower line
    grid.text(sprintf("%.2f", rho), x, y - unit(0.6, "mm"),
              gp = gpar(fontsize = 6, fontfamily = my_font, col = txt_col))
    # stars stacked on upper line
    if (nchar(star) > 0)
      grid.text(star, x, y + unit(1.6, "mm"),
                gp = gpar(fontsize = 6.5, fontfamily = my_font,
                          fontface = "bold", col = txt_col))
  },
  heatmap_legend_param = list(
    at = c(-0.6, -0.3, 0, 0.3, 0.6),
    title_gp = gpar(fontsize = 7, fontfamily = my_font),
    labels_gp = gpar(fontsize = 6, fontfamily = my_font)
  ),
  width  = unit(ncol(m_f) * 7, "mm"),
  height = unit(nrow(m_f) * 7, "mm")
)
.cellsize_f <- 7
# width = grid + left rownames (~28mm Neu labels) + right legend (~25mm)
.cw_f <- ncol(m_f) * .cellsize_f + 55
.ch_f <- nrow(m_f) * .cellsize_f + 22
pdf(file.path(OUT, "fig5f_neu_mp_correlation_heatmap.pdf"),
    width = .cw_f/25.4, height = .ch_f/25.4)
ComplexHeatmap::draw(ht_f, heatmap_legend_side = "right",
                     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
png(file.path(OUT, "fig5f_neu_mp_correlation_heatmap.png"),
    width = .cw_f, height = .ch_f, units = "mm", res = 300)
ComplexHeatmap::draw(ht_f, heatmap_legend_side = "right",
                     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# 5G — LIANA simplified dot plot (2 senders × MP1/MP2 × top 15 LR)
cat("[5G] LIANA dot plot (simplified)\n")
df_g <- fread(file.path(DAT, "fig5g_liana_focus_all_senders.csv"))
df_g$lr <- paste0(df_g$ligand, " > ", df_g$receptor)
df_g$receiver_short <- sub("Mal_", "", df_g$receiver)

# Per spec: 2 senders × MP1/MP2
df_g_slim <- df_g[df_g$sender %in% c("Neu_Inflammatory", "Neu_Metastatic") &
                   df_g$receiver_short %in% c("MP1", "MP2"), ]
df_g_slim <- df_g_slim %>%
  group_by(sender, receiver_short) %>%
  slice_min(magnitude_rank, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.table()

# pathway color (NPG-aligned)
pw_colors <- c(
  "TGFb"        = "#E64B35",
  "IL1"         = "#F39B7F",
  "TNF"         = "#FFCB5C",
  "OSM"         = "#00A087",
  "CXCL8"       = "#3C5488",
  "CXCL_other"  = "#8491B4",
  "CCL"         = "#4DBBD5",
  "PLAU"        = "#7E57C2",
  "VEGFA"       = "#9467BD",
  "EGFR_ligand" = "#B09C85",
  "Other"       = "#D9D9D9"
)
df_g_slim$pathway <- factor(df_g_slim$pathway, levels = names(pw_colors))

# Compose facet labels — short to avoid clipping
sender_short <- c("Neu_Inflammatory" = "Inflam", "Neu_Metastatic" = "Metastatic")
df_g_slim$facet_label <- paste0(sender_short[df_g_slim$sender], "  >>  ",
                                  df_g_slim$receiver_short)
facet_order <- c("Inflam  >>  MP1", "Inflam  >>  MP2",
                 "Metastatic  >>  MP1", "Metastatic  >>  MP2")
df_g_slim$facet_label <- factor(df_g_slim$facet_label, levels = facet_order)

# Order LR pairs within each facet by magnitude_rank
df_g_slim <- df_g_slim %>%
  group_by(facet_label) %>%
  arrange(magnitude_rank, .by_group = TRUE) %>%
  ungroup() %>%
  as.data.table()
df_g_slim$lr <- factor(df_g_slim$lr,
                       levels = rev(unique(df_g_slim$lr[order(df_g_slim$facet_label,
                                                              df_g_slim$magnitude_rank)])))

p5g <- ggplot(as.data.frame(df_g_slim), aes(x = 1, y = lr)) +
  geom_point(aes(size = -log10(pmax(magnitude_rank, 1e-3)),
                 fill = pathway),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(range = c(1.5, 5), name = "-log10(rank)") +
  scale_fill_manual(values = pw_colors, name = "Pathway", drop = FALSE) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_pub(8) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank()
  )

ggsave(file.path(OUT, "fig5g_liana.pdf"), p5g,
       width = 160, height = 120, units = "mm")
ggsave(file.path(OUT, "fig5g_liana.png"), p5g,
       width = 160, height = 120, units = "mm", dpi = 300)

# 5H — TCGA correlation heatmap
cat("[5H] TCGA correlation heatmap\n")
df_h <- fread(file.path(DAT, "fig5h_tcga_correlation_matrix.csv"))
m_h <- df_h %>%
  pivot_wider(id_cols = neu_subtype, names_from = MP, values_from = spearman_rho) %>%
  as.data.frame()
rownames(m_h) <- m_h$neu_subtype; m_h$neu_subtype <- NULL
m_h <- as.matrix(m_h)

# row order: Neu subtypes first, EMT_Hallmark last
row_order_h <- c(intersect(neu_order_no_un, rownames(m_h)),
                 intersect("EMT_Hallmark", rownames(m_h)))
m_h <- m_h[row_order_h, , drop = FALSE]

p_h <- df_h %>%
  pivot_wider(id_cols = neu_subtype, names_from = MP, values_from = pvalue) %>%
  as.data.frame()
rownames(p_h) <- p_h$neu_subtype; p_h$neu_subtype <- NULL
p_h <- as.matrix(p_h[row_order_h, , drop = FALSE])

# rename row labels
disp_rows <- ifelse(rownames(m_h) == "EMT_Hallmark", "EMT (Hallmark)",
                    neu_labels[rownames(m_h)])
rownames(m_h) <- disp_rows
rownames(p_h) <- disp_rows

cor_pal_h <- colorRamp2(c(-0.8, 0, 0.8), c("#3C5488", "#F7F7F7", "#E64B35"))
ht_h <- Heatmap(
  m_h, name = "rho", col = cor_pal_h,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_gp = gpar(fontsize = 8, fontfamily = my_font),
  column_names_rot = 0, column_names_centered = TRUE,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    rho <- m_h[i, j]; pv <- p_h[i, j]
    star <- ifelse(is.na(pv) | pv >= 0.05, "",
            ifelse(pv < 0.001, "***",
            ifelse(pv < 0.01, "**", "*")))
    txt_col <- ifelse(abs(rho) > 0.45, "white", "black")
    grid.text(sprintf("%.2f", rho), x, y - unit(0.6, "mm"),
              gp = gpar(fontsize = 6, fontfamily = my_font, col = txt_col))
    if (nchar(star) > 0)
      grid.text(star, x, y + unit(1.6, "mm"),
                gp = gpar(fontsize = 6.5, fontfamily = my_font,
                          fontface = "bold", col = txt_col))
  },
  heatmap_legend_param = list(
    at = c(-0.6, 0, 0.6),
    title_gp = gpar(fontsize = 7, fontfamily = my_font),
    labels_gp = gpar(fontsize = 6, fontfamily = my_font)
  ),
  width  = unit(ncol(m_h) * 7, "mm"),
  height = unit(nrow(m_h) * 7, "mm")
)
.cellsize_h <- 7
# width = grid + left rownames (~32mm "EMT (Hallmark)") + right legend (~25mm)
.cw_h <- ncol(m_h) * .cellsize_h + 60
.ch_h <- nrow(m_h) * .cellsize_h + 22
pdf(file.path(OUT, "fig5h_tcga_correlation_heatmap.pdf"),
    width = .cw_h/25.4, height = .ch_h/25.4)
ComplexHeatmap::draw(ht_h, heatmap_legend_side = "right",
                     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
png(file.path(OUT, "fig5h_tcga_correlation_heatmap.png"),
    width = .cw_h, height = .ch_h, units = "mm", res = 300)
ComplexHeatmap::draw(ht_h, heatmap_legend_side = "right",
                     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# 5I — KM 4-group (Neu_Metastatic × MP2), NPG palette + short labels
cat("[5I] KM 4-group\n")
df_i <- fread(file.path(DAT, "fig5i_km_data.csv"))
df_i[, OS_years := OS_days / 365.25]

short_lab <- function(x) {
  x <- sub("Low_NeuMet/Low_MP2",   "Met-/MP2-", x)
  x <- sub("Low_NeuMet/High_MP2",  "Met-/MP2+", x)
  x <- sub("High_NeuMet/Low_MP2",  "Met+/MP2-", x)
  x <- sub("High_NeuMet/High_MP2", "Met+/MP2+", x)
  x
}
df_i$combined_group <- factor(short_lab(as.character(df_i$combined_group)),
                              levels = c("Met-/MP2-", "Met-/MP2+",
                                         "Met+/MP2-", "Met+/MP2+"))
# Ordinal palette: grey (best, ref) → blue → salmon → red (worst)
km_palette <- c(
  "Met-/MP2-" = "#7E7E7E",
  "Met-/MP2+" = "#3C5488",
  "Met+/MP2-" = "#F39B7F",
  "Met+/MP2+" = "#E64B35"
)

fit_i <- survfit(Surv(OS_years, OS_status) ~ combined_group, data = df_i)
# Position p-value in upper-right quadrant (away from curves)
.pval_x <- max(df_i$OS_years, na.rm = TRUE) * 0.62
km <- ggsurvplot(
  fit_i, data = df_i,
  palette = km_palette,
  pval = TRUE, pval.size = 3,
  pval.coord = c(.pval_x, 0.92),
  risk.table = TRUE, risk.table.height = 0.25,
  risk.table.fontsize = 2.5, risk.table.title = "No. at risk",
  risk.table.y.text = FALSE,
  legend = "right", legend.title = "",
  legend.labs = levels(df_i$combined_group),
  xlab = "Time (years)", ylab = "Overall survival",
  ggtheme = theme_pub(8),
  font.family = my_font,
  font.x = c(8, "plain"), font.y = c(8, "plain"),
  font.tickslab = c(7, "plain"), font.legend = c(7, "plain"),
  surv.median.line = "hv",
  break.x.by = 2
)
# survminer sets plot.title by default — clear it explicitly
km$plot <- km$plot + theme(plot.title = element_blank())
km$table <- km$table + theme(plot.title = element_text(size = 7, face = "bold"))

pdf(file.path(OUT, "fig5i_km.pdf"), width = 120/25.4, height = 100/25.4)
print(km)
dev.off()
png(file.path(OUT, "fig5i_km.png"),
    width = 120, height = 100, units = "mm", res = 300)
print(km)
dev.off()

# 5J — Cox forest plot (uni + multi stacked vertically)
cat("[5J] Cox forest plot\n")
df_j <- fread(file.path(DAT, "fig5j_cox_forest_data.csv"))
df_j[, CI_lower_plot := pmax(CI_lower, 0.01)]
df_j[, CI_upper_plot := pmin(CI_upper, 100)]

# Univariate panel (top): all variables, single forest
df_j_uni <- df_j[model == "univariate"]
df_j_uni[, variable := factor(variable, levels = rev(unique(variable[order(-HR)])))]

# Multivariate panel (bottom): pick the headline model NeuInflam+MP1+covar
# (per Axis 1 narrative, this is the publication-grade adjustment)
mv_model_pick <- "NeuInflam+MP1+covar"
df_j_multi <- df_j[model == "multivariate" & covariate == mv_model_pick]
df_j_multi[, variable := factor(variable, levels = rev(unique(variable)))]

# Tag for stacking
df_j_uni[, panel := "Univariate"]
df_j_multi[, panel := paste0("Multivariate (", mv_model_pick, ")")]
df_j_plot <- rbind(df_j_uni, df_j_multi, fill = TRUE)
df_j_plot[, panel := factor(panel,
                            levels = c("Univariate",
                                       paste0("Multivariate (", mv_model_pick, ")")))]

p5j <- ggplot(as.data.frame(df_j_plot), aes(HR, variable, color = effect)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(xmin = CI_lower_plot, xmax = CI_upper_plot),
                width = 0, linewidth = 0.4, orientation = "y") +
  geom_point(size = 2.5, shape = 15) +
  scale_color_manual(values = c("risk" = "#E64B35",
                                "protective" = "#3C5488",
                                "n.s." = "grey50"), name = NULL) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10),
                labels = c("0.1", "0.5", "1", "2", "5", "10")) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_pub(8) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7))

ggsave(file.path(OUT, "fig5j_cox_forest.pdf"), p5j,
       width = 140, height = 130, units = "mm")
ggsave(file.path(OUT, "fig5j_cox_forest.png"), p5j,
       width = 140, height = 130, units = "mm", dpi = 300)

# 5K — Pseudotime trajectory (UMAP colored by DPT + PAGA edges)
cat("[5K] pseudotime trajectory + PAGA overlay\n")
fpt <- file.path(DAT, "fig5k_pseudotime_umap.csv.gz")
fpos <- file.path(DAT, "fig5k_paga_positions.csv")
fconn <- file.path(DAT, "fig5k_paga_connectivity.csv")

if (file.exists(fpt) && file.exists(fpos) && file.exists(fconn)) {
  dpt <- fread(fpt)
  pos <- fread(fpos)
  conn <- fread(fconn)

  # finite pseudotime: scanpy may emit Inf for unreachable cells
  dpt$dpt_pseudotime[!is.finite(dpt$dpt_pseudotime)] <- NA
  pt_cap <- quantile(dpt$dpt_pseudotime, 0.99, na.rm = TRUE)
  dpt$pt_clip <- pmin(dpt$dpt_pseudotime, pt_cap)

  set.seed(42)
  dpt <- dpt[sample(nrow(dpt)), ]

  # Mean pseudotime per subtype → direction (low → high)
  pt_by_sub <- dpt[!is.na(dpt_pseudotime),
                   .(pt_mean = mean(dpt_pseudotime, na.rm = TRUE)),
                   by = neu_subtype]

  # PAGA edges with positions
  edge_df <- merge(conn,
                   pos[, .(source_subtype = subtype, x_src = x, y_src = y)],
                   by = "source_subtype", all.x = TRUE)
  edge_df <- merge(edge_df,
                   pos[, .(target_subtype = subtype, x_tgt = x, y_tgt = y)],
                   by = "target_subtype", all.x = TRUE)
  edge_df <- merge(edge_df,
                   pt_by_sub[, .(source_subtype = neu_subtype, pt_src = pt_mean)],
                   by = "source_subtype", all.x = TRUE)
  edge_df <- merge(edge_df,
                   pt_by_sub[, .(target_subtype = neu_subtype, pt_tgt = pt_mean)],
                   by = "target_subtype", all.x = TRUE)

  # Filter: drop weak edges (bottom 50% by connectivity), keep low-pt → high-pt only
  edge_df <- edge_df[!is.na(connectivity) & !is.na(pt_src) & !is.na(pt_tgt)]
  conn_thresh <- quantile(edge_df$connectivity, 0.5, na.rm = TRUE)
  edge_df <- edge_df[connectivity >= conn_thresh & pt_src < pt_tgt]

  # Shorten edges so arrows land on node boundary, not center
  edge_df[, c("dx", "dy") := .(x_tgt - x_src, y_tgt - y_src)]
  edge_df[, len := sqrt(dx^2 + dy^2)]
  edge_df[, c("ux", "uy") := .(dx / len, dy / len)]
  shrink <- 0.45
  edge_df[, c("x_src2", "y_src2", "x_tgt2", "y_tgt2") := .(
    x_src + ux * shrink,  y_src + uy * shrink,
    x_tgt - ux * shrink,  y_tgt - uy * shrink
  )]

  p5k <- ggplot() +
    ggrastr::rasterise(
      geom_point(data = as.data.frame(dpt),
                 aes(UMAP1, UMAP2, color = pt_clip),
                 size = 0.45, stroke = 0, alpha = 0.55, shape = 16),
      dpi = 600
    ) +
    scale_color_gradientn(
      colours = c("#3C5488", "#4DBBD5", "#F7F7F7", "#F39B7F", "#E64B35"),
      name = "Pseudotime",
      na.value = "#D9D9D9"
    ) +
    # Directed PAGA edges with arrows
    geom_segment(data = as.data.frame(edge_df),
                 aes(x = x_src2, y = y_src2, xend = x_tgt2, yend = y_tgt2,
                     linewidth = connectivity),
                 color = "black", alpha = 0.7,
                 arrow = arrow(length = unit(2, "mm"),
                               type = "closed", angle = 22)) +
    scale_linewidth_continuous(range = c(0.3, 1.4), name = "PAGA conn.") +
    # Nodes
    geom_point(data = as.data.frame(pos),
               aes(x = x, y = y, size = n_cells),
               shape = 21, fill = "white", color = "black", stroke = 0.6) +
    scale_size_continuous(range = c(2.5, 6.5), name = "n cells",
                          guide = "none") +
    # Labels
    geom_text(data = as.data.frame(pos),
              aes(x = x, y = y,
                  label = ifelse(subtype %in% names(neu_labels),
                                 neu_labels[subtype], subtype)),
              size = 2.0, family = my_font, fontface = "bold",
              color = "black", nudge_y = 0.4) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = NULL) +
    theme_pub(8) +
    theme(plot.background = element_blank(), panel.background = element_blank()) +
    umap_arrow_axes(as.data.frame(dpt), "UMAP1", "UMAP2")

  ggsave(file.path(OUT, "fig5k_trajectory.pdf"), p5k,
         width = 120, height = 95, units = "mm")
  ggsave(file.path(OUT, "fig5k_trajectory.png"), p5k,
         width = 120, height = 95, units = "mm", dpi = 300)
} else {
  cat("[5K] missing data files — skip\n")
}

# 5L/5M/5N — GO BP + Hallmark enrichment per anchor subtype
# ── Combined 3-anchor enrichment dotplot (replaces 3 separate bar panels) ──
panels_LMN <- list(
  list(letter = "l", file = "fig5l_enrichment_neu_inflammatory.csv",
       anchor = "Neu_Inflammatory"),
  list(letter = "m", file = "fig5m_enrichment_neu_metastatic.csv",
       anchor = "Neu_Metastatic"),
  list(letter = "n", file = "fig5n_enrichment_neu_ecm_remodeling.csv",
       anchor = "Neu_ECM_remodeling")
)

# Pull top-6 terms per (anchor, library) to define which terms enter the plot
top_per_anchor <- rbindlist(lapply(panels_LMN, function(pl) {
  fp <- file.path(DAT, pl$file)
  if (!file.exists(fp)) return(NULL)
  e <- fread(fp)
  if (nrow(e) == 0) return(NULL)
  e[, anchor := pl$anchor]
  e
}), fill = TRUE)

if (nrow(top_per_anchor) > 0) {
  cat("[5LMN] combined dotplot from", uniqueN(top_per_anchor$anchor),
      "anchors\n")

  top_per_anchor[, lib_short := ifelse(library == "MSigDB_Hallmark_2020",
                                       "Hallmark", "GO BP")]
  top_per_anchor[, lib_short := factor(lib_short, levels = c("GO BP", "Hallmark"))]
  top_per_anchor[, term_short := sub(" \\(GO:[0-9]+\\)", "", term)]
  top_per_anchor[, term_short := sub(" Homo sapiens.*$", "", term_short)]
  top_per_anchor[, term_short := ifelse(nchar(term_short) > 50,
                                        paste0(substr(term_short, 1, 47), "..."),
                                        term_short)]

  # Take top 6 per (anchor, lib) for term set, then pull cross-anchor stats
  top_terms <- top_per_anchor %>%
    group_by(anchor, lib_short) %>%
    arrange(desc(combined_score)) %>%
    slice_head(n = 6) %>%
    ungroup() %>%
    distinct(lib_short, term_short, .keep_all = TRUE) %>%
    arrange(lib_short, anchor, desc(combined_score)) %>%
    pull(term_short)
  top_terms <- unique(top_terms)

  plot_dat <- top_per_anchor[term_short %in% top_terms]
  # if a term is missing in some anchor, fill with NA so dot is invisible
  plot_dat <- plot_dat %>%
    group_by(lib_short, term_short, anchor) %>%
    slice_max(combined_score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.table()

  # Order: terms grouped by lib_short, then by maximum combined_score across anchors.
  # If the same term string appears in two libraries, prefix with lib code so the
  # factor levels stay unique (the prefix is stripped for display).
  plot_dat[, term_uniq := paste0(as.integer(lib_short), "::", term_short)]
  term_order_dt <- plot_dat[, .(score = max(combined_score, na.rm = TRUE)),
                            by = .(lib_short, term_uniq, term_short)]
  term_order_dt <- term_order_dt[order(lib_short, score)]
  plot_dat[, term_uniq := factor(term_uniq, levels = term_order_dt$term_uniq)]
  # display labels (strip prefix)
  plot_dat[, term_label := sub("^[0-9]+::", "", as.character(term_uniq))]

  # Anchor display order + display labels
  anchor_levels <- c("Neu_Inflammatory", "Neu_Metastatic", "Neu_ECM_remodeling")
  anchor_disp <- neu_labels[anchor_levels]
  plot_dat[, anchor_disp := factor(neu_labels[anchor], levels = anchor_disp)]

  # Map factor → display label for y-axis
  y_break_lvls <- levels(plot_dat$term_uniq)
  y_break_labs <- sub("^[0-9]+::", "", y_break_lvls)

  p_lmn <- ggplot(as.data.frame(plot_dat),
                  aes(x = anchor_disp, y = term_uniq)) +
    geom_point(aes(size = minus_log10_padj, fill = combined_score),
               shape = 21, color = "black", stroke = 0.25) +
    scale_size_continuous(name = expression(-log[10]~italic(P)["adj"]),
                          range = c(1, 6),
                          breaks = pretty(plot_dat$minus_log10_padj, 4)) +
    scale_fill_gradientn(name = "Combined\nscore",
                         colours = c("#F7F7F7", "#FEE0D2", "#FCAE91",
                                     "#FB6A4A", "#CB181D"),
                         na.value = "white",
                         guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                barheight = unit(15, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.3,
                                                ticks.colour = "black")) +
    scale_y_discrete(breaks = y_break_lvls, labels = y_break_labs) +
    facet_grid(lib_short ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_pub(7) +
    theme(axis.text.y = element_text(size = 6, lineheight = 0.85),
          axis.text.x = element_text(size = 7, face = "italic", angle = 20,
                                     hjust = 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.3),
          strip.text.y = element_text(face = "bold", size = 7, angle = -90),
          legend.box = "vertical",
          legend.spacing.y = unit(2, "mm"))

  ggsave(file.path(OUT, "fig5lmn_enrichment.pdf"), p_lmn,
         width = 150, height = 165, units = "mm")
  ggsave(file.path(OUT, "fig5lmn_enrichment.png"), p_lmn,
         width = 150, height = 165, units = "mm", dpi = 300)

  # remove obsolete single-panel bar files
  for (pl in panels_LMN) {
    for (ext in c("pdf", "png")) {
      f <- file.path(OUT, sprintf("fig5%s_enrichment.%s", pl$letter, ext))
      if (file.exists(f)) file.remove(f)
    }
  }
  cat("  fig5lmn (3-anchor combined dotplot) saved\n")
} else {
  cat("[5LMN] no enrichment data — skip\n")
}

cat("\n=== DONE — outputs in", OUT, "===\n")
print(list.files(OUT))
