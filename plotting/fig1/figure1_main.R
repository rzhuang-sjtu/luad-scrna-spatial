# Figure 1: LUAD
# Panel layout adapted from Fig. 1 of the HCC study cited in the manuscript
# Arial | : |

required_cran <- c("ggplot2", "dplyr", "tidyr", "data.table", "R.utils",
                   "patchwork", "RColorBrewer", "scales", "ggnewscale",
                   "ggrepel", "showtext", "ragg")
required_bioc <- c("ComplexHeatmap", "circlize")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(R.utils)
library(patchwork)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(showtext)
library(ragg)
library(ComplexHeatmap)
library(circlize)
library(grid)

font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf")
showtext_auto()
showtext_opts(dpi = 300)

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig1" else
  "${WORK_ROOT}/luad_figures/fig1")

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

# ── ( Fig.1C ) ──
ct_colors <- c(
  "Epithelial"  = "#E64B35FF",
  "Endothelial" = "#4DBBD5FF",
  "Fibroblast"  = "#00A087FF",
  "Myeloid"     = "#3C5488FF",
  "Mast"        = "#F39B7FFF",
  "T_NK"        = "#8491B4FF",
  "B"           = "#91D1C2FF",
  "Plasma"      = "#B09C85FF"
)

ct_order <- c("Epithelial", "Endothelial", "Fibroblast",
              "Myeloid", "Mast", "T_NK", "B", "Plasma")

tissue_order <- c("Normal_Lung", "Adjacent_Normal", "Normal_LN",
                  "Precancerous", "Primary_Tumor", "LN_Metastasis",
                  "Brain_Metastasis", "Distant_Metastasis", "Pleural_Effusion")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") %+replace%
    theme(
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "pt"),
      axis.text = element_text(color = "black", size = rel(1)),
      axis.title = element_text(size = rel(1.1)),
      legend.title = element_text(size = rel(1), face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      legend.spacing.y = unit(1, "pt"),
      plot.margin = margin(4, 4, 4, 4, "pt"),
      plot.title = element_text(size = rel(1.3), face = "bold", hjust = 0),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = rel(1), face = "bold")
    )
}


# Fig 1B: Dot Plot + +

cat("── Fig 1B: Dot Plot ──\n")

dot <- fread("dotplot_markers.csv")

# y x marker →
ct_order_stair <- c("Epithelial", "Endothelial", "Fibroblast",
                    "Myeloid", "Mast", "T_NK", "B", "Plasma")

# marker cell type
marker_block_order <- c(
  "Epithelial", "Epithelial_prolif", "Supplementary",
  "Endothelial", "Fibroblast",
  "Myeloid", "Mast", "T_NK", "B", "Plasma"
)

dot$marker_group <- factor(dot$marker_group, levels = marker_block_order)
dot <- dot[order(dot$marker_group), ]
gene_order <- unique(dot$gene)

# x y cell type
dot$gene <- factor(dot$gene, levels = gene_order)
dot$celltype <- factor(dot$celltype, levels = ct_order_stair)

p1b <- ggplot(dot, aes(x = gene, y = celltype)) +
  geom_point(aes(size = frac_expressed * 100, fill = mean_expression),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(
    name = "Fraction of cells\nin group (%)",
    range = c(0.3, 5.5),
    breaks = c(25, 50, 75),
    limits = c(0, 100)
  ) +
  scale_fill_gradient2(
    name = "Mean expression\nin group",
    low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
    midpoint = 1,
    limits = c(0, 3),
    oob = scales::squish,
    guide = guide_colorbar(barwidth = unit(3.5, "mm"), barheight = unit(20, "mm"),
                           frame.colour = "black", frame.linewidth = 0.3,
                           ticks.colour = "black")
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                               size = 6, face = "italic"),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(4, "mm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.3)
  )

ggsave("Fig1B_dotplot.pdf", p1b, width = 180, height = 55, units = "mm",
       device = cairo_pdf)
ggsave("Fig1B_dotplot.png", p1b, width = 180, height = 55, units = "mm",
       dpi = 300, device = ragg::agg_png)
cat("  Fig1B saved\n")

# Fig 1C: UMAP scatter

cat("── Fig 1C: UMAP ──\n")

meta <- fread("cell_metadata.csv.gz")
meta$celltype_coarse <- factor(meta$celltype_coarse, levels = ct_order)

set.seed(42)
meta <- meta[sample(.N), ]
ct_freq <- meta[, .N, by = celltype_coarse][order(-N)]
meta$celltype_coarse <- factor(meta$celltype_coarse,
                               levels = ct_freq$celltype_coarse)
meta <- meta[order(meta$celltype_coarse, decreasing = TRUE), ]

# 30 85 PDF
# 30 85 PDF
set.seed(42)
if (nrow(meta) > 300000) {
  total_n <- nrow(meta)
  meta_plot <- meta[, {
    target <- min(.N, as.integer(ceiling(300000 * .N / total_n)) + 500L)
    .SD[sample(.N, target)]
  }, by = celltype_coarse]
} else {
  meta_plot <- meta
}

p1c <- ggplot(meta_plot, aes(x = UMAP_1, y = UMAP_2, color = celltype_coarse)) +
  geom_point(size = 0.25, alpha = 0.55, stroke = 0, shape = 16) +
  scale_color_manual(
    values = ct_colors,
    name = "celltype",
    guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  coord_fixed(ratio = 1) +
  umap_arrow_axes(meta_plot, "UMAP_1", "UMAP_2")

ggsave("Fig1C_umap.pdf", p1c, width = 100, height = 85, units = "mm",
       device = cairo_pdf)
ggsave("Fig1C_umap.png", p1c, width = 100, height = 85, units = "mm",
       dpi = 300, device = ragg::agg_png)
cat("  Fig1C saved\n")


# Fig 1D: Enrichment Heatmap Fig.1D

cat("── Fig 1D: Enrichment Heatmap ──\n")

enrich <- read.csv("enrichment_heatmap.csv", row.names = 1, check.names = FALSE)
enrich_p <- read.csv("enrichment_pvalues.csv", row.names = 1, check.names = FALSE)

tissues_avail <- intersect(tissue_order, rownames(enrich))
ct_avail <- intersect(ct_order, colnames(enrich))
mat_e <- as.matrix(enrich[tissues_avail, ct_avail])
mat_p <- as.matrix(enrich_p[tissues_avail, ct_avail])

sig_marks <- matrix("", nrow(mat_p), ncol(mat_p))
sig_marks[mat_p < 0.0001] <- "****"
sig_marks[mat_p >= 0.0001 & mat_p < 0.001] <- "***"
sig_marks[mat_p >= 0.001 & mat_p < 0.01] <- "**"
sig_marks[mat_p >= 0.01 & mat_p < 0.05] <- "*"

mat_e_clamp <- pmax(pmin(mat_e, 2), -2)

col_fun <- colorRamp2(
  seq(-2, 2, length.out = 9),
  rev(brewer.pal(9, "RdBu"))
)

tissue_group <- ifelse(grepl("Normal|Adjacent", tissues_avail), "Normal",
                       ifelse(grepl("Precancerous", tissues_avail), "Precancerous",
                              ifelse(grepl("Primary", tissues_avail), "Primary", "Metastatic")))
tissue_group_colors <- c(
  "Normal" = "#4DAF4A", "Precancerous" = "#FF7F00",
  "Primary" = "#E41A1C", "Metastatic" = "#984EA3"
)

ha_row <- rowAnnotation(
  Stage = tissue_group,
  col = list(Stage = tissue_group_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Stage = list(title_gp = gpar(fontsize = 7, fontfamily = "Arial", fontface = "bold"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Arial"))
  ),
  simple_anno_size = unit(3, "mm")
)

# celltype
ha_col <- HeatmapAnnotation(
  CellType = ct_avail,
  col = list(CellType = ct_colors[ct_avail]),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(3, "mm")
)

ht_d <- Heatmap(
  mat_e_clamp,
  name = "log2(O/E)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_marks[i, j], x, y,
              gp = gpar(fontsize = 5, fontfamily = "Arial", col = "black"))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "Arial", fontface = "bold"),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    legend_height = unit(25, "mm"),
    legend_width = unit(3, "mm"),
    grid_width = unit(3, "mm")
  ),
  width = unit(ncol(mat_e_clamp) * 8, "mm"),
  height = unit(nrow(mat_e_clamp) * 6, "mm")
)

pdf("Fig1D_enrichment_heatmap.pdf", width = 7, height = 4)
ComplexHeatmap::draw(ht_d, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()

png("Fig1D_enrichment_heatmap.png", width = 7, height = 4, units = "in", res = 300)
ComplexHeatmap::draw(ht_d, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()
cat("  Fig1D saved\n")


# Fig 1E: Fig.1E

cat("── Fig 1E: Stacked Bar ──\n")

prop <- fread("proportion_by_tissue.csv")
prop <- prop[tissue_type %in% tissue_order]
prop$tissue_type <- factor(prop$tissue_type, levels = tissue_order)
prop$celltype <- factor(prop$celltype, levels = rev(ct_order))

tissue_labels <- c(
  "Normal_Lung" = "Normal", "Adjacent_Normal" = "Adjacent",
  "Normal_LN" = "Normal LN", "Precancerous" = "Precancerous",
  "Primary_Tumor" = "Tumor", "LN_Metastasis" = "LN Met",
  "Brain_Metastasis" = "Brain Met", "Distant_Metastasis" = "Distant Met",
  "Pleural_Effusion" = "Pleural Effusion"
)

p1e <- ggplot(prop, aes(x = tissue_type, y = percent / 100, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.75, linewidth = 0.1, color = "white") +
  scale_fill_manual(values = ct_colors, name = "Cell Type",
                    breaks = ct_order) +
  scale_x_discrete(labels = tissue_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cell percent ratio") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.key.size = unit(3, "mm"),
    # Force the plotting panel itself to be square regardless of canvas
    aspect.ratio = 1
  )

ggsave("Fig1E_stacked_bar.pdf", p1e, width = 110, height = 90, units = "mm",
       device = cairo_pdf)
ggsave("Fig1E_stacked_bar.png", p1e, width = 110, height = 90, units = "mm",
       dpi = 300, device = ragg::agg_png)
cat("  Fig1E saved\n")


# Fig 1F: Fig.1F

cat("── Fig 1F: Correlation Heatmap ──\n")

corr_mat <- read.csv("correlation_matrix.csv", row.names = 1, check.names = FALSE)
corr_p <- read.csv("correlation_pvalues.csv", row.names = 1, check.names = FALSE)

ct_avail_f <- intersect(ct_order, rownames(corr_mat))
mat_c <- as.matrix(corr_mat[ct_avail_f, ct_avail_f])
mat_cp <- as.matrix(corr_p[ct_avail_f, ct_avail_f])

sig_f <- matrix("", nrow(mat_cp), ncol(mat_cp))
sig_f[mat_cp < 0.0001] <- "****"
sig_f[mat_cp >= 0.0001 & mat_cp < 0.001] <- "***"
sig_f[mat_cp >= 0.001 & mat_cp < 0.01] <- "**"
sig_f[mat_cp >= 0.01 & mat_cp < 0.05] <- "*"
diag(sig_f) <- ""

col_fun_corr <- colorRamp2(
  seq(-1, 1, length.out = 11),
  rev(brewer.pal(11, "RdBu"))
)

# CellType color strips on both axes matches 1D's top-annotation pattern
# anchors each row/column to the celltype palette used elsewhere in Fig 1.
ha_col_f <- HeatmapAnnotation(
  CellType = ct_avail_f,
  col = list(CellType = ct_colors[ct_avail_f]),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(3, "mm")
)
ha_row_f <- rowAnnotation(
  CellType = ct_avail_f,
  col = list(CellType = ct_colors[ct_avail_f]),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(3, "mm")
)

ht_f <- Heatmap(
  mat_c,
  name = "Correlation",
  col = col_fun_corr,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha_col_f,
  left_annotation = ha_row_f,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- mat_c[i, j]
    mark <- sig_f[i, j]
    grid.text(
      sprintf("%.1f", val), x, y - unit(0.5, "mm"),
      gp = gpar(fontsize = 5.5, fontfamily = "Arial",
                col = ifelse(abs(val) > 0.6, "white", "black"))
    )
    if (nchar(mark) > 0) {
      grid.text(
        mark, x, y + unit(2, "mm"),
        gp = gpar(fontsize = 4, fontfamily = "Arial",
                  col = ifelse(abs(val) > 0.6, "white", "black"))
      )
    }
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = "Arial", fontface = "bold"),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    legend_height = unit(25, "mm"),
    grid_width = unit(3, "mm")
  ),
  width = unit(ncol(mat_c) * 9, "mm"),
  height = unit(nrow(mat_c) * 9, "mm")
)

pdf("Fig1F_correlation_heatmap.pdf", width = 6, height = 4.5)
ComplexHeatmap::draw(ht_f, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()

png("Fig1F_correlation_heatmap.png", width = 6, height = 4.5, units = "in", res = 300)
ComplexHeatmap::draw(ht_f, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()
cat("  Fig1F saved\n")


cat("\n══════════════════════════════\n")
cat("  \n")
for (f in list.files(pattern = "^Fig1.*\\.(pdf|png)$")) {
  sz <- round(file.info(f)$size / 1e6, 1)
  cat(sprintf("  %s (%.1f MB)\n", f, sz))
}
cat("══════════════════════════════\n")