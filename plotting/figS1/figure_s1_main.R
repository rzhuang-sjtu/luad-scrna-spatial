# Supplementary Figure 1: QC, Batch Correction, Cell Type Validation
# Panel layout adapted from Fig. S1 of the HCC study cited in the manuscript

required_cran <- c("ggplot2", "dplyr", "data.table", "R.utils",
                   "patchwork", "RColorBrewer", "scales",
                   "showtext", "ragg", "ggrastr", "tidyr")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(data.table)
library(R.utils)
library(patchwork)
library(RColorBrewer)
library(scales)
library(showtext)
library(ragg)
library(ggrastr)
library(tidyr)

font_add("Arial",
         regular = "C:/Windows/Fonts/arial.ttf",
         bold    = "C:/Windows/Fonts/arialbd.ttf",
         italic  = "C:/Windows/Fonts/ariali.ttf")
showtext_auto()
showtext_opts(dpi = 300)

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig_s1" else
  "${WORK_ROOT}/luad_figures/fig_s1")

my_font <- "Arial"
# ── UMAP corner-arrow helper ──
umap_arrow_axes <- function(data, x_col, y_col,
                            frac = 0.16, label_x = "UMAP 1", label_y = "UMAP 2",
                            text_size = 2.0, line_size = 0.4, arrow_mm = 1.5,
                            inset_frac = 0.02) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr); y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr);          y1 <- y0 + frac * diff(yr)
  arr <- grid::arrow(length = grid::unit(arrow_mm, "mm"), ends = "last", type = "closed")
  list(
    annotate("segment", x=x0, xend=x1, y=y0, yend=y0, arrow=arr, linewidth=line_size, color="black"),
    annotate("segment", x=x0, xend=x0, y=y0, yend=y1, arrow=arr, linewidth=line_size, color="black"),
    annotate("text", x=(x0+x1)/2, y=y0, label=label_x, vjust=2.4, size=text_size, family=my_font),
    annotate("text", x=x0, y=(y0+y1)/2, label=label_y, angle=90, vjust=-1.4, size=text_size, family=my_font),
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          panel.grid=element_blank())
  )
}

ct_colors <- c(
  "Epithelial"  = "#E64B35FF", "Endothelial" = "#4DBBD5FF",
  "Fibroblast"  = "#00A087FF", "Myeloid"     = "#3C5488FF",
  "Mast"        = "#F39B7FFF", "T_NK"        = "#8491B4FF",
  "B"           = "#91D1C2FF", "Plasma"      = "#B09C85FF"
)
ct_order <- c("Epithelial", "Endothelial", "Fibroblast",
              "Myeloid", "Mast", "T_NK", "B", "Plasma")

tissue_group_colors <- c(
  "Normal" = "#4DAF4A", "Precancerous" = "#FF7F00",
  "Primary" = "#E41A1C", "Metastatic" = "#984EA3"
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") %+replace%
    theme(
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(1.5, "pt"),
      axis.text         = element_text(color = "black", size = rel(1)),
      axis.title        = element_text(size = rel(1.1)),
      legend.title      = element_text(size = rel(1), face = "bold"),
      legend.text       = element_text(size = rel(0.85)),
      legend.key.size   = unit(2.5, "mm"),
      legend.background = element_blank(),
      plot.margin       = margin(4, 4, 4, 4, "pt"),
      plot.title        = element_text(size = 9, face = "bold", hjust = 0,
                                       family = "Arial"),
      panel.border      = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(size = 6, face = "bold")
    )
}

cat("...\n")
meta <- fread("s1_cell_metadata.csv.gz")
cat(sprintf("  %s cells\n", format(nrow(meta), big.mark = ",")))

meta[, tissue_group := fifelse(
  tissue_type %in% c("Normal_Lung", "Adjacent_Normal", "Normal_LN"), "Normal",
  fifelse(tissue_type == "Precancerous", "Precancerous",
          fifelse(tissue_type == "Primary_Tumor", "Primary", "Metastatic"))
)]
meta$tissue_group <- factor(meta$tissue_group,
                            levels = c("Normal", "Precancerous", "Primary", "Metastatic"))

# dataset sample
n_samples <- meta[, uniqueN(sample_id)]
cat(sprintf("  %d samples across %d datasets\n", n_samples, meta[, uniqueN(dataset)]))

# Panel A+B+C: GSE C

cat("── Panel A+B+C: QC per sample ──\n")

set.seed(42)
meta_qc <- meta[, .SD[sample(.N, min(.N, 2000))], by = sample_id]
sample_order <- meta_qc[, .(med = median(n_counts)), by = .(dataset, sample_id)
][order(dataset, med)]$sample_id
meta_qc$sample_id <- factor(meta_qc$sample_id, levels = sample_order)

ds_pal <- brewer.pal(min(7, meta[, uniqueN(dataset)]), "Set2")
ds_names <- sort(unique(meta$dataset))
names(ds_pal) <- ds_names

pA <- ggplot(meta_qc, aes(x = sample_id, y = n_counts, fill = dataset)) +
  rasterise(geom_violin(scale = "width", linewidth = 0, alpha = 0.7, width = 0.9), dpi = 300) +
  geom_boxplot(width = 0.3, outlier.shape = NA, linewidth = 0.15,
               fill = "white", alpha = 0.6, coef = 0) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = ds_pal) +
  facet_grid(~ dataset, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "nCount_RNA") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1, "pt"),
    strip.text = element_blank()
  )

pB <- ggplot(meta_qc, aes(x = sample_id, y = n_genes, fill = dataset)) +
  rasterise(geom_violin(scale = "width", linewidth = 0, alpha = 0.7, width = 0.9), dpi = 300) +
  geom_boxplot(width = 0.3, outlier.shape = NA, linewidth = 0.15,
               fill = "white", alpha = 0.6, coef = 0) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = ds_pal) +
  facet_grid(~ dataset, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "nFeature_RNA") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1, "pt"),
    strip.text = element_blank()
  )

pC <- ggplot(meta_qc, aes(x = sample_id, y = pct_mito, fill = dataset)) +
  rasterise(geom_violin(scale = "width", linewidth = 0, alpha = 0.7, width = 0.9), dpi = 300) +
  geom_boxplot(width = 0.3, outlier.shape = NA, linewidth = 0.15,
               fill = "white", alpha = 0.6, coef = 0) +
  scale_fill_manual(values = ds_pal, name = "Dataset") +
  scale_y_continuous(limits = c(0, 30), oob = scales::squish) +
  facet_grid(~ dataset, scales = "free_x", space = "free_x") +
  labs(x = "Samples (grouped by dataset)", y = "percent.mt (%)") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "pt"),
    strip.text = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm")
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))

pABC <- pA / pB / pC + plot_layout(heights = c(1, 1, 1.15)) &
        theme(plot.margin = margin(1, 2, 1, 2))

ggsave("FigS1_ABC.pdf", pABC, width = 183, height = 90, units = "mm", device = cairo_pdf)
ggsave("FigS1_ABC.png", pABC, width = 183, height = 90, units = "mm",
       dpi = 300, device = ragg::agg_png)
cat("  FigS1_ABC saved\n")

# Panel D: QC by tissue group
# Normal / Tumor / PVTT / MLN → LUAD: Normal / Precancerous / Primary / Metastatic

cat("── Panel D: QC by tissue group ──\n")

set.seed(42)
meta_d <- meta[sample(.N, min(.N, 200000))]

pD1 <- ggplot(meta_d, aes(x = tissue_group, y = n_counts, fill = tissue_group)) +
  rasterise(geom_violin(scale = "width", linewidth = 0.2, alpha = 0.7), dpi = 300) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, linewidth = 0.2,
               fill = "white", alpha = 0.8) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = tissue_group_colors) +
  labs(x = NULL, y = "nCount_RNA") +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
        legend.position = "none")

pD2 <- ggplot(meta_d, aes(x = tissue_group, y = n_genes, fill = tissue_group)) +
  rasterise(geom_violin(scale = "width", linewidth = 0.2, alpha = 0.7), dpi = 300) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, linewidth = 0.2,
               fill = "white", alpha = 0.8) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = tissue_group_colors) +
  labs(x = NULL, y = "nFeature_RNA") +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
        legend.position = "none")

pD3 <- ggplot(meta_d, aes(x = tissue_group, y = pct_mito, fill = tissue_group)) +
  rasterise(geom_violin(scale = "width", linewidth = 0.2, alpha = 0.7), dpi = 300) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, linewidth = 0.2,
               fill = "white", alpha = 0.8) +
  scale_fill_manual(values = tissue_group_colors) +
  labs(x = NULL, y = "percent.mt (%)") +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
        legend.position = "none")

pD <- (pD1 | pD2 | pD3)


# Panel E: Leiden clusters UMAP
# 33 clusters, UMAP

cat("── Panel E: Clusters UMAP ──\n")

n_clusters <- meta[, uniqueN(`leiden_1.0`)]
cluster_pal <- colorRampPalette(
  c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(9, "Pastel1"))
)(n_clusters)
cluster_ids_sorted <- as.character(sort(as.numeric(unique(meta$`leiden_1.0`))))
names(cluster_pal) <- cluster_ids_sorted

# UMAP
set.seed(42)
total_n <- nrow(meta)
meta_umap <- meta[, {
  target <- min(.N, as.integer(ceiling(300000 * .N / total_n)) + 500L)
  .SD[sample(.N, target)]
}, by = celltype_coarse]

pE <- ggplot(meta_umap, aes(x = UMAP_1, y = UMAP_2,
                            color = factor(`leiden_1.0`, levels = cluster_ids_sorted))) +
  rasterise(geom_point(size = 0.25, alpha = 0.55, stroke = 0, shape = 16), dpi = 300) +
  scale_color_manual(values = cluster_pal, name = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5)) +
  coord_fixed() +
  umap_arrow_axes(meta_umap, "UMAP_1", "UMAP_2")


# Panel F: Before / After Harmony UMAP

cat("── Panel F: Harmony before/after ──\n")

cat("  meta_umap :", paste(names(meta_umap), collapse = ", "), "\n")
cat("   UMAP_pre_1:", "UMAP_pre_1" %in% names(meta_umap), "\n")

pre1_col <- grep("UMAP_pre.*1", names(meta_umap), value = TRUE)[1]
pre2_col <- grep("UMAP_pre.*2", names(meta_umap), value = TRUE)[1]

if (!is.na(pre1_col) && !is.na(pre2_col)) {
  setnames(meta_umap, c(pre1_col, pre2_col), c("UMAP_pre_1", "UMAP_pre_2"))
}

if ("UMAP_pre_1" %in% names(meta_umap)) {
  
  # geom_point_rast() both failed silently at this aspect ratio with
  # ~300k points). Sub-sample to 80k cells per panel and use plain
  # vector geom_point file size stays acceptable (~1-2 MB) and
  # rendering is reliable.
  set.seed(7)
  meta_pf <- meta_umap[sample(.N, min(.N, 80000))]
  pF_pre <- ggplot(meta_pf, aes(x = UMAP_pre_1, y = UMAP_pre_2, color = dataset)) +
    geom_point(size = 0.18, alpha = 0.5, stroke = 0, shape = 16) +
    scale_color_manual(values = ds_pal) +
    labs(x = NULL, y = NULL, subtitle = "Before Harmony") +
    theme_pub(base_size = 7) +
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5, size = 7, face = "italic")) +
    coord_fixed() +
    umap_arrow_axes(meta_pf, "UMAP_pre_1", "UMAP_pre_2")

  pF_post <- ggplot(meta_pf, aes(x = UMAP_1, y = UMAP_2, color = dataset)) +
    geom_point(size = 0.18, alpha = 0.5, stroke = 0, shape = 16) +
    scale_color_manual(values = ds_pal, name = "Dataset") +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = NULL, y = NULL, subtitle = "After Harmony") +
    theme_pub(base_size = 7) +
    theme(legend.position = "right",
          legend.key.size = unit(2.5, "mm"),
          plot.subtitle = element_text(hjust = 0.5, size = 7, face = "italic")) +
    coord_fixed() +
    umap_arrow_axes(meta_pf, "UMAP_1", "UMAP_2")

  pF <- (pF_pre | pF_post)

} else {
  cat("  ️ UMAP_pre  :", paste(names(meta_umap), collapse = ", "), "\n")
  pF <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                            label = "Panel F: UMAP_pre columns not found", size = 3) +
    theme_void()
}
# Panel G: LISI ( kBET
# kBET rejection rate LISI violin

cat("── Panel G: LISI ──\n")

lisi <- fread("s1_lisi_per_cell.csv.gz")
lisi_long <- data.table(
  LISI  = c(lisi$LISI_pre, lisi$LISI_post),
  stage = rep(c("Before Harmony", "After Harmony"), each = nrow(lisi))
)
lisi_long$stage <- factor(lisi_long$stage, levels = c("Before Harmony", "After Harmony"))

set.seed(42)
lisi_sub <- lisi_long[sample(.N, min(.N, 100000))]

lisi_summary <- fread("s1_lisi_scores.csv")

pG <- ggplot(lisi_sub, aes(x = stage, y = LISI, fill = stage)) +
  rasterise(geom_violin(linewidth = 0.2, alpha = 0.7, scale = "width"), dpi = 300) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, linewidth = 0.2,
               fill = "white", alpha = 0.8) +
  scale_fill_manual(values = c("Before Harmony" = "#D73027", "After Harmony" = "#4575B4")) +
  geom_hline(yintercept = 7, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  annotate("text", x = 2.4, y = 6.7, label = "Perfect mixing (7 datasets)",
           size = 2, color = "grey40", family = "Arial", hjust = 1) +
  annotate("text", x = 1, y = 5.5,
           label = sprintf("median = %.2f", lisi_summary$median_LISI[1]),
           size = 2, family = "Arial", fontface = "bold", color = "#D73027") +
  annotate("text", x = 2, y = 5.5,
           label = sprintf("median = %.2f", lisi_summary$median_LISI[2]),
           size = 2, family = "Arial", fontface = "bold", color = "#4575B4") +
  scale_y_continuous(limits = c(0.8, 7.5), breaks = 1:7) +
  labs(x = NULL, y = "LISI (dataset integration)") +
  theme_pub(base_size = 7) +
  theme(legend.position = "none")
# Panel H: Cluster × sample facet

cat("── Panel H: Cluster × sample ──\n")

cluster_sample <- fread("s1_cluster_by_sample.csv", header = FALSE)
real_header <- as.character(cluster_sample[1, ])
cluster_sample <- cluster_sample[-1, ]
setnames(cluster_sample, real_header)
cluster_sample[, (real_header[-1]) := lapply(.SD, as.numeric), .SDcols = real_header[-1]]

sample_ds <- unique(meta[, .(sample_id, dataset)])
cluster_sample <- merge(cluster_sample, sample_ds, by = "sample_id")

cluster_long <- melt(cluster_sample, id.vars = c("sample_id", "dataset"),
                     variable.name = "cluster", value.name = "proportion")
cluster_long <- cluster_long[!is.na(proportion)]

actual_clusters <- real_header[-1]
n_cl <- length(actual_clusters)
cluster_pal_h <- colorRampPalette(
  c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(9, "Pastel1"))
)(n_cl)
names(cluster_pal_h) <- actual_clusters
cluster_long$cluster <- factor(cluster_long$cluster, levels = actual_clusters)

sample_order_h <- cluster_sample[order(dataset, sample_id)]$sample_id
cluster_long$sample_id <- factor(cluster_long$sample_id, levels = sample_order_h)

pH_main <- ggplot(cluster_long, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, linewidth = 0) +
  scale_fill_manual(values = cluster_pal_h, name = "Cluster", drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Proportion") +
  theme_pub(base_size = 7) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(size = 4),
    plot.margin = margin(4, 4, 0, 4, "pt")
  ) +
  guides(fill = guide_legend(ncol = 2))

ds_strip <- data.table(
  sample_id = factor(sample_order_h, levels = sample_order_h),
  dataset = cluster_sample[order(dataset, sample_id)]$dataset
)

pH_strip <- ggplot(ds_strip, aes(x = sample_id, y = 1, fill = dataset)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = ds_pal, name = "Dataset") +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6, face = "bold"),
    plot.margin = margin(0, 4, 4, 4, "pt")
  ) +
  guides(fill = guide_legend(nrow = 1))

pH <- pH_main / pH_strip + plot_layout(heights = c(12, 1))

ggsave("FigS1_H.pdf", pH, width = 183, height = 70, units = "mm", device = cairo_pdf)
ggsave("FigS1_H.png", pH, width = 183, height = 70, units = "mm",
       dpi = 300, device = ragg::agg_png)
cat("  FigS1_H saved\n")
# Panel I: Before / After regressing out nCount UMAP

cat("── Panel I: Regress comparison ──\n")

has_regress <- file.exists("s1_regress_umap.csv.gz")

if (has_regress) {
  regress <- fread("s1_regress_umap.csv.gz")
  setnames(regress, 1, "cell_id")
  
  # meta_umap
  cat("  meta_umap :", paste(names(meta_umap), collapse = ", "), "\n")
  
  # barcode barcode +
  known_non_barcode <- c("celltype_coarse", "dataset", "sample_id", "patient_id",
                         "tissue_type", "tissue_group", "celltype_confidence",
                         "celltype_celltypist", "celltype_ct_coarse",
                         "celltype_original_mapped",
                         "leiden_1.0", "leiden_2.0",
                         "UMAP_1", "UMAP_2", "UMAP_pre_1", "UMAP_pre_2",
                         "n_counts", "n_genes", "pct_mito", "pct_ribo",
                         "doublet_score")
  barcode_col_meta <- NULL
  for (col in names(meta_umap)) {
    if (col %in% known_non_barcode) next
    val <- as.character(meta_umap[[col]][1])
    if (nchar(val) > 10 && grepl("GSE|\\d{6,}", val)) {
      barcode_col_meta <- col
      break
    }
  }
  
  if (is.null(barcode_col_meta)) {
    if ("V1" %in% names(meta_umap)) barcode_col_meta <- "V1"
  }
  
  if (!is.null(barcode_col_meta)) {
    cat("   barcode :", barcode_col_meta, "\n")
    cat("  barcode :", head(meta_umap[[barcode_col_meta]], 3), "\n")
    cat("  regress :", head(regress$cell_id, 3), "\n")
    
    regress_idx <- match(meta_umap[[barcode_col_meta]], regress$cell_id)
    n_matched <- sum(!is.na(regress_idx))
    cat(sprintf("  : %d / %d\n", n_matched, nrow(meta_umap)))
    
    meta_umap_i <- copy(meta_umap)
    meta_umap_i$UMAP_regress_1 <- regress$UMAP_regress_1[regress_idx]
    meta_umap_i$UMAP_regress_2 <- regress$UMAP_regress_2[regress_idx]
    
  } else {
    # meta merge
    cat("  ️  meta_umap  barcode   meta \n")
    
    meta_full <- fread("s1_cell_metadata.csv.gz")
    setnames(meta_full, 1, "cell_id")
    
    total_n_full <- nrow(meta_full)
    set.seed(42)
    meta_sub_i <- meta_full[, {
      target <- min(.N, as.integer(ceiling(300000 * .N / total_n_full)) + 500L)
      .SD[sample(.N, target)]
    }, by = celltype_coarse]
    
    regress_idx <- match(meta_sub_i$cell_id, regress$cell_id)
    n_matched <- sum(!is.na(regress_idx))
    cat(sprintf("   meta : %d / %d\n", n_matched, nrow(meta_sub_i)))
    
    meta_umap_i <- copy(meta_sub_i)
    meta_umap_i$UMAP_regress_1 <- regress$UMAP_regress_1[regress_idx]
    meta_umap_i$UMAP_regress_2 <- regress$UMAP_regress_2[regress_idx]
  }
  
  n_final <- sum(!is.na(meta_umap_i$UMAP_regress_1))
  cat(sprintf("  : %d cells\n", n_final))
  
  if (n_final > 1000) {
    meta_umap_i$celltype_coarse <- factor(meta_umap_i$celltype_coarse, levels = ct_order)
    
    ct_freq_i <- meta_umap_i[!is.na(UMAP_regress_1), .N, by = celltype_coarse][order(-N)]
    meta_umap_i <- meta_umap_i[order(match(celltype_coarse, ct_freq_i$celltype_coarse),
                                     decreasing = TRUE)]
    
    pI_before <- ggplot(meta_umap_i, aes(x = UMAP_1, y = UMAP_2, color = celltype_coarse)) +
      rasterise(geom_point(size = 0.25, alpha = 0.55, stroke = 0, shape = 16), dpi = 300) +
      scale_color_manual(values = ct_colors, breaks = ct_order) +
      labs(x = NULL, y = NULL, subtitle = "Before regressing out nCount") +
      theme_pub(base_size = 7) +
      theme(legend.position = "none",
            plot.subtitle = element_text(hjust = 0.5, size = 6, face = "italic")) +
      coord_fixed() +
      umap_arrow_axes(meta_umap_i, "UMAP_1", "UMAP_2")

    pI_after <- ggplot(meta_umap_i[!is.na(UMAP_regress_1)],
                       aes(x = UMAP_regress_1, y = UMAP_regress_2, color = celltype_coarse)) +
      rasterise(geom_point(size = 0.25, alpha = 0.55, stroke = 0, shape = 16), dpi = 300) +
      scale_color_manual(values = ct_colors, name = "Cell Type", breaks = ct_order) +
      guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
      labs(x = NULL, y = NULL, subtitle = "After regressing out nCount") +
      theme_pub(base_size = 7) +
      theme(legend.position = "right", legend.key.size = unit(3, "mm"),
            plot.subtitle = element_text(hjust = 0.5, size = 6, face = "italic")) +
      coord_fixed() +
      umap_arrow_axes(meta_umap_i[!is.na(UMAP_regress_1)], "UMAP_regress_1", "UMAP_regress_2")
    
    pI <- (pI_before | pI_after)

    # Vertical (stacked) variant same two UMAPs but top/bottom, with a
    # single shared legend at the bottom.
    pI_after_v <- pI_after +
      theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.key.size = unit(2.5, "mm"),
            legend.text = element_text(size = 5),
            legend.title = element_text(size = 6)) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                  override.aes = list(size = 2, alpha = 1)))
    pI_v <- (pI_before / pI_after_v) + plot_layout(heights = c(1, 1.25))

  } else {
    cat("    Panel I \n")
    pI <- ggplot() + theme_void() + labs(title = "I: insufficient barcode match")
    pI_v <- pI
  }

} else {
  pI <- ggplot() + theme_void() + labs(title = "I: data not found")
  pI_v <- pI
}


# Panel J: Dot plot + +

cat("── Panel J: Dot plot ──\n")

dot_path <- "../fig1/dotplot_markers.csv"
if (file.exists(dot_path)) {
  dot <- fread(dot_path)
  
  # celltype marker
  ct_order_stair <- c("Epithelial", "Endothelial", "Fibroblast",
                      "Myeloid", "Mast", "T_NK", "B", "Plasma")
  
  marker_block_order <- c(
    "Epithelial", "Epithelial_prolif", "Supplementary",
    "Endothelial", "Fibroblast",
    "Myeloid", "Mast", "T_NK", "B", "Plasma"
  )
  dot$marker_group <- factor(dot$marker_group, levels = marker_block_order)
  dot <- dot[order(dot$marker_group)]
  gene_order <- unique(dot$gene)
  
  # x = gene y = celltype
  dot$gene <- factor(dot$gene, levels = gene_order)
  dot$celltype <- factor(dot$celltype, levels = ct_order_stair)
  
  pJ <- ggplot(dot, aes(x = gene, y = celltype)) +
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
      midpoint = 1, limits = c(0, 3), oob = scales::squish,
      guide = guide_colorbar(barwidth = unit(3.5, "mm"), barheight = unit(18, "mm"),
                             frame.colour = "black", frame.linewidth = 0.3,
                             ticks.colour = "black")
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = 5, face = "italic"),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(3, "mm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.grid = element_blank(),
      axis.line = element_blank()
    )
} else {
  cat("  ️ dotplot_markers.csv \n")
  pJ <- ggplot() + theme_void() + labs(title = "J: data not found")
}

cat("──  ──\n")

# panel

# ── panel (PDF + PNG) ──
save_panel <- function(plot, name, w, h) {
  ggsave(sprintf("FigS1_%s.pdf", name), plot, width = w, height = h,
         units = "mm", device = cairo_pdf)
  ggsave(sprintf("FigS1_%s.png", name), plot, width = w, height = h,
         units = "mm", dpi = 300, device = ragg::agg_png)
  cat(sprintf("  FigS1_%s pdf+png\n", name))
}

save_panel(pA, "A", 183, 45)
save_panel(pB, "B", 183, 45)
save_panel(pC, "C", 183, 55)
save_panel(pD, "D", 140, 50)
save_panel(pE, "E", 110, 90)
save_panel(pF, "F", 170, 75)
save_panel(pG, "G", 70,  65)
save_panel(pH, "H", 183, 60)
save_panel(pI, "I", 170, 75)
save_panel(pI_v, "I_vertical", 95, 165)
save_panel(pJ, "J", 170, 55)
cat("   panel PDF \n")

full <- (pA / pB / pC / pD /
           (pE | pF + plot_layout(widths = c(1, 2))) /
           (pG | pH + plot_layout(widths = c(1, 2.5))) /
           pI / pJ) +
  plot_layout(heights = c(0.8, 0.8, 0.9, 0.6, 1.3, 1, 1.2, 0.9))

ggsave("FigS1_full.pdf", full, width = 183, height = 380, units = "mm",
       device = cairo_pdf, limitsize = FALSE)
ggsave("FigS1_full.png", full, width = 183, height = 380, units = "mm",
       dpi = 300, device = ragg::agg_png, limitsize = FALSE)
cat("  FigS1_full saved\n")

cat("\n══════════════════════════════\n")
cat(" \n")
for (f in list.files(pattern = "^FigS1.*\\.(pdf|png)$")) {
  sz <- round(file.info(f)$size / 1e6, 1)
  cat(sprintf("  %s (%.1f MB)\n", f, sz))
}
cat("══════════════════════════════\n")