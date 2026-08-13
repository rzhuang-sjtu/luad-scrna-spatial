# Supplementary Fig S7E (Hallmark heatmap) + S7F (GO BP bar plots)
# Data: ~/luad/results/fig5_plot_data/figS7e_*.csv  &  figS7f_*.csv
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork); library(scales)
  library(ComplexHeatmap); library(circlize); library(grid)
})

DATA_DIR  <- "${PROJECT_ROOT}/results/fig5_plot_data"
OUT_DIR   <- "${PROJECT_ROOT}/results/fig5_panels"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SUBTYPE_ORDER <- c(
  "Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low","Neu_IFN_response",
  "Neu_Angiogenic","Neu_ECM_remodeling","Neu_Metastatic"
)
SUBTYPE_COLOR <- c(
  Neu_Inflammatory   = "#D62728",
  Neu_OSM_priming    = "#FF7F0E",
  Neu_OSM_low        = "#FDBF6F",
  Neu_IFN_response   = "#9467BD",
  Neu_Angiogenic     = "#1F77B4",
  Neu_ECM_remodeling = "#2CA02C",
  Neu_Metastatic     = "#8C564B"
)

# S7E — Hallmark z-score heatmap
z <- fread(file.path(DATA_DIR, "figS7e_hallmark_zscore.csv"))
mat <- as.matrix(z[, -1, with = FALSE])
rownames(mat) <- z$subtype
mat <- mat[SUBTYPE_ORDER, , drop = FALSE]      # subtype rows, hallmark cols
# Column names in this CSV are already pretty (mixed case, e.g. "Hypoxia")

HL <- c("Epithelial Mesenchymal Transition",
        "Inflammatory Response",
        "TNF-alpha Signaling via NF-kB",
        "IL-6/JAK/STAT3 Signaling",
        "Reactive Oxygen Species Pathway",
        "Hypoxia")
HL <- intersect(HL, colnames(mat))

# Prune to top-25 hallmarks: union of (highlight set) and (top-by-range)
range_score <- apply(mat, 2, function(x) max(x) - min(x))
keep_top    <- names(sort(range_score, decreasing = TRUE))
keep_set    <- unique(c(HL, head(keep_top, 25)))[1:25]
mat <- mat[, keep_set, drop = FALSE]

# column order: by which subtype each hallmark scores highest
ord_cols <- order(apply(mat, 2, which.max), -apply(mat, 2, max))
mat <- mat[, ord_cols, drop = FALSE]

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

# Subtype identity is already conveyed by the row names; skip the left
# coloured strip to reclaim horizontal space and fit A4 portrait width.

col_face <- ifelse(colnames(mat) %in% HL, "bold", "plain")
col_col  <- ifelse(colnames(mat) %in% HL, "#B2182B", "black")

# Wrap only the genuinely long hallmark names (> 28 chars) at word
# boundaries; short labels like "IL-6/JAK/STAT3 Signaling" are kept on
# a single line so the rotated text reads cleanly.
colnames(mat) <- vapply(colnames(mat), function(s) {
  if (nchar(s) <= 28) return(s)
  if (requireNamespace("stringr", quietly = TRUE))
    stringr::str_wrap(s, width = 28)
  else s
}, character(1))

ht <- Heatmap(
  mat,
  name = "z-score",
  col  = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp  = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 6, fontface = col_face, col = col_col),
  column_names_rot = 45,
  column_names_max_height = unit(40, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.25),
  # 25 cols x 5 mm = 125 mm body, 7 rows x 4.5 mm = 31.5 mm body
  width  = unit(ncol(mat) * 5,   "mm"),
  height = unit(nrow(mat) * 4.5, "mm"),
  heatmap_legend_param = list(title = "z-score", at = c(-2, 0, 2),
                              legend_height = unit(18, "mm"),
                              legend_width  = unit(3, "mm"))
)

# Canvas tuned to fit A4 portrait usable width (~190 mm), still with
# generous padding so rotated column names + row names never clip.
pdf(file.path(OUT_DIR, "figS7e_hallmark_heatmap.pdf"),
    width = 7.5, height = 4.5)
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(12, 6, 6, 8), "mm"))
dev.off()

png(file.path(OUT_DIR, "figS7e_hallmark_heatmap.png"),
    width = 7.5, height = 4.5, units = "in", res = 300)
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(12, 6, 6, 8), "mm"))
dev.off()
cat("S7E saved (A4-fit, padded edges)\n")

# S7F — top-10 GO BP bar plots, one panel per subtype
go <- fread(file.path(DATA_DIR, "figS7f_go_bp_enrichr.csv"))
# pretty term: drop "(GO:0001234)" suffix, hard-wrap long terms
go[, term_clean := str_squish(sub("\\(GO:[0-9]+\\)$", "", term))]

# top-5 per subtype (was 10) so the 7-panel grid fits the A4-portrait bottom half
top5 <- go[order(p_adj), head(.SD, 5), by = subtype]
top5[, neg_log10_p := -log10(pmax(p_adj, 1e-30))]

# Lock the GeneRatio fill scale so the single shared legend is one-to-one
gr_lim <- range(top5$gene_ratio, na.rm = TRUE)

build_panel <- function(sub) {
  d <- top5[subtype == sub]
  if (nrow(d) == 0) return(NULL)
  d <- d[order(neg_log10_p)]
  d[, term_clean := factor(term_clean, levels = unique(d$term_clean))]
  # Adaptive y-axis font: shrink only when terms get long, never below 4.2pt.
  # Bar geometry / row spacing unchanged.
  max_len <- max(nchar(as.character(d$term_clean)))
  y_size  <- max(4.2, min(6.5, 6.5 * 35 / max(max_len, 1)))
  fill_col <- unname(SUBTYPE_COLOR[sub])
  ggplot(d, aes(x = neg_log10_p, y = term_clean, fill = gene_ratio)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(low = scales::alpha(fill_col, 0.30),
                        high = fill_col,
                        name = "GeneRatio", limits = gr_lim) +
    labs(title = sub, x = expression(-log[10]~italic(p[adj])), y = NULL) +
    theme_classic(base_size = 7.5) +
    theme(plot.title  = element_text(face = "bold", size = 8,
                                     color = unname(SUBTYPE_COLOR[sub])),
          axis.text.y = element_text(size = y_size, color = "black"),
          axis.text.x = element_text(size = 6.5, color = "black"),
          axis.title.x = element_text(size = 7),
          legend.position = "right",
          legend.key.size = unit(2.5, "mm"),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 5.5),
          plot.margin = margin(2, 4, 2, 2, "pt"))
}

panels <- lapply(SUBTYPE_ORDER, build_panel)
panels <- panels[!vapply(panels, is.null, logical(1))]
panels[[length(panels) + 1]] <- patchwork::plot_spacer()  # 8th slot for 2x4

# Each panel keeps its own subtype-coloured GeneRatio legend on the right.
combined <- wrap_plots(panels, ncol = 2, nrow = 4, byrow = TRUE)

pdf(file.path(OUT_DIR, "figS7f_go_enrichment.pdf"),
    width = 7.1, height = 6.0)
print(combined)
dev.off()
png(file.path(OUT_DIR, "figS7f_go_enrichment.png"),
    width = 7.1, height = 6.0, units = "in", res = 300)
print(combined)
dev.off()
cat("S7F (7 subtypes x top-5 GO, 2x4 grid, A4 bottom half) saved\n")
cat("\nAll S7 panels written to ", OUT_DIR, "\n", sep = "")
