# Figure S9 — global vs subset-level CellChat heatmaps (Tumor / Metastasis)
# S9A: major-class weight   S9B: major-class count
# S9C: subset count         S9D: subset weight
suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

my_font <- "sans"
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  suppressPackageStartupMessages({ library(showtext); library(sysfonts) })
  sysfonts::font_add("Arial", regular = arial_p,
    bold = "~/.local/share/fonts/arialbd.ttf",
    italic = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300)
  my_font <- "Arial"
}

OUT_DIR <- "${PROJECT_ROOT}/results/fig6_panels"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cc_tumor <- readRDS("${PROJECT_ROOT}/data/processed/cellchat_rds/cellchat_Tumor.rds")
cc_met   <- readRDS("${PROJECT_ROOT}/data/processed/cellchat_rds/cellchat_Metastasis.rds")

# ---- major-class mapping ---------------------------------------------------
to_major <- function(x) {
  m <- character(length(x))
  for (i in seq_along(x)) {
    s <- as.character(x[i])
    m[i] <- if (startsWith(s, "Mal_"))            "Malignant"
       else if (startsWith(s, "Macro_"))          "Macrophage"
       else if (startsWith(s, "Neu_"))            "Neutrophil"
       else if (startsWith(s, "cDC") || s == "pDC") "DC"
       else if (startsWith(s, "Mono"))            "Monocyte"
       else if (s == "Endothelial")               "Endothelial"
       else if (s == "Fibroblast")                "Fibroblast"
       else if (s == "Epithelial_Normal")         "Epithelial"
       else if (s == "T_NK")                      "T_NK"
       else if (s == "B")                         "B"
       else if (s == "Plasma")                    "Plasma"
       else if (s == "Mast")                      "Mast"
       else                                       s
  }
  m
}

major_levels <- c("Malignant","Epithelial","Macrophage","Monocyte","DC","Mast",
                  "Neutrophil","T_NK","B","Plasma","Endothelial","Fibroblast")

# ---- aggregate fine matrix -> major matrix ---------------------------------
agg_mat <- function(M, fine, major_map, levels_use) {
  maj_send <- major_map[match(rownames(M), fine)]
  maj_recv <- major_map[match(colnames(M), fine)]
  out <- matrix(0, nrow = length(levels_use), ncol = length(levels_use),
                dimnames = list(levels_use, levels_use))
  for (i in seq_along(maj_send)) {
    for (j in seq_along(maj_recv)) {
      out[maj_send[i], maj_recv[j]] <- out[maj_send[i], maj_recv[j]] + M[i, j]
    }
  }
  out
}

# build a coarse-ident "shadow" CellChat object that netVisual_heatmap accepts
make_major_cc <- function(cc, levels_use) {
  fine <- levels(cc@idents)
  mapping <- to_major(fine)
  W <- agg_mat(cc@net$weight, fine, mapping, levels_use)
  N <- agg_mat(cc@net$count,  fine, mapping, levels_use)
  cc2 <- cc
  cc2@net$weight <- W
  cc2@net$count  <- N
  # rebuild @idents at major level so heatmap labels/colors match
  meta_major <- factor(mapping[match(as.character(cc@idents), fine)],
                       levels = levels_use)
  cc2@idents <- factor(meta_major, levels = levels_use)
  cc2@meta$labels <- meta_major
  cc2
}

# trim a major-level cc to only classes that appear in *both* T and M
restrict_levels <- function(cc, levels_use) {
  cc2 <- cc
  cc2@net$weight <- cc@net$weight[levels_use, levels_use]
  cc2@net$count  <- cc@net$count [levels_use, levels_use]
  cc2@idents     <- factor(as.character(cc@idents), levels = levels_use)
  cc2@meta$labels <- cc2@idents
  cc2
}

# ---- panel renderer (CellChat netVisual_heatmap, original colored style) ---
# Two heatmaps side-by-side share a single right-side legend slot, which
# gets squeezed/overlapped when both have marginal barplots. Bottom legend
# overflowed the canvas. Cleanest fix: write each condition to its own
# PDF/PNG file (suffix _tumor / _met). Each file has its own full legend
# area on the right. The user can compose them side-by-side in Illustrator.
draw_pair <- function(cc_l, cc_r, measure, title_l, title_r,
                       file_pdf, file_png, w_in, h_in, font_size = 9) {
  per_w <- w_in / 2 + 0.4   # each single panel a touch wider than half
  per_h <- h_in
  draw_one <- function(cc, title, pdf_path, png_path) {
    ht <- netVisual_heatmap(cc, measure = measure, title.name = title,
                            font.size = font_size, font.size.title = 11)
    pdf(pdf_path, width = per_w, height = per_h)
    draw(ht, heatmap_legend_side = "right",
         annotation_legend_side = "right",
         padding = unit(c(2, 4, 2, 6), "mm"))
    dev.off()
    png(png_path, width = per_w, height = per_h, units = "in", res = 300)
    draw(ht, heatmap_legend_side = "right",
         annotation_legend_side = "right",
         padding = unit(c(2, 4, 2, 6), "mm"))
    dev.off()
  }
  pdf_t <- sub("\\.pdf$", "_tumor.pdf", file_pdf)
  png_t <- sub("\\.png$", "_tumor.png", file_png)
  pdf_m <- sub("\\.pdf$", "_met.pdf",   file_pdf)
  png_m <- sub("\\.png$", "_met.png",   file_png)
  draw_one(cc_l, title_l, pdf_t, png_t)
  draw_one(cc_r, title_r, pdf_m, png_m)
  # remove combined-file leftovers from previous runs
  for (p in c(file_pdf, file_png)) if (file.exists(p)) file.remove(p)
}

# ---- new square-cell renderer (ComplexHeatmap, no CellChat dep) ------------
# Draws Tumor + Met side-by-side, square cellsize-mm cells, marginal
# barplots showing column sums (top, "Incoming") and row sums (right,
# "Outgoing"). Diagonal zeroed so self-loops don't dominate the colour scale.
draw_pair_square <- function(cc_l, cc_r, measure, title_l, title_r,
                              file_pdf, file_png,
                              cellsize_mm = 7, font_size = 8) {
  M_l <- cc_l@net[[measure]]
  M_r <- cc_r@net[[measure]]
  diag(M_l) <- 0; diag(M_r) <- 0
  vmax <- max(c(M_l, M_r), na.rm = TRUE)
  if (vmax == 0) vmax <- 1
  col_fun <- colorRamp2(c(0, vmax * 0.5, vmax),
                        c("#F7F7F7", "#FCAE91", "#CB181D"))

  out_l <- rowSums(M_l); inc_l <- colSums(M_l)
  out_r <- rowSums(M_r); inc_r <- colSums(M_r)
  bar_max <- max(c(out_l, inc_l, out_r, inc_r), na.rm = TRUE) * 1.05

  build_ht <- function(M, out_v, inc_v, name_id, show_legend = TRUE) {
    top_anno <- HeatmapAnnotation(
      `Incoming` = anno_barplot(inc_v, ylim = c(0, bar_max),
                  gp = gpar(fill = "grey55", col = NA),
                  bar_width = 0.85,
                  height = unit(13, "mm"),
                  axis_param = list(gp = gpar(fontsize = font_size - 2,
                                              fontfamily = my_font))),
      annotation_name_gp = gpar(fontsize = font_size - 1, fontfamily = my_font),
      annotation_name_side = "left",
      annotation_name_rot = 0
    )
    right_anno <- rowAnnotation(
      `Outgoing` = anno_barplot(out_v, ylim = c(0, bar_max),
                  gp = gpar(fill = "grey55", col = NA),
                  bar_width = 0.85,
                  width = unit(13, "mm"),
                  axis_param = list(gp = gpar(fontsize = font_size - 2,
                                              fontfamily = my_font))),
      annotation_name_gp = gpar(fontsize = font_size - 1, fontfamily = my_font),
      annotation_name_rot = -90
    )
    Heatmap(
      M, name = name_id, col = col_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      top_annotation = top_anno, right_annotation = right_anno,
      column_names_gp = gpar(fontsize = font_size, fontfamily = my_font),
      row_names_gp    = gpar(fontsize = font_size, fontfamily = my_font),
      column_names_rot = 45,
      rect_gp = gpar(col = "white", lwd = 0.4),
      border = TRUE, border_gp = gpar(col = "black", lwd = 0.6),
      width  = unit(ncol(M) * cellsize_mm, "mm"),
      height = unit(nrow(M) * cellsize_mm, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = font_size, fontfamily = my_font, fontface = "bold"),
        labels_gp = gpar(fontsize = font_size - 1, fontfamily = my_font),
        legend_height = unit(20, "mm"), grid_width = unit(2.5, "mm")
      ),
      column_title = NULL, row_title = NULL,
      show_heatmap_legend = show_legend
    )
  }

  ht_l <- build_ht(M_l, out_l, inc_l, name_id = measure, show_legend = FALSE)
  ht_r <- build_ht(M_r, out_r, inc_r, name_id = measure, show_legend = TRUE)

  n  <- nrow(M_l)
  cw <- (n * cellsize_mm * 2 + 100) / 25.4   # 100mm = labels + right barplot + legend gap
  ch <- (n * cellsize_mm + 65) / 25.4        # 65mm = top barplot + colnames rotation + title

  pdf(file_pdf, width = cw, height = ch)
  draw(ht_l + ht_r, ht_gap = unit(12, "mm"),
       column_title = paste0(title_l, "    |    ", title_r),
       column_title_gp = gpar(fontsize = font_size + 2, fontface = "bold",
                              fontfamily = my_font))
  dev.off()
  png(file_png, width = cw, height = ch, units = "in", res = 300)
  draw(ht_l + ht_r, ht_gap = unit(12, "mm"),
       column_title = paste0(title_l, "    |    ", title_r),
       column_title_gp = gpar(fontsize = font_size + 2, fontface = "bold",
                              fontfamily = my_font))
  dev.off()
}

# S9A / S9B — major-class heatmaps
cc_tumor_M <- make_major_cc(cc_tumor, major_levels)
cc_met_M   <- make_major_cc(cc_met,   major_levels)

# keep only classes present in BOTH conditions (avoid empty rows/cols)
present_T <- rownames(cc_tumor_M@net$weight)[rowSums(cc_tumor_M@net$weight) + colSums(cc_tumor_M@net$weight) > 0]
present_M <- rownames(cc_met_M@net$weight)  [rowSums(cc_met_M@net$weight)   + colSums(cc_met_M@net$weight)   > 0]
keep_major <- intersect(major_levels, intersect(present_T, present_M))
cat("S9 major-class kept:", paste(keep_major, collapse = ", "), "\n")

cc_tumor_M <- restrict_levels(cc_tumor_M, keep_major)
cc_met_M   <- restrict_levels(cc_met_M,   keep_major)

draw_pair(cc_tumor_M, cc_met_M, measure = "weight",
          title_l = "Tumor — interaction strength",
          title_r = "Metastasis — interaction strength",
          file_pdf = file.path(OUT_DIR, "figS9a_global_weight.pdf"),
          file_png = file.path(OUT_DIR, "figS9a_global_weight.png"),
          w_in = 16, h_in = 6, font_size = 10)
cat("S9A done.\n")

draw_pair(cc_tumor_M, cc_met_M, measure = "count",
          title_l = "Tumor — interaction number",
          title_r = "Metastasis — interaction number",
          file_pdf = file.path(OUT_DIR, "figS9b_global_count.pdf"),
          file_png = file.path(OUT_DIR, "figS9b_global_count.png"),
          w_in = 16, h_in = 6, font_size = 10)
cat("S9B done.\n")

# S9C / S9D — subset-level heatmaps (fine labels common to T and M)
# helper: hard-subset @net matrices and @idents to a chosen group set
subset_cc_groups <- function(cc, groups) {
  cc2 <- cc
  for (slot in c("weight","count","sum")) {
    if (!is.null(cc2@net[[slot]])) {
      M <- cc2@net[[slot]]
      keep <- intersect(groups, rownames(M))
      cc2@net[[slot]] <- M[keep, keep, drop = FALSE]
    }
  }
  # netVisual_heatmap only uses levels(@idents) for sidebar colors;
  # collapse @idents to a one-per-level factor to keep dims consistent
  cc2@idents <- factor(groups, levels = groups)
  cc2
}

common_fine <- intersect(levels(cc_tumor@idents), levels(cc_met@idents))
common_fine <- sort(common_fine)
cat("S9 subset-level kept (", length(common_fine), "):", paste(common_fine, collapse = ", "), "\n", sep = "")

cc_tumor_S <- subset_cc_groups(cc_tumor, common_fine)
cc_met_S   <- subset_cc_groups(cc_met,   common_fine)

draw_pair(cc_tumor_S, cc_met_S, measure = "count",
          title_l = "Tumor — interaction number",
          title_r = "Metastasis — interaction number",
          file_pdf = file.path(OUT_DIR, "figS9c_subset_count.pdf"),
          file_png = file.path(OUT_DIR, "figS9c_subset_count.png"),
          w_in = 20, h_in = 8, font_size = 7)
cat("S9C done.\n")

draw_pair(cc_tumor_S, cc_met_S, measure = "weight",
          title_l = "Tumor — interaction strength",
          title_r = "Metastasis — interaction strength",
          file_pdf = file.path(OUT_DIR, "figS9d_subset_weight.pdf"),
          file_png = file.path(OUT_DIR, "figS9d_subset_weight.png"),
          w_in = 20, h_in = 8, font_size = 7)
cat("S9D done.\n")

cat("\nAll S9 panels written to", OUT_DIR, "\n")
