## Figure 7 v3 — Spatial Multi-Omics Analysis of LUAD
## Critical fixes vs v2:
##   • H&E pathology image as background for every spatial panel
##   • Strong-contrast white→signature-color ramps (no faint viridis)
##   • Explicit cohort + sample annotations (no overlapping titles)
##   • Square heatmap cells with borders + visible row/column names
##   • ROI panels: H&E underlay + ROI-only colored spots, zoomed
##   • COMMOT vector field: visible arrow length + bright signal color

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(scales); library(dplyr); library(tidyr)
  library(viridis); library(grid); library(png); library(jsonlite)
  library(ggnewscale); library(showtext); library(sysfonts)
})

R_DATA   <- "${DATA_ROOT}/ST/results/r_data"
PER      <- file.path(R_DATA, "per_section")
HE_DIR   <- file.path(R_DATA, "he")
FIG_DIR  <- "${WORK_ROOT}/luad_figures/fig7"
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

## Fonts
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p,
    bold = "~/.local/share/fonts/arialbd.ttf",
    italic = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300)
  my_font <- "Arial"
} else { my_font <- "sans" }

SEC_A <- list(stem = "EMTAB13530__P16_T1", label = "P16_T1 · E-MTAB-13530 · LUAD")
SEC_B <- list(stem = "EMTAB13530__P10_T1", label = "P10_T1 · E-MTAB-13530 · LUAD")
SEC_R <- list(stem = "EMTAB13530__P16_T2", label = "P16_T2 · E-MTAB-13530 · LUAD")  # ROI panels

neu_colors <- c(
  Neu_Inflammatory   = "#E64B35", Neu_Angiogenic   = "#F39B7F",
  Neu_Metastatic     = "#3C5488", Neu_ECM_remodeling = "#4DBBD5",
  Neu_OSM_priming    = "#00A087", Neu_OSM_low      = "#8491B4",
  Neu_IFN_response   = "#91D1C2")
mp_colors <- c(MP1 = "#E64B35", MP2 = "#4DBBD5", MP3 = "#00A087",
               MP4 = "#3C5488", MP5 = "#999999")
macro_colors <- c(Macro_C1QC = "#4DBBD5", Macro_FCN1 = "#E64B35",
  Macro_FOLR2 = "#00A087", Macro_MARCO = "#3C5488",
  Macro_SPP1  = "#F39B7F", Macro_general = "#8491B4",
  Macro_prolif = "#91D1C2")
ct_colors <- c(Fibroblast = "#B09C85", Endothelial = "#651FFF",
  T_NK = "#FFCB5C", B = "#9C27B0", Plasma = "#D9D9D9", Mast = "#F0A23B",
  Epithelial = "#D2B48C", Mono_nonclassical = "#999999",
  cDC1 = "#80CBC4", cDC2 = "#7570B3", cDC_LAMP3 = "#984EA3", pDC = "#B2B2B2",
  Malignant = "#C73E2A")
all_colors <- c(neu_colors, macro_colors, ct_colors)

## map cell-type → high contrast gradient color
ramp_high <- function(ct) {
  if (ct %in% names(all_colors)) all_colors[[ct]] else "#B2182B"
}

theme_panel <- function(base_size = 8) {
  theme_void(base_family = my_font, base_size = base_size) +
    theme(
      plot.subtitle = element_text(face = "bold", size = rel(1), color = "black",
                                    margin = margin(b = 1.5)),
      plot.caption  = element_text(face = "plain", size = rel(0.75),
                                   color = "grey25", hjust = 0.5,
                                   margin = margin(t = 1)),
      legend.title = element_text(size = rel(0.8)),
      legend.text  = element_text(size = rel(0.75)),
      legend.key.width  = unit(0.18, "cm"),
      legend.key.height = unit(0.42, "cm"),
      legend.margin = margin(0,0,0,1),
      plot.margin = margin(2, 2, 2, 2)
    )
}

load_section <- function(stem) fread(file.path(PER, paste0(stem, ".csv")))

load_he <- function(stem) {
  he_dir <- file.path(HE_DIR, stem)
  img_path <- file.path(he_dir, "tissue_hires_image.png")
  sf_path  <- file.path(he_dir, "scalefactors_json.json")
  img <- png::readPNG(img_path)
  sf  <- jsonlite::fromJSON(sf_path)
  list(img = img, sf = sf$tissue_hires_scalef,
       w = dim(img)[2], h = dim(img)[1])
}

prep_xy <- function(d, he) {
  d$x_he <-  d$spatial1 * he$sf
  d$y_he <- -d$spatial2 * he$sf      # flip for ggplot Y-up
  d
}

DA <- load_section(SEC_A$stem); HEA <- load_he(SEC_A$stem); DA <- prep_xy(DA, HEA)
DB <- load_section(SEC_B$stem); HEB <- load_he(SEC_B$stem); DB <- prep_xy(DB, HEB)
DR <- load_section(SEC_R$stem); HER <- load_he(SEC_R$stem); DR <- prep_xy(DR, HER)
cat(sprintf("[A] %s  spots=%d  he=%dx%d  sf=%.4f\n", SEC_A$stem, nrow(DA), HEA$w, HEA$h, HEA$sf))
cat(sprintf("[B] %s  spots=%d  he=%dx%d  sf=%.4f\n", SEC_B$stem, nrow(DB), HEB$w, HEB$h, HEB$sf))
cat(sprintf("[R] %s  spots=%d  he=%dx%d  sf=%.4f  ROI=%d\n",
            SEC_R$stem, nrow(DR), HER$w, HER$h, HER$sf, sum(DR$roi == 1)))

## Spatial heatmap on H&E underlay.
## - low values → transparent (H&E shows through)
## - high values → bright signature color
plot_spatial_he <- function(d, color_col, he, panel_name,
                             high_color = "#B2182B",
                             low_color  = "white",
                             value_clip = c(0.02, 0.98),
                             point_size = 0.85,
                             alpha_max  = 0.92,
                             diverging  = FALSE) {
  v <- d[[color_col]]
  if (diverging) {
    a <- max(abs(quantile(v, value_clip[1], na.rm = TRUE)),
             abs(quantile(v, value_clip[2], na.rm = TRUE)))
    lim <- c(-a, a)
  } else {
    lim <- c(quantile(v, value_clip[1], na.rm = TRUE),
             quantile(v, value_clip[2], na.rm = TRUE))
  }
  ## alpha by quantile rank within positive direction
  if (diverging) {
    a_norm <- pmin(1, abs(v) / lim[2])
  } else {
    rng <- diff(lim)
    a_norm <- pmin(1, pmax(0, (v - lim[1]) / ifelse(rng>0, rng, 1)))
  }
  d$alpha_pt <- alpha_max * a_norm
  p <- ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d, aes(x = x_he, y = y_he, color = .data[[color_col]],
                              alpha = alpha_pt),
               size = point_size, stroke = 0)
  if (diverging) {
    p <- p + scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                    midpoint = 0, limits = lim,
                                    oob = scales::squish, name = "")
  } else {
    p <- p + scale_color_gradient(low = low_color, high = high_color,
                                   limits = lim, oob = scales::squish, name = "")
  }
  p + scale_alpha_identity() +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_panel() +
    labs(subtitle = panel_name)
}

make_grid <- function(d, he, sample_label, panels_list) {
  ps <- lapply(panels_list, function(spec) {
    plot_spatial_he(d, spec$col, he, spec$name,
                    high_color = spec$color,
                    point_size = 0.85, alpha_max = 0.95)
  })
  wrap_plots(ps, nrow = 2, ncol = 3) +
    plot_annotation(
      title = sample_label,
      theme = theme(plot.title = element_text(face = "bold", size = 9,
                                                family = my_font, hjust = 0.02))
    )
}

cell_panels_main <- list(
  list(col = "ct_Neu_Inflammatory", name = "Neu_Inflammatory", color = unname(neu_colors["Neu_Inflammatory"])),
  list(col = "ct_Neu_OSM_priming",  name = "Neu_OSM_priming",  color = unname(neu_colors["Neu_OSM_priming"])),
  list(col = "ct_Macro_SPP1",       name = "Macro_SPP1",       color = unname(macro_colors["Macro_SPP1"])),
  list(col = "ct_Fibroblast",       name = "Fibroblast",       color = "#5C4830"),
  list(col = "ct_Endothelial",      name = "Endothelial",      color = unname(ct_colors["Endothelial"])),
  list(col = "ct_Malignant",        name = "Malignant",        color = unname(ct_colors["Malignant"]))
)

P_7A <- make_grid(DA, HEA, SEC_A$label, cell_panels_main)
ggsave(file.path(PANEL_DIR, "7A.pdf"), P_7A, width = 130, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7A.png"), P_7A, width = 130, height = 90, units = "mm", dpi = 300)

P_7B <- make_grid(DB, HEB, SEC_B$label, cell_panels_main)
ggsave(file.path(PANEL_DIR, "7B.pdf"), P_7B, width = 130, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7B.png"), P_7B, width = 130, height = 90, units = "mm", dpi = 300)

plot_composite_he <- function(d, he, sample_label, focus_cts, palette,
                                top_q = 0.5, point_size = 1.1, alpha_pt = 0.9) {
  ab_cols <- paste0("ct_", focus_cts)
  ab_cols <- ab_cols[ab_cols %in% names(d)]
  m <- as.matrix(d[, ..ab_cols])
  colnames(m) <- sub("^ct_", "", ab_cols)
  best_idx <- apply(m, 1, which.max)
  d2 <- copy(d)
  d2$dom_ct  <- factor(colnames(m)[best_idx], levels = colnames(m))
  d2$dom_val <- m[cbind(seq_len(nrow(m)), best_idx)]
  high_mask <- rep(FALSE, nrow(d2))
  for (ct in colnames(m)) {
    idx <- d2$dom_ct == ct
    if (sum(idx) < 5) next
    th <- quantile(d2$dom_val[idx], top_q, na.rm = TRUE)
    high_mask <- high_mask | (idx & d2$dom_val >= th)
  }
  d2$is_high <- high_mask
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d2[is_high == TRUE],
               aes(x = x_he, y = y_he, color = dom_ct),
               size = point_size, stroke = 0, alpha = alpha_pt) +
    scale_color_manual(values = palette, name = NULL,
                       guide = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_panel() +
    theme(legend.position = "right",
          legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 7)) +
    labs(subtitle = sprintf("Composite (top-50%% per cell type)"),
         caption  = sample_label)
}

focus_cts <- c("Neu_Inflammatory","Neu_OSM_priming","Macro_SPP1",
               "Fibroblast","Endothelial","Malignant")
focus_palette <- c(
  Neu_Inflammatory = unname(neu_colors["Neu_Inflammatory"]),
  Neu_OSM_priming  = unname(neu_colors["Neu_OSM_priming"]),
  Macro_SPP1       = unname(macro_colors["Macro_SPP1"]),
  Fibroblast       = "#5C4830",
  Endothelial      = unname(ct_colors["Endothelial"]),
  Malignant        = unname(ct_colors["Malignant"]))

P_7C <- plot_composite_he(DA, HEA, SEC_A$label, focus_cts, focus_palette)
ggsave(file.path(PANEL_DIR, "7C.pdf"), P_7C, width = 100, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7C.png"), P_7C, width = 100, height = 95, units = "mm", dpi = 300)
P_7D <- plot_composite_he(DB, HEB, SEC_B$label, focus_cts, focus_palette)
ggsave(file.path(PANEL_DIR, "7D.pdf"), P_7D, width = 100, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7D.png"), P_7D, width = 100, height = 95, units = "mm", dpi = 300)

plot_commot_field_he <- function(d, he, sample_label, pathway = "OSM",
                                  arrow_topN = 60, target_arrow_frac = 0.05) {
  send_col <- sprintf("commot_s_%s", pathway)
  tot_col  <- sprintf("commot_total_%s", pathway)
  vfx <- sprintf("vf_s_%s_dx", pathway); vfy <- sprintf("vf_s_%s_dy", pathway)
  d <- d[!is.na(get(vfx)) & !is.na(get(vfy))]
  d$tot <- d[[tot_col]]
  d$tot[is.na(d$tot)] <- 0
  ## auto scale arrow length to fraction of section span
  span <- max(diff(range(d$x_he)), diff(range(d$y_he)))
  vmag <- sqrt(d[[vfx]]^2 + d[[vfy]]^2)
  med  <- median(vmag[vmag > 0], na.rm = TRUE)
  scale_factor <- if (is.na(med) || med == 0) 0 else span * target_arrow_frac / med
  ord <- order(d[[send_col]], decreasing = TRUE)
  arr <- d[ord[seq_len(min(arrow_topN, nrow(d)))]]
  arr$ax <- arr$x_he; arr$ay <- arr$y_he
  ## scale_factor is computed in he-space; vector field is in commot-output units.
  ## After scale_factor * vfx, the arrow is already span * target_arrow_frac long.
  arr$bx <- arr$x_he + scale_factor * arr[[vfx]]
  arr$by <- arr$y_he - scale_factor * arr[[vfy]]
  ## use signal magnitude as alpha so dim spots fade into H&E
  d$alpha_pt <- pmin(1, d$tot / quantile(d$tot[d$tot>0], 0.95, na.rm = TRUE)) * 0.95
  d$alpha_pt[is.na(d$alpha_pt)] <- 0
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d, aes(x = x_he, y = y_he, color = tot, alpha = alpha_pt),
               size = 0.8, stroke = 0) +
    scale_color_gradient(low = "#FFFAF0", high = "#B2182B",
                          limits = c(0, quantile(d$tot[d$tot>0], 0.95, na.rm=TRUE)),
                          oob = scales::squish, name = sprintf("%s total", pathway)) +
    scale_alpha_identity() +
    geom_segment(data = arr, aes(x = ax, y = ay, xend = bx, yend = by),
                 arrow = arrow(length = unit(0.10, "cm"), type = "closed"),
                 color = "#003C30", linewidth = 0.32) +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_panel() +
    labs(subtitle = sprintf("COMMOT %s sender vectors", pathway),
         caption = sample_label)
}

P_7E <- plot_commot_field_he(DA, HEA, SEC_A$label, "OSM",
                              arrow_topN = 70, target_arrow_frac = 0.10)
ggsave(file.path(PANEL_DIR, "7E.pdf"), P_7E, width = 100, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7E.png"), P_7E, width = 100, height = 95, units = "mm", dpi = 300)

P_7F <- plot_spatial_he(DA, "gex_OSM",  HEA, "OSM expression",
                         high_color = "#00BCD4", point_size = 1.0)
P_7F <- P_7F + labs(caption = SEC_A$label)
ggsave(file.path(PANEL_DIR, "7F.pdf"), P_7F, width = 80, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7F.png"), P_7F, width = 80, height = 80, units = "mm", dpi = 300)

P_7G <- plot_spatial_he(DA, "gex_IL1B", HEA, "IL1B expression",
                         high_color = "#FFB300", point_size = 1.0)
P_7G <- P_7G + labs(caption = SEC_A$label)
ggsave(file.path(PANEL_DIR, "7G.pdf"), P_7G, width = 80, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7G.png"), P_7G, width = 80, height = 80, units = "mm", dpi = 300)

imp <- fread(file.path(R_DATA, "misty_aggregated_importance.csv"))
imp_emtab <- imp[cohort == "EMTAB13530" & view == "intra"]
imp_emtab[, target := sub("^progeny_", "", Target)]
imp_emtab[, target := gsub("\\.", "-", target)]
imp_emtab <- imp_emtab[!grepl("^progeny_", Predictor)]
target_keep <- c("NFkB","JAK-STAT","TNFa","TGFb","Hypoxia","MAPK","EGFR","p53","VEGF")
imp_emtab <- imp_emtab[target %in% target_keep]
mat <- dcast(imp_emtab, Predictor ~ target, value.var = "mean_importance", fill = 0)
mat_m <- as.matrix(mat[, -1, with = FALSE])
rownames(mat_m) <- mat$Predictor
mat_m <- mat_m[order(-rowMeans(abs(mat_m))), ]

vrange <- max(abs(mat_m), na.rm = TRUE)
col_fun <- colorRamp2(c(-vrange, 0, vrange),
                      c("#2166AC", "white", "#67001F"))

ht <- Heatmap(
  mat_m, name = "MISTy\nimportance", col = col_fun,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_names_side = "left", column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_rot = 45,
  ## square cells with white borders
  rect_gp = gpar(col = "white", lwd = 1.1),
  width  = unit(0.8, "cm") * ncol(mat_m),
  height = unit(0.45, "cm") * nrow(mat_m),
  row_title = "Cell type (Predictor)",
  column_title = "PROGENy pathway (Target)",
  row_title_gp = gpar(fontsize = 8, fontfamily = my_font),
  column_title_gp = gpar(fontsize = 8, fontfamily = my_font),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = my_font),
    labels_gp = gpar(fontsize = 7, fontfamily = my_font),
    grid_height = unit(0.32, "cm"), grid_width = unit(0.32, "cm"),
    title_position = "topleft"),
  border = TRUE
)

pdf(file.path(PANEL_DIR, "7H.pdf"), width = 5.6, height = 7.0)
draw(ht, padding = unit(c(2, 6, 2, 4), "mm"))
dev.off()
png(file.path(PANEL_DIR, "7H.png"), width = 5.6, height = 7.0, units = "in", res = 300)
draw(ht, padding = unit(c(2, 6, 2, 4), "mm"))
dev.off()

P_7I <- plot_spatial_he(DR, "progeny_NFkB", HER, "PROGENy NFkB",
                         diverging = TRUE, point_size = 0.9, alpha_max = 0.92)
P_7I <- P_7I + labs(caption = SEC_R$label)
ggsave(file.path(PANEL_DIR, "7I.pdf"), P_7I, width = 80, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7I.png"), P_7I, width = 80, height = 80, units = "mm", dpi = 300)

jak_col <- if ("progeny_JAK-STAT" %in% names(DR)) "progeny_JAK-STAT" else "progeny_JAK.STAT"
P_7J <- plot_spatial_he(DR, jak_col, HER, "PROGENy JAK-STAT",
                         diverging = TRUE, point_size = 0.9, alpha_max = 0.92)
P_7J <- P_7J + labs(caption = SEC_R$label)
ggsave(file.path(PANEL_DIR, "7J.pdf"), P_7J, width = 80, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7J.png"), P_7J, width = 80, height = 80, units = "mm", dpi = 300)

plot_roi_outline <- function(d, he, sample_label) {
  d2 <- copy(d)
  d2$roi_lab <- factor(ifelse(d2$roi == 1, "ROI", "non-ROI"),
                       levels = c("non-ROI", "ROI"))
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d2[d2$roi == 0],
               aes(x = x_he, y = y_he), color = "grey85",
               size = 0.5, stroke = 0, alpha = 0.4) +
    geom_point(data = d2[d2$roi == 1],
               aes(x = x_he, y = y_he), color = "#B2182B",
               size = 0.95, stroke = 0, alpha = 0.95) +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_panel() +
    labs(subtitle = "ROI: NFkB-high & Neutrophil-high",
         caption  = sample_label)
}
P_7K <- plot_roi_outline(DR, HER, SEC_R$label)
ggsave(file.path(PANEL_DIR, "7K.pdf"), P_7K, width = 90, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7K.png"), P_7K, width = 90, height = 80, units = "mm", dpi = 300)

plot_roi_zoom <- function(d, he, sample_label, color_col, panel_name,
                           palette = "rocket", high_color = "#B2182B",
                           diverging = FALSE,
                           point_size_roi = 1.6, point_size_bg = 0.45,
                           margin_frac = 0.10) {
  d2 <- copy(d)
  roi_d <- d2[d2$roi == 1]
  if (nrow(roi_d) == 0) return(ggplot() + theme_panel() + labs(subtitle = panel_name))
  xr <- range(roi_d$x_he); yr <- range(roi_d$y_he)
  pad <- max(diff(xr), diff(yr)) * margin_frac
  xr2 <- xr + c(-pad, pad); yr2 <- yr + c(-pad, pad)
  ## clip raster to the zoom window so H&E only shows in zoomed region
  v <- roi_d[[color_col]]
  if (diverging) {
    a <- max(abs(quantile(v, 0.05, na.rm = TRUE)), abs(quantile(v, 0.95, na.rm = TRUE)))
    lim <- c(-a, a)
  } else {
    lim <- c(quantile(v, 0.02, na.rm = TRUE), quantile(v, 0.98, na.rm = TRUE))
  }
  p <- ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d2[d2$roi == 0],
               aes(x = x_he, y = y_he),
               color = "grey85", size = point_size_bg, stroke = 0, alpha = 0.4) +
    geom_point(data = roi_d, aes(x = x_he, y = y_he,
                                  color = .data[[color_col]]),
               size = point_size_roi, stroke = 0, alpha = 0.95)
  if (diverging) {
    p <- p + scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                    midpoint = 0, limits = lim, oob = scales::squish, name = "")
  } else {
    p <- p + scale_color_gradient(low = "white", high = high_color,
                                   limits = lim, oob = scales::squish, name = "")
  }
  p + coord_fixed(xlim = xr2, ylim = yr2, expand = FALSE) +
    theme_panel() +
    labs(subtitle = panel_name, caption = sample_label)
}

## 7L: MP1-4 (2x2) — high-contrast palette against magenta H&E (gold/cyan/purple/green)
P_7L <- (plot_roi_zoom(DR, HER, SEC_R$label, "MP1_score", "MP1 (Stress/AP-1)",  high_color = "#FFB300") |
         plot_roi_zoom(DR, HER, SEC_R$label, "MP2_score", "MP2 (Proliferative)", high_color = "#00BCD4")) /
        (plot_roi_zoom(DR, HER, SEC_R$label, "MP3_score", "MP3 (EMT/IFN)",       high_color = "#651FFF") |
         plot_roi_zoom(DR, HER, SEC_R$label, "MP4_score", "MP4 (AT2-like)",      high_color = "#43A047"))
ggsave(file.path(PANEL_DIR, "7L.pdf"), P_7L, width = 130, height = 110, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7L.png"), P_7L, width = 130, height = 110, units = "mm", dpi = 300)

## 7M: TF expression (2x2) — 4 distinct hues to avoid the prior all-magenta blend
P_7M <- (plot_roi_zoom(DR, HER, SEC_R$label, "gex_ATF3",   "ATF3",   high_color = "#FFB300") |
         plot_roi_zoom(DR, HER, SEC_R$label, "gex_FOSB",   "FOSB",   high_color = "#00BCD4")) /
        (plot_roi_zoom(DR, HER, SEC_R$label, "gex_JUN",    "JUN",    high_color = "#651FFF") |
         plot_roi_zoom(DR, HER, SEC_R$label, "gex_NFKBIA", "NFKBIA", high_color = "#43A047"))
ggsave(file.path(PANEL_DIR, "7M.pdf"), P_7M, width = 130, height = 110, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7M.png"), P_7M, width = 130, height = 110, units = "mm", dpi = 300)

## 7N: Neu + Macro + Endo (2x2)
P_7N <- (plot_roi_zoom(DR, HER, SEC_R$label, "ct_Neu_Inflammatory",
                       "Neu_Inflammatory", high_color = unname(neu_colors["Neu_Inflammatory"])) |
         plot_roi_zoom(DR, HER, SEC_R$label, "ct_Neu_OSM_priming",
                       "Neu_OSM_priming",  high_color = unname(neu_colors["Neu_OSM_priming"]))) /
        (plot_roi_zoom(DR, HER, SEC_R$label, "ct_Macro_SPP1",
                       "Macro_SPP1", high_color = unname(macro_colors["Macro_SPP1"])) |
         plot_roi_zoom(DR, HER, SEC_R$label, "ct_Endothelial",
                       "Endothelial", high_color = unname(ct_colors["Endothelial"])))
ggsave(file.path(PANEL_DIR, "7N.pdf"), P_7N, width = 130, height = 110, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7N.png"), P_7N, width = 130, height = 110, units = "mm", dpi = 300)

## 7O: NFkB + IL1B (1x2)
P_7O <- plot_roi_zoom(DR, HER, SEC_R$label, "progeny_NFkB", "NFkB activity", diverging = TRUE) |
        plot_roi_zoom(DR, HER, SEC_R$label, "gex_IL1B",     "IL1B expression", high_color = "#FFB300")
ggsave(file.path(PANEL_DIR, "7O.pdf"), P_7O, width = 130, height = 65, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7O.png"), P_7O, width = 130, height = 65, units = "mm", dpi = 300)

## 7P: JAK-STAT + OSM (1x2)
P_7P <- plot_roi_zoom(DR, HER, SEC_R$label, jak_col, "JAK-STAT activity", diverging = TRUE) |
        plot_roi_zoom(DR, HER, SEC_R$label, "gex_OSM", "OSM expression", high_color = "#00BCD4")
ggsave(file.path(PANEL_DIR, "7P.pdf"), P_7P, width = 130, height = 65, units = "mm", device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "7P.png"), P_7P, width = 130, height = 65, units = "mm", dpi = 300)


## per-spot Mann-Whitney + BH-FDR via step09e_roi_significance.py
##   x = cohort,  y = metric,  size = |delta|,  fill = signed delta
##   stars overlaid in white at the bubble center
roi_pv_path <- file.path(R_DATA, "roi_vs_nonroi_aggregate_pvalues.csv")
if (file.exists(roi_pv_path)) {
  roi_q <- fread(roi_pv_path)
  roi_q[is.na(sig) | sig == "", sig := ""]
  in_both <- roi_q[, .(n_coh = uniqueN(cohort)), by = metric][n_coh == 2, metric]
  roi_q <- roi_q[metric %in% in_both]
  rk <- roi_q[, .(absd = mean(abs(delta), na.rm = TRUE)), by = metric]
  top12 <- rk[order(-absd)][1:12, metric]
  roi_q <- roi_q[metric %in% top12]
  roi_q[, metric := factor(metric, levels = rev(top12))]
  roi_q[, cohort := factor(cohort, levels = c("EMTAB13530","Okamura"))]
  lim_d <- max(abs(roi_q$delta), na.rm = TRUE)

  P_7Q <- ggplot(roi_q, aes(x = cohort, y = metric)) +
    geom_point(aes(size = abs(delta), fill = delta),
               shape = 21, color = "grey25", stroke = 0.25) +
    geom_text(aes(label = sig), vjust = 0.55, hjust = 0.5,
              size = 1.7, family = my_font, fontface = "bold", color = "white") +
    scale_x_discrete(labels = c(EMTAB13530 = "E-MTAB-13530",
                                Okamura    = "Takano 2024"),
                     position = "top",
                     expand = expansion(add = c(0.7, 0.7))) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                         midpoint = 0, limits = c(-lim_d, lim_d),
                         oob = scales::squish, name = "Delta",
                         guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                barheight = unit(20, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.25)) +
    scale_size_continuous(range = c(1.6, 7), name = "|Delta|",
                          guide = guide_legend(override.aes = list(fill = "grey55"))) +
    labs(subtitle = "ROI vs non-ROI · top-12 |delta| · per-spot Mann-Whitney + BH-FDR\n*** <0.001  ** <0.01  * <0.05",
         x = NULL, y = NULL) +
    theme_minimal(base_family = my_font, base_size = 7) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
          panel.background   = element_rect(fill = "white", color = NA),
          axis.text   = element_text(color = "black"),
          axis.text.x = element_text(size = 7, face = "bold", margin = margin(t = 1)),
          axis.ticks  = element_blank(),
          axis.line   = element_blank(),
          plot.subtitle = element_text(size = 6, lineheight = 1.1, margin = margin(b = 2)),
          legend.position  = "right",
          legend.box       = "vertical",
          legend.key.size  = unit(3, "mm"),
          legend.text      = element_text(size = 6),
          legend.title     = element_text(size = 6, face = "bold"),
          plot.margin = margin(2, 4, 2, 4))
  ggsave(file.path(PANEL_DIR, "7Q.pdf"), P_7Q, width = 95, height = 95, units = "mm", device = cairo_pdf)
  ggsave(file.path(PANEL_DIR, "7Q.png"), P_7Q, width = 95, height = 95, units = "mm", dpi = 300)
} else {
  cat("[skip] 7Q: aggregate p-value CSV not found at", roi_pv_path, "\n")
}

cat("[done] Fig7 v3 panels saved -> ", PANEL_DIR, "\n")
