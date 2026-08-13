#!/usr/bin/env Rscript
# Compose Figure S7 (8 panels) into a single PDF/PNG via PNG raster.
# Heatmap panels (S7D/F) come from ComplexHeatmap and don't easily live in
# patchwork, so we read each panel's saved PNG and arrange via cowplot.
suppressPackageStartupMessages({
  library(png); library(grid); library(ggplot2); library(patchwork)
})

BASE <- if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures" else
  "${WORK_ROOT}/luad_figures"

panel_from_png <- function(path) {
  if (!file.exists(path)) return(patchwork::plot_spacer())
  img <- png::readPNG(path)
  g <- grid::rasterGrob(img, width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"),
                        interpolate = TRUE)
  ggplot() + annotation_custom(g) + theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2, "pt"))
}

panels <- list(
  panel_from_png(file.path(BASE, "fig_s7/figS7a_umap_tissue.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7b_umap_dataset.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7c_scanvi_violin.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7d_confusion.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7e_mp1_categories.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7f_tcga_partial.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7g_celltype_bar.png")),
  panel_from_png(file.path(BASE, "fig_s7/figS7h_km_3panels.png"))
)

combo <- wrap_plots(panels, ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

out <- file.path(BASE, "fig_s7", "FigureS7.pdf")
ggsave(out, combo, width = 360, height = 480, units = "mm",
       limitsize = FALSE)
ggsave(sub("\\.pdf$", ".png", out), combo,
       width = 360, height = 480, units = "mm", dpi = 300, limitsize = FALSE)
cat(" Saved", out, "\n")
