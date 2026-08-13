#!/usr/bin/env Rscript
# Combine Supplementary Figures S4 + S5 + S6 into one mega figure.
# Reads each panel PNG (saved by individual scripts) and lays them out
# in a single composite via patchwork + grid::rasterGrob.

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(ggplot2)
  library(patchwork)
})

BASE <- if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures" else
  "${WORK_ROOT}/luad_figures"

# Helper: load a PNG as a ggplot wrapping its rasterGrob
panel_from_png <- function(path, label = NULL) {
  if (!file.exists(path)) {
    warning("missing: ", path)
    return(patchwork::plot_spacer())
  }
  img <- png::readPNG(path)
  g <- grid::rasterGrob(img, width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"),
                        interpolate = TRUE)
  p <- ggplot() + annotation_custom(g) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2, "pt"))
  if (!is.null(label)) {
    p <- p + labs(tag = label) +
      theme(plot.tag = element_text(size = 12, face = "bold",
                                    family = if (capabilities("cairo")) "sans" else "sans"))
  }
  p
}

s4 <- list(
  panel_from_png(file.path(BASE, "fig_s4/fig_s4a_MP1.png"), "a"),
  panel_from_png(file.path(BASE, "fig_s4/fig_s4b_MP2.png"), "b"),
  panel_from_png(file.path(BASE, "fig_s4/fig_s4c_MP3.png"), "c"),
  panel_from_png(file.path(BASE, "fig_s4/fig_s4d_MP4.png"), "d")
)
S4 <- wrap_plots(s4, ncol = 2) +
  plot_annotation(title = "Figure S4 — Functional enrichment per MP",
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold")))

# ── Build S5 (4 NES + ES per MP) ──
# S5-1 produces fig_s5a_MP1..d_MP4.{pdf,png} — NES bar chart
# Use those as canonical S5 panels.
s5 <- list(
  panel_from_png(file.path(BASE, "fig_s5/fig_s5a_MP1.png"), "a"),
  panel_from_png(file.path(BASE, "fig_s5/fig_s5b_MP2.png"), "b"),
  panel_from_png(file.path(BASE, "fig_s5/fig_s5c_MP3.png"), "c"),
  panel_from_png(file.path(BASE, "fig_s5/fig_s5d_MP4.png"), "d")
)
S5 <- wrap_plots(s5, ncol = 2) +
  plot_annotation(title = "Figure S5 — Hallmark GSEA per MP",
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold")))

s6 <- list(
  panel_from_png(file.path(BASE, "fig_s6/figS6a_dotplot.png"), "a"),
  panel_from_png(file.path(BASE, "fig_s6/figS6b_macro_dotplot.png"), "b"),
  panel_from_png(file.path(BASE, "fig_s6/figS6c_general_go.png"), "c"),
  panel_from_png(file.path(BASE, "fig_s6/figS6d_prolif_go.png"), "d")
)
S6 <- wrap_plots(s6, ncol = 2) +
  plot_annotation(title = "Figure S6 — Macrophage subsets",
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold")))

# ── Save individual mega-files for each ──
out_dir <- file.path(BASE, "_combined")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "FigS4_combined.pdf"), S4, width = 16, height = 12)
ggsave(file.path(out_dir, "FigS4_combined.png"), S4, width = 16, height = 12, dpi = 300)
ggsave(file.path(out_dir, "FigS5_combined.pdf"), S5, width = 16, height = 12)
ggsave(file.path(out_dir, "FigS5_combined.png"), S5, width = 16, height = 12, dpi = 300)
ggsave(file.path(out_dir, "FigS6_combined.pdf"), S6, width = 16, height = 12)
ggsave(file.path(out_dir, "FigS6_combined.png"), S6, width = 16, height = 12, dpi = 300)

# ── Mega-supp: S4 + S5 + S6 stacked vertically ──
mega <- (S4) / (S5) / (S6) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(title = "Supplementary Figures S4–S6",
                  theme = theme(plot.title = element_text(size = 18,
                                                          face = "bold",
                                                          hjust = 0.5)))
ggsave(file.path(out_dir, "FigS4_S5_S6_combined.pdf"), mega,
       width = 16, height = 36, limitsize = FALSE)
ggsave(file.path(out_dir, "FigS4_S5_S6_combined.png"), mega,
       width = 16, height = 36, dpi = 300, limitsize = FALSE)

cat("\n=== Combined supp figures saved in:", out_dir, "===\n")
print(list.files(out_dir))
