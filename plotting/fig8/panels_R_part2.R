#!/usr/bin/env Rscript
# Fig 8 v2 — spatial panels: 8H co-expression, 8I/J/K/L (R/NR sections × 3 genes),
# S11A/B/C (extra sections). Spatial panels follow Fig 7 stepB_fig7_plot_v3.R style:
# H&E underlay via annotation_raster, per-section tissue_hires_scalef, white→signature
# color gradient, alpha proportional to value.

suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(png); library(grid); library(jsonlite); library(scales)
})

DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

# Per-cohort H&E source dirs (different layouts).
HE_ROOTS <- list(
  `E-MTAB-13530` = list(dir = "${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530",
                        sub = function(s) paste0(s, "-spatial")),
  Okamura        = list(dir = "${DATA_ROOT}/ST/results/r_data/he",
                        sub = function(s) paste0("Okamura__", s))
)

# Per-gene high color: 3 maximally-distinct hues, all high-contrast against the
# magenta-pink H&E background (avoid red family). Lead (SEC61G) = gold for max pop.
GENE_HIGH <- c(SEC61G = "#FFB300", SRSF9 = "#00BCD4", ANGPTL4 = "#43A047")

cat("\n[8H] co-expression embedding\n")
ce  <- read_csv(file.path(DATA, "8H_coexpr_embedding.csv"), show_col_types = FALSE)
ei  <- read_csv(file.path(DATA, "8H_embedding_info.csv"),    show_col_types = FALSE)
ce$gene <- factor(ce$gene, levels = LEADS)
# subsample to keep PDF size sane (~8k spots per gene panel)
set.seed(42)
ce_sub <- ce %>% group_by(gene) %>%
  group_modify(~ if (nrow(.x) > 8000) slice_sample(.x, n = 8000) else .x) %>%
  ungroup()
p_8H <- ggplot(ce_sub, aes(x, y, color = expr)) +
  geom_point(size = 0.25, alpha = 0.7, shape = 16) +
  scale_st_color(name = "expr") +
  facet_wrap(~ gene, nrow = 1) +
  coord_fixed() +
  labs(x = sprintf("%s-1", ei$embedding[1]), y = sprintf("%s-2", ei$embedding[1]),
       subtitle = sprintf("%s spot embedding · n=%d (subsampled to 8k/gene)",
                          ei$embedding[1], length(unique(ce$x)))) +
  theme_pub(8) +
  theme(legend.position = "right",
        legend.key.height = unit(0.45, "cm"),
        legend.key.width  = unit(0.18, "cm"))
save_panel(p_8H, file.path(OUT, "8H_coexpr_embedding"), 6.0, 2.4)

cat("\n[spatial] full coverage: per-section representatives + per-cohort atlases\n")
library(patchwork)
spot   <- read_csv(file.path(DATA, "8I_spot_long.csv"),         show_col_types = FALSE)
panels <- read_csv(file.path(DATA, "8I_panel_assignments.csv"), show_col_types = FALSE)
spot$gene <- factor(spot$gene, levels = LEADS)

# Per-section H&E + scalefactor loader. Returns list(img, w, h, sf) or NULL.
load_he_section <- function(cohort, section) {
  cfg <- HE_ROOTS[[cohort]]
  if (is.null(cfg)) return(NULL)
  base <- file.path(cfg$dir, cfg$sub(section))
  img_p <- file.path(base, "tissue_hires_image.png")
  sf_p  <- file.path(base, "scalefactors_json.json")
  if (!file.exists(img_p) || !file.exists(sf_p)) return(NULL)
  img <- png::readPNG(img_p)
  sf  <- jsonlite::fromJSON(sf_p)$tissue_hires_scalef
  list(img = img, w = dim(img)[2], h = dim(img)[1], sf = sf)
}

# Single-gene H&E spatial map (Fig 7 plot_spatial_he pattern):
# white -> signature color gradient, alpha proportional to value, image row-0 at top
# via negative y_he so coord_fixed runs in natural order.
plot_gene_he <- function(d, he, gene, high_color,
                         value_clip = c(0.02, 0.98),
                         fixed_lim = NULL,
                         point_size = 0.85, alpha_max = 0.92,
                         show_legend = TRUE, subtitle_size = 7,
                         legend_name = NULL) {
  v <- d$expr
  if (!is.null(fixed_lim)) {
    lo <- fixed_lim[1]; hi <- fixed_lim[2]
  } else {
    lo <- quantile(v, value_clip[1], na.rm = TRUE)
    hi <- quantile(v, value_clip[2], na.rm = TRUE)
  }
  rng <- hi - lo
  a_norm <- pmin(1, pmax(0, (v - lo) / ifelse(rng > 0, rng, 1)))
  d$alpha_pt <- alpha_max * a_norm
  if (is.null(legend_name)) legend_name <- sprintf("%s\nexpression", gene)
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d, aes(x = x_he, y = y_he, color = expr, alpha = alpha_pt),
               size = point_size, stroke = 0, shape = 16) +
    scale_color_gradient(low = "white", high = high_color,
                         limits = c(lo, hi), oob = scales::squish,
                         name = legend_name,
                         guide = if (show_legend)
                           guide_colorbar(barwidth  = unit(3, "mm"),
                                          barheight = unit(28, "mm"),
                                          frame.colour = "black",
                                          frame.linewidth = 0.3,
                                          ticks.colour = "black",
                                          title.position = "top")
                         else "none") +
    scale_alpha_identity() +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_pub(8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          legend.position = "right",
          legend.title = element_text(size = 9, family = FAM, face = "bold",
                                       color = "black", margin = margin(b = 2)),
          legend.text  = element_text(size = 7, family = FAM, color = "black"),
          plot.subtitle = element_text(size = subtitle_size, color = "black",
                                       face = "bold", margin = margin(b = 1.5))) +
    labs(subtitle = gene)
}

# ---- (1) Individual representative panels: 1 section x 3 genes per PDF ----
draw_section <- function(cohort, section, panel_label, kind, out_stem) {
  sub <- spot %>% filter(cohort == !!cohort, sample == section)
  if (!nrow(sub)) { cat(sprintf("  WARN: no spots for %s/%s\n", cohort, section)); return(invisible(NULL)) }
  he <- load_he_section(cohort, section)
  if (is.null(he)) {
    cat(sprintf("  WARN: no H&E for %s/%s, skipping\n", cohort, section)); return(invisible(NULL))
  }
  sub <- sub %>% mutate(x_he = spatial1 * he$sf, y_he = -spatial2 * he$sf)
  ps <- lapply(LEADS, function(g) {
    plot_gene_he(sub %>% filter(gene == g), he, g, unname(GENE_HIGH[g]))
  })
  p <- wrap_plots(ps, nrow = 1) +
        plot_annotation(
          subtitle = sprintf("%s   %s   %s   %s", panel_label, cohort, section, kind),
          theme = theme(plot.subtitle = element_text(size = 8, family = FAM,
                                                     face = "bold", hjust = 0.02,
                                                     margin = margin(b = 2))))
  save_panel(p, out_stem, 6.6, 2.4)
}

cat("\n[part2.1] individual representative panels (8I-L, S11A-C, S11I-J)\n")
for (i in seq_len(nrow(panels))) {
  draw_section(panels$cohort[i], panels$sample[i], panels$panel[i], panels$kind[i],
               file.path(OUT, sprintf("%s_spatial_%s", panels$panel[i], panels$sample[i])))
  cat(sprintf("  %s - %s/%s (%s)\n", panels$panel[i], panels$cohort[i],
              panels$sample[i], panels$kind[i]))
}

# ---- (2) Per-cohort, per-gene full atlas ----
# Each atlas = grid of all sections in a cohort. We share one colour scale
# across the cohort (cohort-wide quantiles) and merge per-panel legends into
# one shared colorbar with `plot_layout(guides = "collect")`.
make_atlas <- function(cohort, gene, ncol_grid, panel_w_in, panel_h_in, out_stem) {
  high <- unname(GENE_HIGH[gene])
  sections <- sort(unique(spot$sample[spot$cohort == cohort]))
  # cohort-wide quantile clip so every panel uses the same colour scale
  vals <- spot$expr[spot$cohort == cohort & spot$gene == gene]
  vals <- vals[is.finite(vals)]
  fixed_lim <- if (length(vals)) {
    as.numeric(quantile(vals, c(0.02, 0.98), na.rm = TRUE))
  } else NULL

  cat(sprintf("  atlas %s/%s : %d sections (shared range %.2f-%.2f)\n",
              cohort, gene, length(sections),
              if (!is.null(fixed_lim)) fixed_lim[1] else NA,
              if (!is.null(fixed_lim)) fixed_lim[2] else NA))

  ps <- list()
  for (s in sections) {
    he <- load_he_section(cohort, s)
    if (is.null(he)) {
      cat(sprintf("    skip %s (no H&E)\n", s)); next
    }
    sub <- spot %>% filter(cohort == !!cohort, sample == s, gene == !!gene)
    if (!nrow(sub)) next
    sub <- sub %>% mutate(x_he = spatial1 * he$sf, y_he = -spatial2 * he$sf)
    p_one <- plot_gene_he(sub, he, s, high,
                          fixed_lim = fixed_lim,
                          point_size = 0.95, alpha_max = 0.92,
                          show_legend = TRUE,
                          legend_name = sprintf("%s\nexpression", gene),
                          subtitle_size = 6.5) +
             theme(plot.margin = margin(1, 1, 1, 1))
    ps[[length(ps) + 1]] <- p_one
  }
  if (!length(ps)) return(invisible(NULL))
  nr <- ceiling(length(ps) / ncol_grid)
  comp <- wrap_plots(ps, ncol = ncol_grid, nrow = nr) +
          plot_layout(guides = "collect") +
          plot_annotation(
            subtitle = sprintf(
              "%s atlas - %s expression across %d sections",
              ifelse(cohort == "Okamura", "Takano 2024", cohort),
              gene, length(ps)),
            theme = theme(plot.subtitle = element_text(size = 9, family = FAM,
                                                       face = "bold", hjust = 0.02,
                                                       margin = margin(b = 3)))) &
          theme(legend.position = "right",
                legend.justification = "center",
                legend.title = element_text(size = 9, family = FAM, face = "bold",
                                             color = "black",
                                             margin = margin(b = 2)),
                legend.text  = element_text(size = 7, family = FAM, color = "black"),
                legend.key.height = unit(28, "mm"),
                legend.key.width  = unit(3, "mm"),
                legend.box.margin = margin(0, 0, 0, 2))
  # Add ~0.7" to width to host the shared legend column.
  save_panel(comp, out_stem, ncol_grid * panel_w_in + 0.7, nr * panel_h_in + 0.5)
}

cat("\n[part2.2] per-cohort, per-gene atlases (S11K-S11P)\n")
# Atlas labels: S11K/L/M = E-MTAB-13530 SEC61G/SRSF9/ANGPTL4; S11N/O/P = Takano same
# Filename uses display cohort name ("Takano") not the internal key ("Okamura")
# so the supplementary figure naming matches the published cohort label.
atlas_specs <- list(
  list(label = "S11K", cohort = "E-MTAB-13530", display = "EMTAB13530", gene = "SEC61G",  ncol = 4, w = 1.6, h = 1.6),
  list(label = "S11L", cohort = "E-MTAB-13530", display = "EMTAB13530", gene = "SRSF9",   ncol = 4, w = 1.6, h = 1.6),
  list(label = "S11M", cohort = "E-MTAB-13530", display = "EMTAB13530", gene = "ANGPTL4", ncol = 4, w = 1.6, h = 1.6),
  list(label = "S11N", cohort = "Okamura",      display = "Takano",     gene = "SEC61G",  ncol = 4, w = 1.6, h = 1.6),
  list(label = "S11O", cohort = "Okamura",      display = "Takano",     gene = "SRSF9",   ncol = 4, w = 1.6, h = 1.6),
  list(label = "S11P", cohort = "Okamura",      display = "Takano",     gene = "ANGPTL4", ncol = 4, w = 1.6, h = 1.6)
)
for (sp in atlas_specs) {
  make_atlas(sp$cohort, sp$gene, sp$ncol, sp$w, sp$h,
             file.path(OUT, sprintf("%s_atlas_%s_%s",
                                    sp$label, sp$display, sp$gene)))
}

cat("\nALL part2 panels DONE.\n")
