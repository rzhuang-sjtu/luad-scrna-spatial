#!/usr/bin/env Rscript
# Tighter margins + remove panel-letter from subtitle. Keep sample name only.
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(patchwork); library(png); library(grid); library(jsonlite)
})

DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT_MAIN <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
OUT_S11  <- "${WORK_ROOT}/luad_figures/fig_s11"
SECTIONS_DIR <- file.path(OUT_S11, "spatial_sections_extended")
dir.create(SECTIONS_DIR, showWarnings = FALSE, recursive = TRUE)
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")
GENE_HIGH <- c(SEC61G = "#FFB300", SRSF9 = "#00BCD4", ANGPTL4 = "#43A047")

spot <- read_csv(file.path(DATA, "8I_spot_long.csv"), show_col_types = FALSE)
spot$gene <- factor(spot$gene, levels = LEADS)

HE_ROOTS <- list(
  `E-MTAB-13530` = list(
    dir = "${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530",
    sub = function(s) sprintf("%s-spatial", s)
  ),
  Okamura = list(
    dir = "${DATA_ROOT}/ST/results/step09_okamura_validation/raw",
    sub = function(s) sprintf("%s/spatial", s)
  )
)
load_he_section <- function(cohort, section) {
  cfg <- HE_ROOTS[[cohort]]; if (is.null(cfg)) return(NULL)
  base <- file.path(cfg$dir, cfg$sub(section))
  img_p <- file.path(base, "tissue_hires_image.png")
  sf_p  <- file.path(base, "scalefactors_json.json")
  if (!file.exists(img_p) || !file.exists(sf_p)) return(NULL)
  img <- png::readPNG(img_p)
  sf  <- jsonlite::fromJSON(sf_p)$tissue_hires_scalef
  list(img = img, w = dim(img)[2], h = dim(img)[1], sf = sf)
}

plot_gene_he <- function(d, he, gene, high_color,
                          point_size = 0.85, alpha_max = 0.92,
                          subtitle_size = 7) {
  v <- d$expr
  lo <- quantile(v, 0.02, na.rm = TRUE)
  hi <- quantile(v, 0.98, na.rm = TRUE)
  rng <- hi - lo
  a_norm <- pmin(1, pmax(0, (v - lo) / ifelse(rng > 0, rng, 1)))
  d$alpha_pt <- alpha_max * a_norm
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d, aes(x = x_he, y = y_he, color = expr, alpha = alpha_pt),
               size = point_size, stroke = 0, shape = 16) +
    scale_color_gradient(low = "white", high = high_color,
                         limits = c(lo, hi), oob = scales::squish,
                         name = sprintf("%s\nexpression", gene),
                         guide = guide_colorbar(barwidth = unit(2.2, "mm"),
                                                barheight = unit(16, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.25,
                                                title.position = "top")) +
    scale_alpha_identity() +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_pub(8) +
    theme(axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          axis.line     = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          panel.spacing = unit(0, "mm"),
          plot.margin   = margin(1, 1, 1, 1),
          legend.position = "right",
          legend.margin   = margin(0, 0, 0, 1),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.title  = element_text(size = 7, family = FAM, face = "bold",
                                       color = "black"),
          legend.text   = element_text(size = 6, family = FAM, color = "black"),
          plot.subtitle = element_text(size = subtitle_size, color = "black",
                                       face = "bold", margin = margin(b = 1))) +
    labs(subtitle = gene)
}

draw_section <- function(cohort, section, kind, out_stem) {
  sub <- spot %>% filter(cohort == !!cohort, sample == section)
  if (!nrow(sub)) { cat(sprintf("  WARN: no spots %s/%s\n", cohort, section)); return(invisible(NULL)) }
  he <- load_he_section(cohort, section)
  if (is.null(he)) { cat(sprintf("  WARN: no H&E %s/%s\n", cohort, section)); return(invisible(NULL)) }
  sub <- sub %>% mutate(x_he = spatial1 * he$sf, y_he = -spatial2 * he$sf)
  ps <- lapply(LEADS, function(g) {
    plot_gene_he(sub %>% filter(gene == g), he, g, unname(GENE_HIGH[g]))
  })
  cohort_disp <- ifelse(cohort == "Okamura", "Takano 2024", cohort)
  # subtitle: cohort | sample | kind  (no panel letter)
  hdr <- sprintf("%s   %s   %s", cohort_disp, section, kind)
  p <- wrap_plots(ps, nrow = 1) +
        plot_layout(guides = "keep") +
        plot_annotation(
          subtitle = hdr,
          theme = theme(plot.subtitle = element_text(size = 7.5, family = FAM,
                                                     face = "bold", hjust = 0.02,
                                                     margin = margin(b = 1)),
                        plot.margin = margin(1, 1, 1, 1)))
  # tighter canvas: 6.0 x 2.1 (was 6.6 x 2.4)
  save_panel(p, out_stem, 6.0, 2.1)
  cat(sprintf("  saved %s\n", basename(out_stem)))
}

# === MAIN figure: 8I-L ===
cat("\n[MAIN]\n")
draw_section("E-MTAB-13530", "P10_T1",     "R-surrogate",  file.path(OUT_MAIN, "8I_spatial_P10_T1"))
draw_section("E-MTAB-13530", "P15_T1",     "NR-surrogate", file.path(OUT_MAIN, "8J_spatial_P15_T1"))
draw_section("Okamura",       "LUAD_No_4", "R-surrogate (Takano cross-cohort)", file.path(OUT_MAIN, "8K_spatial_LUAD_No_4"))
draw_section("Okamura",       "LUAD_No_1", "NR-surrogate (Takano cross-cohort)", file.path(OUT_MAIN, "8L_spatial_LUAD_No_1"))

# === S11: ALL extra sections ===
EM_MAIN  <- c("P10_T1","P15_T1")
TK_MAIN  <- c("LUAD_No_4","LUAD_No_1")
em_all <- sort(unique(spot$sample[spot$cohort == "E-MTAB-13530"]))
tk_all <- sort(unique(spot$sample[spot$cohort == "Okamura"]))
em_other <- setdiff(em_all, EM_MAIN)
tk_other <- setdiff(tk_all, TK_MAIN)

cat(sprintf("\n[S11] E-MTAB extended (%d):\n", length(em_other)))
for (s in em_other) {
  draw_section("E-MTAB-13530", s, "extended (E-MTAB)",
                file.path(SECTIONS_DIR, sprintf("EMTAB_%s", s)))
}
cat(sprintf("\n[S11] Takano extended (%d):\n", length(tk_other)))
for (s in tk_other) {
  draw_section("Okamura", s, "extended (Takano)",
                file.path(SECTIONS_DIR, sprintf("Takano_%s", s)))
}
cat("\nDONE.\n")
