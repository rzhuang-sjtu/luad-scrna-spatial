#!/usr/bin/env Rscript
# Three-ROI analysis: the third region is labelled Compositional-mixing in the figure. Directory and column names keep the older invasive_front spelling, which is what the stored tables contain.
# Composite S11H panel + individual sub-panels.
# Composite S11H: left = H&E + 3-ROI overlay for representative section per cohort
#                 (E-MTAB on top, Okamura below); right = 7-gene × 3-ROI Δ heatmap
#                 faceted by cohort.
# Individual sub-panels saved to step_invasive_front/plots/ as backup data.

suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(png); library(grid); library(patchwork); library(jsonlite)
})

DATA <- "${DATA_ROOT}/ST/results/step_invasive_front"
OUT  <- file.path(DATA, "plots")
OUT_S11 <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ROI_PRETTY <- c(stromal_immune = "Stromal-immune",
                tumor_intrinsic = "Tumor-intrinsic",
                invasive_front  = "Compositional-mixing")
ROI_COLS <- c(`Stromal-immune` = "#3C5488",
              `Tumor-intrinsic` = "#E64B35",
              `Compositional-mixing`  = "#F39B7F",
              Multiple = "#2A2A2A",
              None = "grey80")
HE_ROOT <- "${DATA_ROOT}/ST/results/r_data/he"

# Load H&E image + per-section tissue_hires_scalef (mirrors fig7 stepB pattern)
load_he <- function(dataset, section) {
  prefix <- if (dataset == "E-MTAB-13530") "EMTAB13530" else "Okamura"
  base <- file.path(HE_ROOT, paste0(prefix, "__", section))
  img_p <- file.path(base, "tissue_hires_image.png")
  sf_p  <- file.path(base, "scalefactors_json.json")
  if (!file.exists(img_p) || !file.exists(sf_p)) return(NULL)
  img <- png::readPNG(img_p)
  sf  <- jsonlite::fromJSON(sf_p)$tissue_hires_scalef
  list(img = img, w = dim(img)[2], h = dim(img)[1], sf = sf)
}

# ---- helper: build overlay ggplot for one section (V1 H&E-underlay style) ----
build_overlay <- function(dataset, section, show_x = FALSE) {
  spots <- read_csv(file.path(DATA, sprintf("spots_%s.csv", dataset)),
                    show_col_types = FALSE) %>% filter(sample == section)
  if (!nrow(spots)) return(NULL)
  he <- load_he(dataset, section)
  if (is.null(he)) {
    cat(sprintf("  WARN: no H&E found for %s/%s\n", dataset, section))
    return(NULL)
  }
  spots <- spots %>%
    mutate(
      n_hits = stromal_immune_roi + tumor_intrinsic_roi + invasive_front_roi,
      roi_state = case_when(
        n_hits >= 2                                ~ "Multiple",
        tumor_intrinsic_roi                        ~ "Tumor-intrinsic",
        invasive_front_roi                         ~ "Compositional-mixing",
        stromal_immune_roi                         ~ "Stromal-immune",
        TRUE                                       ~ "None"),
      x_he =  x * he$sf,
      y_he = -y * he$sf
    )
  spots$roi_state <- factor(spots$roi_state,
                            levels = c("None","Stromal-immune","Tumor-intrinsic",
                                       "Compositional-mixing","Multiple"))
  spot_alphas <- c(None = 0.18, `Stromal-immune` = 0.95,
                   `Tumor-intrinsic` = 0.95, `Compositional-mixing` = 0.95,
                   Multiple = 0.95)
  n_per_roi <- spots %>% count(roi_state) %>%
    filter(roi_state != "None") %>%
    arrange(match(roi_state, c("Stromal-immune","Tumor-intrinsic","Compositional-mixing","Multiple")))
  roi_label <- paste(sprintf("%s n=%d", n_per_roi$roi_state, n_per_roi$n), collapse = " | ")
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = subset(spots, roi_state == "None"),
               aes(x_he, y_he, color = roi_state, alpha = roi_state),
               size = 0.22, shape = 16) +
    geom_point(data = subset(spots, roi_state != "None"),
               aes(x_he, y_he, color = roi_state, alpha = roi_state),
               size = 0.55, shape = 16) +
    scale_color_manual(values = ROI_COLS, name = NULL,
                       breaks = c("Stromal-immune","Tumor-intrinsic",
                                  "Compositional-mixing","Multiple")) +
    scale_alpha_manual(values = spot_alphas, guide = "none") +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    labs(x = NULL, y = NULL,
         subtitle = sprintf("%s   %s   %s",
                            ifelse(dataset == "Okamura", "Takano 2024", dataset),
                            section, roi_label)) +
    theme_pub(8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_blank())
}

# ---- helper: build heatmap ----
build_heatmap <- function() {
  g <- read_csv(file.path(DATA, "gene_by_roi_stats.csv"), show_col_types = FALSE)
  g$roi_type <- factor(ROI_PRETTY[g$roi_type], levels = unname(ROI_PRETTY))
  g$gene <- factor(g$gene, levels = c("OSM","IL1B","FOSB","ATF3","ANGPTL4","SRSF9","SEC61G"))
  g$dataset <- factor(g$dataset, levels = c("E-MTAB-13530", "Okamura"),
                                   labels = c("E-MTAB-13530", "Takano 2024"))
  g$star <- sapply(g$mw_p, sig_stars)
  g$lbl  <- sprintf("%+.2f\n%s", g$delta_roi_minus_non, g$star)
  mx <- max(abs(g$delta_roi_minus_non), na.rm = TRUE)
  ggplot(g, aes(roi_type, gene, fill = delta_roi_minus_non)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = lbl), family = FAM, size = 1.95, lineheight = 0.85) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                         midpoint = 0, name = "Δ(ROI − non)",
                         limits = c(-mx, mx), breaks = pretty(c(-mx, mx), 5)) +
    facet_wrap(~ dataset, ncol = 1, strip.position = "top") +
    labs(x = NULL, y = NULL,
         subtitle = "7 genes × 3 ROIs · Δ(mean expr in ROI − non) · stars: * p<0.05  ** <1e-2  *** <1e-3  **** <1e-4") +
    theme_pub(8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          legend.key.height = unit(0.45, "cm"),
          legend.key.width  = unit(0.18, "cm"),
          strip.text = element_text(size = 7, face = "bold"))
}

cat("\n[S11H composite] H&E + ROI overlay (per cohort) | gene × ROI heatmap\n")
rep <- read_csv(file.path(DATA, "representative_sections.csv"), show_col_types = FALSE)
# Use [[ ]] vector extraction — chained `$col[mask]` on a readr tibble
# silently picks the wrong column in this readr/tibble version.
sec_E <- rep[["section"]][rep[["dataset"]] == "E-MTAB-13530"][1]
sec_O <- rep[["section"]][rep[["dataset"]] == "Okamura"][1]
ov_E <- build_overlay("E-MTAB-13530", sec_E)
ov_O <- build_overlay("Okamura",      sec_O)
hm   <- build_heatmap()

left  <- (ov_E / ov_O) + plot_layout(guides = "collect")
right <- hm
composite <- (left | right) +
  plot_layout(widths = c(1.2, 1.0)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = NA)))

save_panel(composite, file.path(OUT_S11, "S11H_invasive_front_composite"), 11.5, 6.0)
save_panel(composite, file.path(OUT,     "S11H_invasive_front_composite"), 11.5, 6.0)

cat("\n[backup] individual sub-panels\n")
save_panel(ov_E, file.path(OUT, sprintf("overlay_E-MTAB-13530_%s", sec_E)), 6, 4)
save_panel(ov_O, file.path(OUT, sprintf("overlay_Takano_%s", sec_O)), 6, 4)
save_panel(hm,   file.path(OUT, "gene_x_ROI_heatmap"), 6.5, 5.0)

# co-occurrence curves panel (extra, not in composite)
cat("\n[backup] co-occurrence curves\n")
co <- read_csv(file.path(DATA, "cooccurrence_long.csv"), show_col_types = FALSE)
co$pair <- factor(paste0(co$anchor, "→", co$target),
                  levels = c("Malignant→Neutrophil",
                             "Malignant→Macro_SPP1",
                             "Malignant→Fibroblast"))
co$dataset <- factor(co$dataset, levels = c("E-MTAB-13530", "Okamura"),
                                  labels = c("E-MTAB-13530", "Takano 2024"))
mean_co <- co %>%
  group_by(dataset, pair, dist_mid) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE), .groups = "drop")
p_curves <- ggplot(co, aes(dist_mid, prob, group = sample)) +
  geom_line(aes(color = dataset), linewidth = 0.3, alpha = 0.35) +
  geom_line(data = mean_co, aes(dist_mid, mean_prob, group = 1),
            color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  facet_grid(dataset ~ pair, scales = "free_y") +
  scale_color_manual(values = c(`E-MTAB-13530` = "#3C5488",
                                `Takano 2024` = "#E64B35"),
                     guide = "none") +
  labs(x = "Distance from anchor (full-res pixels)",
       y = "P(target | anchor)",
       subtitle = "Per-section co-occurrence; thin = sections, thick = cohort mean") +
  theme_pub(8) +
  theme(strip.text = element_text(size = 7))
save_panel(p_curves, file.path(OUT, "cooccurrence_curves"), 7.5, 4.0)
save_panel(p_curves, file.path(OUT_S11, "S11H_cooccurrence_curves_supp"), 7.5, 4.0)

cat("\nDONE.\n")
