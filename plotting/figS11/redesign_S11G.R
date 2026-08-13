#!/usr/bin/env Rscript
# Redesign S11G: 2-row composite, top = 2 cohort H&E side-by-side,
# bottom = 7-gene x 3-ROI delta heatmap (cohorts as facet columns).
# Canvas 7.0 x 5.5 in (square-ish, top journal style).
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(png); library(grid); library(patchwork); library(jsonlite)
})

DATA <- "${DATA_ROOT}/ST/results/step_invasive_front"
OUT_S11 <- "${WORK_ROOT}/luad_figures/fig_s11"

ROI_PRETTY <- c(stromal_immune = "Stromal-immune",
                tumor_intrinsic = "Tumor-intrinsic",
                invasive_front  = "Compositional-mixing")
ROI_COLS <- c(`Stromal-immune` = "#3C5488",
              `Tumor-intrinsic` = "#E64B35",
              `Compositional-mixing`  = "#F39B7F",
              Multiple = "#2A2A2A",
              None = "grey80")

HE_ROOT <- "${DATA_ROOT}/ST/results/r_data/he"
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

build_overlay <- function(dataset, section) {
  spots <- read_csv(file.path(DATA, sprintf("spots_%s.csv", dataset)),
                    show_col_types = FALSE) %>% filter(sample == section)
  if (!nrow(spots)) return(NULL)
  he <- load_he(dataset, section)
  if (is.null(he)) return(NULL)
  spots <- spots %>%
    mutate(n_hits = stromal_immune_roi + tumor_intrinsic_roi + invasive_front_roi,
           roi_state = case_when(
             n_hits >= 2                  ~ "Multiple",
             tumor_intrinsic_roi          ~ "Tumor-intrinsic",
             invasive_front_roi           ~ "Compositional-mixing",
             stromal_immune_roi           ~ "Stromal-immune",
             TRUE                         ~ "None"),
           x_he =  x * he$sf,
           y_he = -y * he$sf)
  spots$roi_state <- factor(spots$roi_state,
                            levels = c("None","Stromal-immune","Tumor-intrinsic",
                                       "Compositional-mixing","Multiple"))
  spot_alphas <- c(None = 0.25, `Stromal-immune` = 0.95,
                   `Tumor-intrinsic` = 0.95, `Compositional-mixing` = 0.95,
                   Multiple = 0.95)
  cohort_disp <- ifelse(dataset == "Okamura", "Takano 2024", dataset)
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = subset(spots, roi_state == "None"),
               aes(x_he, y_he), color = "grey80", size = 0.18, shape = 16, alpha = 0.30) +
    geom_point(data = subset(spots, roi_state != "None"),
               aes(x_he, y_he, color = roi_state, alpha = roi_state),
               size = 0.55, shape = 16) +
    scale_color_manual(values = ROI_COLS, name = NULL,
                       breaks = c("Stromal-immune","Tumor-intrinsic",
                                  "Compositional-mixing","Multiple"),
                       guide = guide_legend(override.aes = list(size = 2.0))) +
    scale_alpha_manual(values = spot_alphas, guide = "none") +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    labs(x = NULL, y = NULL,
         subtitle = sprintf("%s   %s", cohort_disp, section)) +
    theme_pub(8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position   = "bottom",
          legend.text       = element_text(size = 6.5),
          legend.key.size   = unit(2.5, "mm"),
          legend.margin     = margin(0, 0, 0, 0),
          legend.box.margin = margin(-2, 0, 0, 0),
          panel.background  = element_rect(fill = "white", color = NA),
          plot.subtitle     = element_text(size = 7.5, face = "bold",
                                            margin = margin(b = 1)),
          plot.margin       = margin(2, 2, 1, 2))
}

build_heatmap_horiz <- function() {
  g <- read_csv(file.path(DATA, "gene_by_roi_stats.csv"), show_col_types = FALSE)
  g$roi_type <- factor(ROI_PRETTY[g$roi_type], levels = unname(ROI_PRETTY))
  g$gene <- factor(g$gene,
                    levels = c("OSM","IL1B","FOSB","ATF3","ANGPTL4","SRSF9","SEC61G"))
  g$dataset <- factor(g$dataset, levels = c("E-MTAB-13530","Okamura"),
                                    labels = c("E-MTAB-13530","Takano 2024"))
  g$star <- sapply(g$mw_p, sig_stars)
  g$lbl <- sprintf("%+.2f%s", g$delta_roi_minus_non,
                    ifelse(g$star == "ns", "", g$star))
  mx <- max(abs(g$delta_roi_minus_non), na.rm = TRUE)
  ggplot(g, aes(gene, roi_type, fill = delta_roi_minus_non)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = lbl), family = FAM, size = 1.9) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                         midpoint = 0, name = "Delta(ROI - non)",
                         limits = c(-mx, mx),
                         breaks = pretty(c(-mx, mx), 5),
                         guide = guide_colorbar(barwidth  = unit(2.5, "mm"),
                                                barheight = unit(15, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.3,
                                                title.position = "top")) +
    facet_wrap(~ dataset, ncol = 1, strip.position = "right") +
    labs(x = NULL, y = NULL,
         subtitle = "7 genes x 3 ROIs - Delta(mean expr in ROI - non);  *p<0.05  **<0.01  ***<0.001  ****<1e-4") +
    theme_pub(8) +
    theme(axis.text.x   = element_text(face = "bold", size = 7.5),
          axis.text.y   = element_text(size = 7),
          strip.text    = element_text(size = 7.5, face = "bold"),
          strip.background = element_rect(fill = "grey94", color = NA),
          legend.position  = "right",
          legend.title  = element_text(size = 7, face = "bold"),
          legend.text   = element_text(size = 6.5),
          plot.subtitle = element_text(size = 7, color = "grey25",
                                        margin = margin(b = 1)),
          plot.margin   = margin(2, 2, 1, 2))
}

# Build sections
rep <- read_csv(file.path(DATA, "representative_sections.csv"), show_col_types = FALSE)
sec_E <- rep[["section"]][rep[["dataset"]] == "E-MTAB-13530"][1]
sec_O <- rep[["section"]][rep[["dataset"]] == "Okamura"][1]
ov_E <- build_overlay("E-MTAB-13530", sec_E)
ov_O <- build_overlay("Okamura",      sec_O)
hm   <- build_heatmap_horiz()

# 2-row composite: top row = 2 H&E side-by-side; bottom = heatmap full-width
top    <- (ov_E | ov_O) + plot_layout(guides = "collect", widths = c(1, 1)) &
          theme(legend.position = "bottom")
bottom <- hm
composite <- (top / bottom) + plot_layout(heights = c(1.2, 1))

save_panel(composite,
           file.path(OUT_S11, "S11G_invasive_front_composite"),
           7.0, 5.6)
cat("S11G refreshed (7.0 x 5.6 in, 2-row balanced layout)\n")
