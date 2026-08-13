#!/usr/bin/env Rscript
# Fig 7 / S10 full-cohort atlases — one PDF per (cohort, feature) showing all
# LUAD sections in a grid, big + lite versions. Mirrors Fig 8 v2_500/S11K-P
# convention: lite version embeds tissue_lowres_image.png (~10x smaller).
#
# Features covered (8): ct_Malignant, ct_Neu_Inflammatory, ct_Macro_SPP1,
# progeny_NFkB, gex_OSM, gex_IL1B, MP3_score, roi.
#
# Output: ${WORK_ROOT}/luad_figures/fig_s10/panels/atlas/

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(dplyr); library(tidyr); library(scales); library(grid)
  library(png); library(jsonlite)
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

# Arial setup (mirror stepB)
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p,
    bold = "~/.local/share/fonts/arialbd.ttf",
    italic = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300)
  FAM <- "Arial"
} else FAM <- "sans"

R_DATA   <- "${DATA_ROOT}/ST/results/r_data"
PER      <- file.path(R_DATA, "per_section")
HE_DIR   <- file.path(R_DATA, "he")
OUT      <- "${WORK_ROOT}/luad_figures/fig_s10/panels/atlas"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Per-cohort CSV stems and HE source dirs
COHORTS <- list(
  `E-MTAB-13530` = list(
    prefix   = "EMTAB13530",
    he_hires = file.path(HE_DIR, "%s"),       # tissue_hires + scalefactors live here
    he_lo    = "${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530/%s-spatial",  # raw, has lowres
    sections = NULL  # filled below
  ),
  Okamura = list(
    prefix   = "Okamura",
    he_hires = file.path(HE_DIR, "%s"),
    he_lo    = "${DATA_ROOT}/ST/results/step09_okamura_validation/raw/%s/spatial",
    sections = NULL
  )
)
# discover sections from per_section/ filenames
all_csv <- list.files(PER, pattern = "\\.csv$")
for (cn in names(COHORTS)) {
  pfx <- COHORTS[[cn]]$prefix
  secs <- sub(paste0("^", pfx, "__"), "",
              sub("\\.csv$", "", grep(paste0("^", pfx, "__"), all_csv, value = TRUE)))
  COHORTS[[cn]]$sections <- sort(secs)
}
cat("E-MTAB sections:", paste(COHORTS[["E-MTAB-13530"]]$sections, collapse=", "), "\n")
cat("Okamura sections:", paste(COHORTS[["Okamura"]]$sections, collapse=", "), "\n")

# ---- HE loader (hires or lowres) ----
load_he <- function(cohort_cfg, section, mode = c("hires", "lowres")) {
  mode <- match.arg(mode)
  base_h <- sprintf(cohort_cfg$he_hires, paste0(cohort_cfg$prefix, "__", section))
  base_l <- sprintf(cohort_cfg$he_lo, section)
  sf_p   <- file.path(base_h, "scalefactors_json.json")
  if (!file.exists(sf_p)) {
    sf_p <- file.path(base_l, "scalefactors_json.json")
    if (!file.exists(sf_p)) return(NULL)
  }
  sf_json <- jsonlite::fromJSON(sf_p)
  if (mode == "hires") {
    img_p <- file.path(base_h, "tissue_hires_image.png")
    sf <- sf_json$tissue_hires_scalef
  } else {
    img_p <- file.path(base_l, "tissue_lowres_image.png")
    sf <- sf_json$tissue_lowres_scalef
  }
  if (!file.exists(img_p)) return(NULL)
  img <- png::readPNG(img_p)
  list(img = img, w = dim(img)[2], h = dim(img)[1], sf = sf)
}

# ---- shared theme ----
theme_atlas <- function(base = 7) {
  theme_void(base_family = FAM, base_size = base) +
    theme(
      plot.subtitle = element_text(size = base, face = "bold", hjust = 0.5,
                                   color = "black", margin = margin(b = 1)),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(1, 1, 1, 1),
      # Keep legends *on* per panel so plot_layout(guides = "collect")
      # has something to merge into a single shared legend on the side.
      legend.position = "right"
    )
}

# ---- single-section single-feature plot ----
# `fixed_lim` (length-2 numeric) = cohort-wide colour limits computed by
# make_atlas(); when NULL each panel falls back to its own quantile range.
plot_one <- function(d, he, feature, opts, fixed_lim = NULL, legend_name = "") {
  v <- d[[feature]]
  d$x_he <-  d$spatial1 * he$sf
  d$y_he <- -d$spatial2 * he$sf
  if (opts$kind == "binary") {
    d$flag <- factor(ifelse(as.logical(v), "ROI", "non-ROI"),
                     levels = c("non-ROI", "ROI"))
    p <- ggplot() +
      annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
      geom_point(data = d,
                 aes(x_he, y_he, color = flag, size = flag, alpha = flag),
                 stroke = 0, shape = 16) +
      scale_color_manual(values = c("non-ROI" = "grey80",
                                    "ROI"     = opts$high_color),
                         name = legend_name,
                         guide = guide_legend(override.aes = list(size = 2))) +
      scale_size_manual(values = c("non-ROI" = opts$pt_bg,
                                   "ROI"     = opts$pt_roi),
                        guide = "none") +
      scale_alpha_manual(values = c("non-ROI" = 0.40,
                                    "ROI"     = 0.95),
                         guide = "none")
  } else if (opts$kind == "diverging") {
    if (!is.null(fixed_lim)) {
      lim <- fixed_lim
    } else {
      a <- max(abs(quantile(v, 0.05, na.rm=TRUE)), abs(quantile(v, 0.95, na.rm=TRUE)))
      lim <- c(-a, a)
    }
    a_norm <- pmin(1, abs(v) / max(lim[2], 1e-9))
    d$alpha_pt <- 0.92 * a_norm
    p <- ggplot() +
      annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
      geom_point(data = d, aes(x_he, y_he, color = .data[[feature]], alpha = alpha_pt),
                 size = opts$pt, stroke = 0, shape = 16) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = opts$high_color,
                            midpoint = 0, limits = lim, oob = scales::squish,
                            name = legend_name,
                            guide = guide_colorbar(barwidth = unit(3, "mm"),
                                                   barheight = unit(28, "mm"),
                                                   frame.colour = "black",
                                                   frame.linewidth = 0.3,
                                                   ticks.colour = "black",
                                                   title.position = "top")) +
      scale_alpha_identity()
  } else { # sequential
    if (!is.null(fixed_lim)) {
      lo <- fixed_lim[1]; hi <- fixed_lim[2]
    } else {
      lo <- quantile(v, 0.02, na.rm = TRUE); hi <- quantile(v, 0.98, na.rm = TRUE)
    }
    rng <- hi - lo
    a_norm <- pmin(1, pmax(0, (v - lo) / ifelse(rng > 0, rng, 1)))
    d$alpha_pt <- 0.92 * a_norm
    p <- ggplot() +
      annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
      geom_point(data = d, aes(x_he, y_he, color = .data[[feature]], alpha = alpha_pt),
                 size = opts$pt, stroke = 0, shape = 16) +
      scale_color_gradient(low = "white", high = opts$high_color,
                           limits = c(lo, hi), oob = scales::squish,
                           name = legend_name,
                           guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                  barheight = unit(20, "mm"),
                                                  frame.colour = "black",
                                                  frame.linewidth = 0.25,
                                                  ticks.colour = "black")) +
      scale_alpha_identity()
  }
  p + coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_atlas() + labs(subtitle = d$sample[1])
}

# ---- atlas builder ----
make_atlas <- function(cohort_name, feature, opts, mode = c("hires", "lowres"),
                       ncol_grid = 4, panel_w = 1.6, panel_h = 1.6, out_stem) {
  mode <- match.arg(mode)
  cfg <- COHORTS[[cohort_name]]

  # First pass: gather every numeric value across this cohort+feature so the
  # whole atlas shares one colour scale (otherwise wrap_plots collects N
  # separate colour-bars and the legend becomes useless).
  all_vals <- numeric(0)
  for (s in cfg$sections) {
    csv_p <- file.path(PER, paste0(cfg$prefix, "__", s, ".csv"))
    if (!file.exists(csv_p)) next
    d <- as.data.frame(fread(csv_p))
    if (feature %in% colnames(d)) {
      x <- suppressWarnings(as.numeric(d[[feature]]))
      all_vals <- c(all_vals, x[is.finite(x)])
    }
  }
  fixed_lim <- NULL
  if (opts$kind == "diverging" && length(all_vals)) {
    a <- max(abs(quantile(all_vals, 0.05, na.rm = TRUE)),
             abs(quantile(all_vals, 0.95, na.rm = TRUE)))
    fixed_lim <- c(-a, a)
  } else if (opts$kind == "sequential" && length(all_vals)) {
    fixed_lim <- as.numeric(quantile(all_vals, c(0.02, 0.98), na.rm = TRUE))
  }

  ps <- list()
  for (s in cfg$sections) {
    csv_p <- file.path(PER, paste0(cfg$prefix, "__", s, ".csv"))
    if (!file.exists(csv_p)) { cat(sprintf("    skip %s (no csv)\n", s)); next }
    d <- as.data.frame(fread(csv_p))
    if (!feature %in% colnames(d)) {
      cat(sprintf("    skip %s (no col %s)\n", s, feature)); next
    }
    he <- load_he(cfg, s, mode)
    if (is.null(he)) { cat(sprintf("    skip %s (no HE/sf for %s)\n", s, mode)); next }
    ps[[length(ps) + 1]] <- plot_one(d, he, feature, opts,
                                     fixed_lim = fixed_lim,
                                     legend_name = opts$pretty)
  }
  if (!length(ps)) return(invisible(NULL))
  nr <- ceiling(length(ps) / ncol_grid)
  comp <- wrap_plots(ps, ncol = ncol_grid, nrow = nr) +
          plot_layout(guides = "collect") +
          plot_annotation(
            subtitle = sprintf("%s atlas%s : %s across %d sections",
                               # Display Takano in legends/captions while
                               # keeping cohort_name="Okamura" as the
                               # internal key for file paths and h5ad fields.
                               ifelse(cohort_name == "Okamura", "Takano",
                                      cohort_name),
                               if (mode == "lowres") " (lowres)" else "",
                               opts$pretty, length(ps)),
            theme = theme(plot.subtitle = element_text(size = 9, family = FAM,
                                                       face = "bold", hjust = 0.02,
                                                       margin = margin(b = 3)))) &
          theme(legend.position = "right",
                legend.justification = "center",
                legend.title = element_text(size = 9, family = FAM, face = "bold",
                                             color = "black",
                                             margin = margin(b = 2)),
                legend.text  = element_text(size = 7, family = FAM, color = "black"),
                legend.key.size = unit(3.5, "mm"),
                legend.box.margin = margin(0, 0, 0, 2))
  for (ext in c("pdf", "png")) {
    # widen by ~0.6" to host the shared legend column
    ggsave(paste0(out_stem, ".", ext), comp,
           width = ncol_grid * panel_w + 0.7, height = nr * panel_h + 0.5,
           units = "in", dpi = 300, device = ext, bg = "white")
  }
  cat(sprintf("    %s -> %d sections, %d cols x %d rows\n", out_stem, length(ps), ncol_grid, nr))
}

# ---- feature catalog ----
# Colours and dot sizes deliberately match Fig 7 main script (stepB v3):
#   MP3 vivid purple = #651FFF (Material Deep Purple A400)
# Larger dots make the spatial signal readable at 1.6"-per-panel grid scale.
FEATS <- list(
  list(col = "ct_Malignant",        pretty = "Malignant abundance",
       kind = "sequential", high_color = "#C73E2A",
       opts = list(kind = "sequential", high_color = "#C73E2A", pt = 0.9)),
  list(col = "ct_Neu_Inflammatory", pretty = "Inflammatory neutrophil abundance",
       kind = "sequential", high_color = "#E64B35",
       opts = list(kind = "sequential", high_color = "#E64B35", pt = 0.9)),
  list(col = "ct_Macro_SPP1",       pretty = "SPP1+ macrophage abundance",
       kind = "sequential", high_color = "#F39B7F",
       opts = list(kind = "sequential", high_color = "#F39B7F", pt = 0.9)),
  list(col = "MP3_score",           pretty = "MP3 (EMT/IFN) score",
       kind = "sequential", high_color = "#651FFF",
       opts = list(kind = "sequential", high_color = "#651FFF", pt = 0.9)),
  list(col = "progeny_NFkB",        pretty = "PROGENy NFkB activity",
       kind = "diverging",  high_color = "#B2182B",
       opts = list(kind = "diverging",  high_color = "#B2182B", pt = 0.9)),
  list(col = "gex_OSM",             pretty = "OSM expression",
       kind = "sequential", high_color = "#00BCD4",
       opts = list(kind = "sequential", high_color = "#00BCD4", pt = 0.9)),
  list(col = "gex_IL1B",            pretty = "IL1B expression",
       kind = "sequential", high_color = "#FFB300",
       opts = list(kind = "sequential", high_color = "#FFB300", pt = 0.9)),
  list(col = "roi",                 pretty = "ROI (NFkB-high & Neu-high)",
       kind = "binary",     high_color = "#B2182B",
       opts = list(kind = "binary",     high_color = "#B2182B",
                   pt_bg = 0.40, pt_roi = 0.95))
)

cat("\n[atlas] generating ...\n")
for (cn in names(COHORTS)) {
  for (ft in FEATS) {
    safe_name <- gsub("[^A-Za-z0-9]+", "_", ft$col)
    pfx_clean <- gsub("-", "", cn)
    base <- file.path(OUT, sprintf("S10atlas_%s_%s", pfx_clean, safe_name))
    cat(sprintf("  >> %s / %s\n", cn, ft$col))
    make_atlas(cn, ft$col, ft$opts, mode = "hires",  ncol_grid = 4,
               panel_w = 1.6, panel_h = 1.6, out_stem = base)
    make_atlas(cn, ft$col, ft$opts, mode = "lowres", ncol_grid = 4,
               panel_w = 1.6, panel_h = 1.6, out_stem = paste0(base, "_lite"))
  }
}
cat("\nDONE.\n")
