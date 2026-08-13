# Fig 8 v2 (500-cell) shared R styling — mirrors Python fig8_style.py
# Source this at top of every R panel script.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(scales)
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

# ---- Arial setup (match other Fig scripts) ----
.find_arial <- function() {
  candidates <- c("~/.local/share/fonts/arial.ttf",
                  "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
                  "/mnt/c/Windows/Fonts/arial.ttf",
                  "/mnt/c/Windows/Fonts/Arial.ttf")
  for (p in candidates) if (file.exists(path.expand(p))) return(path.expand(p))
  return(NA_character_)
}
.arial <- .find_arial()
if (!is.na(.arial) && requireNamespace("sysfonts", quietly = TRUE) &&
    requireNamespace("showtext", quietly = TRUE)) {
  bold_p <- sub("arial.ttf$", "arialbd.ttf", .arial, ignore.case = TRUE)
  it_p   <- sub("arial.ttf$", "ariali.ttf",  .arial, ignore.case = TRUE)
  sysfonts::font_add("Arial", regular = .arial,
                     bold = if (file.exists(bold_p)) bold_p else .arial,
                     italic = if (file.exists(it_p)) it_p else .arial)
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 300)
  FAM <- "Arial"
} else FAM <- "sans"

# ---- theme_pub (no titles by default — set via add_subtitle if needed) ----
theme_pub <- function(base_size = 8) {
  theme_classic(base_family = FAM, base_size = base_size) +
    theme(
      axis.text   = element_text(color = "black", size = base_size - 1),
      axis.title  = element_text(size = base_size),
      axis.line   = element_line(linewidth = 0.4, color = "black"),
      axis.ticks  = element_line(linewidth = 0.4, color = "black"),
      plot.title  = element_blank(),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25",
                                   margin = margin(b = 2)),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 1, color = "black", face = "plain"),
      panel.grid = element_blank()
    )
}

# ---- palette ----
COL <- list(
  Macro_C1QC="#4DBBD5", Macro_FCN1="#E64B35", Macro_FOLR2="#00A087",
  Macro_MARCO="#3C5488", Macro_SPP1="#F39B7F",
  Macro_general="#8491B4", Macro_prolif="#91D1C2",
  MP1="#E64B35", MP2="#4DBBD5", MP3="#00A087", MP4="#3C5488",
  Neu_Inflammatory="#E64B35", Neu_Angiogenic="#F39B7F",
  Neu_Metastatic="#3C5488", Neu_ECM_remodeling="#4DBBD5",
  Neu_OSM_priming="#00A087", Neu_OSM_low="#8491B4",
  Neu_IFN_response="#91D1C2", Neu_unclassified="#D9D9D9",
  tumor="#E64B35",   normal="#4DBBD5",
  high="#E64B35",    low="#4DBBD5",
  NR="#E64B35",      R="#4DBBD5",
  MPR="#4DBBD5",     NMPR="#E64B35",
  ROI="#E64B35",     nonROI="#4DBBD5",
  LUAD="#E64B35",    other="#8491B4",
  hr_up="#E64B35",   hr_down="#4DBBD5",
  ref_red="#cb181d", grid="#cccccc",
  venn_macro="#4DBBD5", venn_mal="#E64B35", venn_neu="#00A087"
)

# Sequential ST gene-expression colormap: blue → red, no white middle
ST_RAMP <- c("#4DBBD5", "#7C95B8", "#A07499", "#C7547A", "#E64B35")
scale_st_color <- function(...) scale_color_gradientn(colors = ST_RAMP, ...)
scale_st_fill  <- function(...) scale_fill_gradientn(colors  = ST_RAMP, ...)

# Significance stars
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# Threshold-text p-value formatter (avoid scientific notation in figures).
# Returns e.g. "p < 0.001 ***", "p = 0.034 *", "p = 0.42 ns".
fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi) || !is.finite(pi)) return("p = NA")
    if (pi < 1e-4) return("p < 0.0001 ****")
    if (pi < 1e-3) return("p < 0.001 ***")
    if (pi < 1e-2) return("p < 0.01 **")
    if (pi < 0.05) return(sprintf("p = %.3f *", pi))
    return(sprintf("p = %.3f ns", pi))
  }, FUN.VALUE = character(1))
}

# Save panel: writes both .pdf and .png at 300 dpi to OUT/<stem>.{pdf,png}
save_panel <- function(p, stem, w, h) {
  for (ext in c("pdf", "png")) {
    fp <- paste0(stem, ".", ext)
    ggsave(fp, p, width = w, height = h, units = "in",
           dpi = 300, device = ext, bg = "white")
  }
  invisible(NULL)
}

# Italic subtitle helper used by Python add_subtitle (small grey25)
add_sub <- function(p, txt) p + labs(subtitle = txt)
