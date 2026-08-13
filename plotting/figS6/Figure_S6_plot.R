# Supplementary Figure S6 — Myeloid additional dot plots + GO
# Strict alignment with Fig 1B / Fig 4H templates
# Outputs: fig4/supp/figS6{a..d}_*.{pdf,png}

required_cran <- c("ggplot2", "dplyr", "tidyr", "data.table",
                   "patchwork", "RColorBrewer", "scales",
                   "showtext", "ragg")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

library(ggplot2); library(dplyr); library(tidyr); library(data.table)
library(patchwork); library(RColorBrewer); library(scales)
library(showtext); library(ragg); library(grid)

# ── Arial (real, from ~/.local/share/fonts on Linux or C:/Windows/Fonts) ──
.find_arial <- function() {
  cands <- c("arial.ttf", "C:/Windows/Fonts/arial.ttf",
             "~/.local/share/fonts/arial.ttf",
             "/mnt/c/Windows/Fonts/arial.ttf")
  for (p in cands) {
    pp <- path.expand(p)
    if (file.exists(pp) || p == "arial.ttf") return(p)
  }
  return("sans")
}
.arial_reg <- .find_arial()
.arial_dir <- dirname(path.expand(.arial_reg))
.arial_bold <- if (file.exists(file.path(.arial_dir, "arialbd.ttf")))
  file.path(.arial_dir, "arialbd.ttf") else .arial_reg
.arial_italic <- if (file.exists(file.path(.arial_dir, "ariali.ttf")))
  file.path(.arial_dir, "ariali.ttf") else .arial_reg
tryCatch({
  font_add("Arial", regular = .arial_reg, bold = .arial_bold, italic = .arial_italic)
  showtext_auto(); showtext_opts(dpi = 300)
}, error = function(e) message("Arial font registration failed: ", conditionMessage(e)))

.fig_dir <- if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig4" else
  "${WORK_ROOT}/luad_figures/fig4"
setwd(.fig_dir)
dir.create("supp", showWarnings = FALSE)


# ── Theme + color ──
theme_pub <- function(base_size = 10) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1))
}

# Major + Macro ordering / colors (same as fig4)
major_order <- c("Macrophage", "Mono_nonclassical", "Neutrophil",
                 "cDC1", "cDC2", "cDC_LAMP3", "pDC")
macro_order <- c("Macro_C1QC", "Macro_FCN1", "Macro_FOLR2", "Macro_MARCO",
                 "Macro_SPP1", "Macro_general", "Macro_prolif")

# ── Helper: per-gene z-score (long form) ──
add_gene_zscore <- function(dt, value_col = "mean_log1p", gene_col = "gene") {
  dt <- as.data.frame(dt)
  dt$z_expr <- ave(dt[[value_col]], dt[[gene_col]],
                   FUN = function(x) {
                     s <- sd(x); if (is.na(s) || s == 0) return(rep(0, length(x)))
                     (x - mean(x)) / s
                   })
  dt
}

build_dotplot <- function(d_long, group_levels, gene_levels, title = NULL,
                          legend_show = TRUE) {
  d <- as.data.frame(d_long)
  d$grp <- factor(d$grp, levels = rev(group_levels))   # rev: bottom-to-top y
  d$gene <- factor(d$gene, levels = gene_levels)
  d$pct100 <- d$pct_expressing * 100
  ggplot(d, aes(x = gene, y = grp)) +
    geom_point(aes(size = pct100, fill = z_expr),
               shape = 21, color = "black", stroke = 0.25) +
    scale_size_continuous(name = "Fraction of cells\nin group (%)",
                          range = c(0.3, 5.5), breaks = c(25, 50, 75),
                          limits = c(0, 100)) +
    scale_fill_gradient2(name = "Mean expression\n(z-score)",
                         low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                         midpoint = 0, limits = c(-2, 2),
                         oob = scales::squish,
                         guide = guide_colorbar(barwidth = unit(3.5, "mm"),
                                                barheight = unit(20, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.3,
                                                ticks.colour = "black")) +
    labs(x = NULL, y = NULL, title = title) +
    theme_pub(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 7, face = "italic"),
          axis.text.y = element_text(size = 8),
          legend.position = if (legend_show) "right" else "none",
          legend.box = "vertical",
          legend.spacing.y = unit(4, "mm"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


#  S6A — 7 major-type dotplot, ~18 markers
cat("\n-- figS6a: major-type dotplot --\n")

dot_sub <- as.data.table(fread("myeloid_dotplot_markers_refined.csv"))
m12 <- as.data.table(fread("myeloid_m1m2_scores_refined.csv"))
setnames(m12, "myeloid_subtype_refined", "subtype")
m12 <- m12[, .(subtype, n_cells = M1_score_count)]

# Map subtype -> major type
maj_map <- c(
  "Macro_C1QC"="Macrophage", "Macro_FCN1"="Macrophage",
  "Macro_FOLR2"="Macrophage", "Macro_MARCO"="Macrophage",
  "Macro_SPP1"="Macrophage", "Macro_general"="Macrophage",
  "Macro_prolif"="Macrophage",
  "Mono_nonclassical"="Mono_nonclassical",
  "Neutrophil"="Neutrophil",
  "cDC1"="cDC1", "cDC2"="cDC2", "cDC_LAMP3"="cDC_LAMP3", "pDC"="pDC"
)
dot_sub[, major := maj_map[subtype]]
dot_sub <- merge(dot_sub, m12, by = "subtype", all.x = TRUE)
dot_sub[is.na(n_cells), n_cells := 1]

# Cell-count weighted aggregation per (major, gene)
dot_major <- dot_sub[, .(
  mean_log1p = sum(mean_log1p * n_cells) / sum(n_cells),
  pct_expressing = sum(pct_expressing * n_cells) / sum(n_cells)
), by = .(major, gene)]
setnames(dot_major, "major", "grp")

# Major-type marker panel — pan-myeloid → Mac → Mono → Neut → DC
markers_A <- c(
  "LYZ", "CD68", "CSF1R",                    # pan-myeloid
  "APOE", "C1QA", "C1QC", "MARCO",            # Macrophage
  "VCAN", "CD14", "S100A8", "S100A9",         # Mono / Neut shared
  "FCGR3A", "CDKN1C",                          # Mono_NC specific
  "CSF3R",                                      # Neutrophil
  "BATF3", "IRF8",                              # cDC1
  "LAMP3", "FSCN1",                             # cDC_LAMP3
  "TCF4", "IL3RA"                               # pDC
)
markers_A <- intersect(markers_A, unique(dot_major$gene))
cat(sprintf("  S6A markers used: %d/%d\n", length(markers_A), 20))

dot_A <- dot_major[gene %in% markers_A & grp %in% major_order]
dot_A <- add_gene_zscore(dot_A, "mean_log1p", "gene")

p6a <- build_dotplot(dot_A, group_levels = major_order, gene_levels = markers_A,
                     title = "Myeloid major-type markers")
ggsave("supp/figS6a_dotplot.pdf", p6a, width = 200, height = 80, units = "mm")
ggsave("supp/figS6a_dotplot.png", p6a, width = 200, height = 80, units = "mm",
       dpi = 300)
cat("  figS6a saved\n")


#  S6B — 7 Macro subset dotplot, ~22 markers
cat("-- figS6b: Macro subset dotplot --\n")

# Use per-subtype dotplot CSV plus extras from panel_GM
extra <- as.data.table(fread("panel_GM_subset_markers.csv"))
extra <- extra[, .(subtype, gene, mean_log1p, pct_expressing)]
dot_sub2 <- rbind(dot_sub[, .(subtype, gene, mean_log1p, pct_expressing)],
                   extra)
# Deduplicate (keep first occurrence per subtype × gene)
dot_sub2 <- dot_sub2[!duplicated(dot_sub2[, .(subtype, gene)])]

markers_B <- c(
  # C1QC
  "C1QA", "C1QB", "C1QC", "AXL", "GPR34",
  # FCN1
  "VCAN", "S100A8", "S100A9", "IL1B", "EREG", "FCGR3A",
  # FOLR2
  "FOLR2", "MRC1", "CD163", "MAF", "STAB1",
  # MARCO
  "MARCO", "MCEMP1", "PPARG", "VSIG4", "MSR1",
  # SPP1
  "SPP1", "TREM2", "MMP9", "VEGFA", "ADAM8", "ENO1", "P4HA1", "LDHA",
  # general
  "CD68", "CSF1R", "APOE",
  # prolif
  "TOP2A", "STMN1", "CDK1", "CCNB1"
)
markers_B <- intersect(markers_B, unique(dot_sub2$gene))
cat(sprintf("  S6B markers used: %d\n", length(markers_B)))

dot_B <- dot_sub2[gene %in% markers_B & subtype %in% macro_order]
setnames(dot_B, "subtype", "grp")
dot_B <- add_gene_zscore(dot_B, "mean_log1p", "gene")

p6b <- build_dotplot(dot_B, group_levels = macro_order, gene_levels = markers_B,
                     title = "Macrophage subset markers")
ggsave("supp/figS6b_macro_dotplot.pdf", p6b,
       width = 240, height = 75, units = "mm")
ggsave("supp/figS6b_macro_dotplot.png", p6b,
       width = 240, height = 75, units = "mm", dpi = 300)
cat("  figS6b saved\n")


#  S6C / S6D — GO BP bar plots (Macro_general, Macro_prolif)
cat("-- figS6c-d: GO bar plots --\n")

go <- as.data.table(fread("myeloid_go_enrichment.csv"))
go <- go[gene_set == "GO_Biological_Process_2023"]
go[, term_short := gsub(" \\(GO:[0-9]+\\)$", "", term)]
# Wrap long terms onto two lines instead of truncating
if (requireNamespace("stringr", quietly = TRUE)) {
  go[, term_short := stringr::str_wrap(term_short, width = 38)]
}

build_go_bar <- function(target_subtype, title, n_top = 10) {
  g <- go[subtype == target_subtype][order(-combined_score)][1:n_top]
  g <- g[!is.na(term_short)]
  g[, term_short := factor(term_short, levels = rev(unique(term_short)))]
  g$mlog_p <- pmin(-log10(g$adj_p_value), 12)
  ggplot(g, aes(x = term_short, y = mlog_p, fill = mlog_p)) +
    geom_col(width = 0.7) +
    # Lock the fill range so a SINGLE collected legend is one-to-one with
    # every panel.  Without this, each panel uses its own range and the
    # shared colour bar misrepresents most panels.
    scale_fill_gradient(name = "-log10(adj.P)",
                        low = "#FEE0D2", high = "#E64B35",
                        limits = c(0, 12), oob = scales::squish,
                        breaks = c(0, 4, 8, 12),
                        guide = guide_colorbar(barwidth = unit(3, "mm"),
                                               barheight = unit(20, "mm"),
                                               frame.colour = "black",
                                               frame.linewidth = 0.3,
                                               ticks.colour = "black")) +
    coord_flip() +
    labs(x = NULL, y = "-log10(adj.P)", title = title) +
    theme_pub(base_size = 8) +
    theme(axis.text.y = element_text(size = 6.5, lineheight = 0.85),
          axis.text.x = element_text(size = 7),
          legend.position = "right",
          plot.title = element_text(size = 9, face = "bold"),
          plot.margin = margin(4, 8, 4, 4, "pt"))
}

# All seven macrophage subtypes laid out as 4 on top and 3 on bottom.
# Single shared colour bar (locked to 0-12 -log10 adj.P) sits to the right.
go_subtypes <- c("Macro_SPP1", "Macro_C1QC", "Macro_FCN1", "Macro_FOLR2",
                 "Macro_MARCO", "Macro_general", "Macro_prolif")
go_panels <- lapply(go_subtypes, function(s)
  build_go_bar(s, sprintf("%s - top 10 GO BP", s), n_top = 10))

top_row    <- patchwork::wrap_plots(go_panels[1:4], nrow = 1)
bottom_row <- patchwork::wrap_plots(go_panels[5:7], nrow = 1)
p6cd <- (top_row / bottom_row) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "right")
ggsave("supp/figS6c_macro_go_grid.pdf", p6cd,
       width = 380, height = 170, units = "mm")
ggsave("supp/figS6c_macro_go_grid.png", p6cd,
       width = 380, height = 170, units = "mm", dpi = 300)
cat("  figS6c-i grid saved (4-on-top + 3-on-bottom, shared legend)\n")

# Keep the old single-panel exports as legacy (not referenced in main supp PDF)
p6c_legacy <- build_go_bar("Macro_general", "Macro_general - top 10 GO BP")
ggsave("supp/figS6_legacy_general_go.pdf", p6c_legacy,
       width = 195, height = 95, units = "mm")
p6d_legacy <- build_go_bar("Macro_prolif", "Macro_prolif - top 10 GO BP")
ggsave("supp/figS6_legacy_prolif_go.pdf", p6d_legacy,
       width = 195, height = 95, units = "mm")


#  Combined 2x2 — FigureS6.pdf with a/b/c/d tag labels
if (requireNamespace("patchwork", quietly = TRUE)) {
  .fam <- if (exists("my_font")) my_font else "sans"
  abrow <- (p6a / p6b)
  combo <- (abrow / p6cd) +
    patchwork::plot_layout(heights = c(1.0, 1.0, 1.8)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 11, face = "bold",
                                  family = .fam))
  ggsave("supp/FigureS6.pdf", combo,
         width = 380, height = 380, units = "mm")
  ggsave("supp/FigureS6.png", combo,
         width = 380, height = 380, units = "mm", dpi = 300)
  cat("  Combined FigureS6 saved (a, b, c-i grid 4+3)\n")
}


#  Summary
cat("\n==============================\n")
cat("Supp Fig S6 outputs:\n")
files <- list.files("supp", pattern = "^figS6.*\\.(pdf|png)$",
                    full.names = TRUE)
for (f in files) {
  sz <- round(file.info(f)$size / 1e6, 2)
  cat(sprintf("  %s (%.2f MB)\n", f, sz))
}
cat("==============================\n")
