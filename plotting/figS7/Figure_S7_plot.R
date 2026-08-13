# Supplementary Figure S7 — Neutrophil supporting validation
# 8 panels (A-H) aligned with fig4 / S6 strict style
# Outputs: fig5/supp/figS7{a..h}_*.{pdf,png}

required_cran <- c("ggplot2", "dplyr", "tidyr", "data.table",
                   "patchwork", "RColorBrewer", "scales",
                   "showtext", "ragg", "survival", "survminer")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

library(ggplot2); library(dplyr); library(tidyr); library(data.table)
library(patchwork); library(RColorBrewer); library(scales)
library(showtext); library(ragg); library(grid)
library(survival); library(survminer)

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
  font_add("Arial", regular = .arial_reg, bold = .arial_bold,
           italic = .arial_italic)
  showtext_auto(); showtext_opts(dpi = 300)
}, error = function(e) message("font_add failed: ", conditionMessage(e)))

.has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
.rasterise_pt <- function(p, dpi = 600) {
  if (.has_ggrastr) ggrastr::rasterise(p, dpi = dpi) else p
}

# ── UMAP corner-arrow helper ──
umap_arrow_axes <- function(data, x_col, y_col,
                            frac = 0.16, label_x = "UMAP 1", label_y = "UMAP 2",
                            text_size = 2.4, line_size = 0.45, arrow_mm = 1.6,
                            inset_frac = 0.02) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr); y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr);          y1 <- y0 + frac * diff(yr)
  arr <- grid::arrow(length = grid::unit(arrow_mm, "mm"), ends = "last", type = "closed")
  fam <- if (exists("my_font")) my_font else "sans"
  list(
    annotate("segment", x=x0, xend=x1, y=y0, yend=y0, arrow=arr, linewidth=line_size, color="black"),
    annotate("segment", x=x0, xend=x0, y=y0, yend=y1, arrow=arr, linewidth=line_size, color="black"),
    annotate("text", x=(x0+x1)/2, y=y0, label=label_x, vjust=2.4, size=text_size, family=fam),
    annotate("text", x=x0, y=(y0+y1)/2, label=label_y, angle=90, vjust=-1.4, size=text_size, family=fam),
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          panel.grid=element_blank())
  )
}

.fig_dir <- if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig5" else
  "${WORK_ROOT}/luad_figures/fig5"
setwd(.fig_dir)
dir.create("supp", showWarnings = FALSE)


theme_pub <- function(base_size = 10) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1))
}

neu_order <- c("Neu_Inflammatory", "Neu_IFN_response", "Neu_Angiogenic",
               "Neu_Metastatic", "Neu_ECM_remodeling", "Neu_OSM_priming",
               "Neu_OSM_low", "Neu_unclassified")
neu_colors <- c(
  "Neu_Inflammatory"   = "#E64B35",
  "Neu_IFN_response"   = "#F39B7F",
  "Neu_Angiogenic"     = "#00A087",
  "Neu_Metastatic"     = "#3C5488",
  "Neu_ECM_remodeling" = "#4DBBD5",
  "Neu_OSM_priming"    = "#7570B3",
  "Neu_OSM_low"        = "#91D1C2",
  "Neu_unclassified"   = "#B09C85"
)

tissue_order <- c("Normal_Lung", "Adjacent_Normal", "Normal_LN",
                  "Precancerous", "Primary_Tumor", "LN_Metastasis",
                  "Brain_Metastasis", "Distant_Metastasis", "Pleural_Effusion")

# Stratified downsample
strat_sample <- function(dt, group_col, target_total, seed = 42) {
  set.seed(seed)
  total_n <- nrow(dt)
  if (total_n <= target_total) return(dt[sample(.N)])
  out <- dt[, {
    target <- min(.N, as.integer(ceiling(as.numeric(target_total) * .N / total_n)) + 200L)
    .SD[sample(.N, target)]
  }, by = c(group_col)]
  out[sample(.N)]
}

stars_fn <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}


#  Load data
cat("\n== Loading S7 inputs ==\n")
meta <- as.data.table(fread("data/fig5a_umap_metadata.csv.gz"))
setnames(meta, c("UMAP1", "UMAP2"), c("UMAP_1", "UMAP_2"))
meta$neu_subtype <- factor(meta$neu_subtype, levels = neu_order)
meta <- meta[!is.na(neu_subtype)]
cat(sprintf("  Neu cells: %d, datasets: %d, tissues: %d\n",
            nrow(meta), uniqueN(meta$dataset), uniqueN(meta$tissue_type)))

mp1_cat   <- as.data.table(fread("~/luad/results/step29_mp1_validation/mp1_gene_categories.csv"))
tcga_pcor <- as.data.table(fread("~/luad/results/step29_mp1_validation/tcga_partial_correlations.csv"))
tcga_ct   <- as.data.table(fread("~/luad/results/step29_mp1_validation/tcga_celltype_correlations.csv"))
rename_map <- as.data.table(fread("~/luad/results/step25h_rename_map.csv"))


#  S7A — Neu UMAP × tissue_type
cat("-- figS7a: Neu UMAP × tissue --\n")
m_a <- meta[!is.na(tissue_type) & tissue_type %in% tissue_order,
            .(UMAP_1, UMAP_2, tissue_type)]
m_a$tissue_type <- factor(m_a$tissue_type, levels = tissue_order)
m_a <- strat_sample(m_a, "tissue_type", 8000L)

tissue_colors <- setNames(
  brewer.pal(9, "Set1"), tissue_order)

p7a <- ggplot(as.data.frame(m_a), aes(UMAP_1, UMAP_2, color = tissue_type)) +
  .rasterise_pt(geom_point(size = 0.7, stroke = 0, alpha = 0.5, shape = 16)) +
  scale_color_manual(values = tissue_colors, breaks = tissue_order, name = "Tissue") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
  coord_fixed(ratio = 1, clip = "off") + labs(x = NULL, y = NULL) +
  theme_pub(9) +
  theme(legend.position = "right",
        plot.background = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        plot.margin = margin(6, 6, 26, 26, "pt")) +
  umap_arrow_axes(as.data.frame(m_a), "UMAP_1", "UMAP_2",
                  inset_frac = -0.07, frac = 0.14,
                  arrow_mm = 1.2, line_size = 0.3)

ggsave("supp/figS7a_umap_tissue.pdf", p7a, width = 110, height = 80, units = "mm")
ggsave("supp/figS7a_umap_tissue.png", p7a, width = 110, height = 80, units = "mm",
       dpi = 300)
cat("  figS7a saved\n")


#  S7B — Neu UMAP × dataset
cat("-- figS7b: Neu UMAP × dataset --\n")
m_b <- meta[, .(UMAP_1, UMAP_2, dataset)]
m_b <- strat_sample(m_b, "dataset", 8000L)

ds_lvls <- sort(unique(m_b$dataset))
# Dark2 + Set1 give saturated colours that all pop on white background
.ds_pool <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))
.ds_pool <- .ds_pool[!.ds_pool %in% c("#FFFFFF", "#F0F0F0", "#FFFF33")]
ds_colors <- setNames(.ds_pool[seq_along(ds_lvls)], ds_lvls)

p7b <- ggplot(as.data.frame(m_b), aes(UMAP_1, UMAP_2, color = dataset)) +
  .rasterise_pt(geom_point(size = 0.7, stroke = 0, alpha = 0.5, shape = 16)) +
  scale_color_manual(values = ds_colors, name = "Dataset") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
  coord_fixed(ratio = 1, clip = "off") + labs(x = NULL, y = NULL) +
  theme_pub(9) +
  theme(legend.position = "right",
        plot.background = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        plot.margin = margin(6, 6, 26, 26, "pt")) +
  umap_arrow_axes(as.data.frame(m_b), "UMAP_1", "UMAP_2",
                  inset_frac = -0.07, frac = 0.14,
                  arrow_mm = 1.2, line_size = 0.3)

ggsave("supp/figS7b_umap_dataset.pdf", p7b, width = 110, height = 80, units = "mm")
ggsave("supp/figS7b_umap_dataset.png", p7b, width = 110, height = 80, units = "mm",
       dpi = 300)
cat("  figS7b saved\n")


#  S7C — scANVI uncertainty violin per neu_subtype
cat("-- figS7c: scANVI uncertainty ridgeline --\n")
m_c <- meta[!is.na(scanvi_uncertainty), .(neu_subtype, scanvi_uncertainty)]
m_c$neu_subtype <- factor(m_c$neu_subtype, levels = rev(neu_order))

# Replace violin with ridgeline (ggridges) — top-journal style for
# multi-group distributions; far cleaner than violin.
suppressPackageStartupMessages({
  if (requireNamespace("ggridges", quietly = TRUE)) library(ggridges)
})

if (requireNamespace("ggridges", quietly = TRUE)) {
  p7c <- ggplot(as.data.frame(m_c),
                aes(x = scanvi_uncertainty, y = neu_subtype,
                    fill = neu_subtype)) +
    ggridges::geom_density_ridges(
      scale = 1.6,
      rel_min_height = 0.01,
      linewidth = 0.3,
      color = "black",
      alpha = 0.9,
      quantile_lines = TRUE,
      quantiles = 0.5,
      vline_color = "white",
      vline_size = 0.3
    ) +
    scale_fill_manual(values = neu_colors, guide = "none") +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = "scANVI prediction uncertainty", y = NULL,
         title = "Per-subtype prediction confidence") +
    theme_pub(8) +
    theme(axis.text.y = element_text(size = 7, face = "italic"),
          axis.text.x = element_text(size = 7),
          plot.title  = element_text(size = 9, face = "bold"),
          panel.grid.major.x = element_line(color = "grey94", linewidth = 0.2))
} else {
  # fallback to violin (still styled)
  p7c <- ggplot(as.data.frame(m_c),
                aes(x = neu_subtype, y = scanvi_uncertainty, fill = neu_subtype)) +
    geom_violin(scale = "width", linewidth = 0.3, color = "black",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = neu_colors, guide = "none") +
    labs(x = NULL, y = "scANVI prediction uncertainty",
         title = "Per-subtype prediction confidence") +
    theme_pub(8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                     face = "italic"),
          axis.text.y = element_text(size = 7),
          plot.title  = element_text(size = 9, face = "bold"))
}

ggsave("supp/figS7c_scanvi_ridges.pdf", p7c, width = 130, height = 100, units = "mm")
ggsave("supp/figS7c_scanvi_ridges.png", p7c, width = 130, height = 100, units = "mm",
       dpi = 300)
# clean up old violin file if present
for (ext in c("pdf", "png"))
  if (file.exists(paste0("supp/figS7c_scanvi_violin.", ext)))
    file.remove(paste0("supp/figS7c_scanvi_violin.", ext))
cat("  figS7c saved (ridgeline)\n")


#  S7D — leiden_0_6 × neu_subtype confusion (cell counts)
cat("-- figS7d: confusion matrix --\n")
xt <- meta[, .N, by = .(leiden_0_6, neu_subtype)]
xt[, leiden_0_6 := factor(leiden_0_6,
                            levels = sort(unique(as.integer(leiden_0_6))))]
xt$neu_subtype <- factor(xt$neu_subtype, levels = neu_order)

# row-normalized fraction (for color), keep N for text
xt[, total := sum(N), by = leiden_0_6]
xt[, frac := N / total]

p7d <- ggplot(as.data.frame(xt),
              aes(x = neu_subtype, y = leiden_0_6, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = N), size = 2.0, family = "Arial") +
  scale_fill_gradient(name = "Row\nfraction",
                      low = "#F7F7F7", high = "#E64B35", limits = c(0, 1),
                      guide = guide_colorbar(barwidth = unit(3, "mm"),
                                             barheight = unit(15, "mm"),
                                             frame.colour = "black",
                                             frame.linewidth = 0.3,
                                             ticks.colour = "black")) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = "leiden_0.6 cluster",
       title = "Data-driven cluster vs. scANVI label") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 7,
                                   face = "italic"),
        axis.text.y = element_text(size = 7),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

# Square cells: width/height proportional to ncol/nrow at 7mm cell size
.s7d_nx <- length(unique(xt$neu_subtype))
.s7d_ny <- length(unique(xt$leiden_0_6))
.s7d_w  <- .s7d_nx * 7 + 60  # legend + label margin
.s7d_h  <- .s7d_ny * 7 + 30
ggsave("supp/figS7d_confusion.pdf", p7d, width = .s7d_w, height = .s7d_h,
       units = "mm")
ggsave("supp/figS7d_confusion.png", p7d, width = .s7d_w, height = .s7d_h,
       units = "mm", dpi = 300)
cat("  figS7d saved\n")


#  S7E — MP1 gene composition by category (bar + table)
cat("-- figS7e: MP1 gene categories --\n")
mp1_cat <- mp1_cat[!is.na(category) & category != "category"]
cat_counts <- mp1_cat[, .N, by = category][order(-N)]
cat_total <- sum(cat_counts$N)
cat_counts[, pct := round(100 * N / cat_total, 1)]
cat_counts$category <- factor(cat_counts$category,
                                levels = rev(cat_counts$category))

cat_palette <- c(
  "AP-1 / stress"         = "#E64B35",
  "Immune marker"         = "#3C5488",
  "Chemokine/Cytokine"    = "#00A087",
  "Heat shock"            = "#F39B7F",
  "Other"                 = "grey70"
)

p7e <- ggplot(as.data.frame(cat_counts),
              aes(x = category, y = N, fill = category)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", N, pct)),
            hjust = -0.1, size = 2.5, family = "Arial") +
  scale_fill_manual(values = cat_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(x = NULL, y = "Number of genes (top 50 MP1)",
       title = "MP1 gene composition") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, face = "bold"))

ggsave("supp/figS7e_mp1_categories.pdf", p7e, width = 130, height = 60, units = "mm")
ggsave("supp/figS7e_mp1_categories.png", p7e, width = 130, height = 60, units = "mm",
       dpi = 300)
cat("  figS7e saved\n")


#  S7F — TCGA partial correlation heatmap
cat("-- figS7f: TCGA partial correlation --\n")
tcga_pcor <- tcga_pcor[!is.na(rho)]
# keep only Neu_* on x (others may exist), MP3 on y
tcga_pcor <- tcga_pcor[grepl("^Neu_", x) & y == "MP3"]
control_order <- c("(none)", "PTPRC_expression", "Generic_Immune",
                   "Macro_general", "T_cell_core", "MP1")
control_labels <- c("(none)" = "Unadjusted",
                    "PTPRC_expression" = "| PTPRC",
                    "Generic_Immune"   = "| Generic_Immune",
                    "Macro_general"    = "| Macro_general",
                    "T_cell_core"      = "| T_cell_core",
                    "MP1"              = "| MP1")
neu_x_levels <- intersect(neu_order, unique(tcga_pcor$x))
tcga_pcor$x <- factor(tcga_pcor$x, levels = neu_x_levels)
tcga_pcor$control <- factor(tcga_pcor$control, levels = control_order)
tcga_pcor$rho_clamp <- pmax(pmin(tcga_pcor$rho, 0.6), -0.6)
tcga_pcor$star <- stars_fn(tcga_pcor$p)

p7f <- ggplot(as.data.frame(tcga_pcor),
              aes(x = control, y = x, fill = rho_clamp)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = 2.0, vjust = 1.4, family = "Arial") +
  geom_text(aes(label = star),
            size = 2.4, vjust = -0.6, fontface = "bold", family = "Arial") +
  scale_fill_gradient2(name = "Partial ρ",
                       low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                       midpoint = 0, limits = c(-0.6, 0.6),
                       oob = scales::squish,
                       guide = guide_colorbar(barwidth = unit(3, "mm"),
                                              barheight = unit(15, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.3,
                                              ticks.colour = "black")) +
  scale_x_discrete(position = "top", labels = control_labels) +
  coord_fixed(ratio = 1) +
  labs(x = "Control variable", y = NULL,
       title = "Neu × MP3 (TCGA, partial Spearman)") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 7),
        axis.text.y = element_text(size = 7, face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

# Square cells via coord_fixed; canvas auto-sized at 7mm/cell
.s7f_nx <- length(unique(tcga_pcor$control))
.s7f_ny <- length(unique(tcga_pcor$x))
.s7f_w  <- .s7f_nx * 8 + 60
.s7f_h  <- .s7f_ny * 8 + 30
ggsave("supp/figS7f_tcga_partial.pdf", p7f, width = .s7f_w, height = .s7f_h,
       units = "mm")
ggsave("supp/figS7f_tcga_partial.png", p7f, width = .s7f_w, height = .s7f_h,
       units = "mm", dpi = 300)
cat("  figS7f saved\n")


#  S7G — Cross-celltype × {MP1, MP3, EMT_Hallmark} bar
cat("-- figS7g: cross-celltype rho bar --\n")
tcga_ct <- tcga_ct[!is.na(rho) & target %in% c("MP1", "MP3", "EMT_Hallmark")]
sig_levels <- tcga_ct[target == "MP3"][order(-rho), signature]
tcga_ct$signature <- factor(tcga_ct$signature, levels = rev(sig_levels))
tcga_ct$target <- factor(tcga_ct$target, levels = c("MP1", "MP3", "EMT_Hallmark"))
tcga_ct$star <- stars_fn(tcga_ct$p)
target_colors <- c("MP1" = "#E64B35", "MP3" = "#00A087", "EMT_Hallmark" = "#3C5488")

p7g <- ggplot(as.data.frame(tcga_ct),
              aes(x = signature, y = rho, fill = target)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = star),
            position = position_dodge(width = 0.75),
            vjust = ifelse(tcga_ct$rho >= 0, -0.2, 1.2),
            size = 2.2, family = "Arial", fontface = "bold") +
  scale_fill_manual(values = target_colors, name = "Target") +
  coord_flip() +
  labs(x = NULL, y = "Spearman ρ (TCGA, n = 528)",
       title = "Cell-type signature × MP / EMT (TCGA)") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        axis.text.x = element_text(size = 7),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        plot.title = element_text(size = 9, face = "bold"))

ggsave("supp/figS7g_celltype_bar.pdf", p7g, width = 150, height = 110, units = "mm")
ggsave("supp/figS7g_celltype_bar.png", p7g, width = 150, height = 110, units = "mm",
       dpi = 300)
cat("  figS7g saved\n")


#  S7H — Three n.s. KM curves (Inflam, OSM, ECM × MP1)
cat("-- figS7h: 3 n.s. KM curves --\n")

build_km <- function(csv_path, sig_name, title) {
  d <- as.data.frame(fread(csv_path))
  fit <- survfit(Surv(OS_days, OS_status) ~ combined_group, data = d)
  lr <- survdiff(Surv(OS_days, OS_status) ~ combined_group, data = d)
  pval <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  # Match Fig 5I palette: grey (best/ref) → blue → salmon → red (worst)
  km_pal <- c("#7E7E7E", "#3C5488", "#F39B7F", "#E64B35")
  p <- ggsurvplot(
    fit, data = d, pval = TRUE, conf.int = FALSE,
    risk.table = FALSE,
    palette = km_pal,
    legend.title = "", legend.labs = sort(unique(d$combined_group)),
    xlab = "Overall survival (days)", ylab = "Survival probability",
    ggtheme = theme_pub(8) +
      theme(plot.title = element_text(size = 9, face = "bold")),
    title = title,
    legend = "right",
    font.family = "Arial",
    pval.size = 3, pval.coord = c(max(d$OS_days, na.rm = TRUE) * 0.62, 0.92)
  )
  p$plot
}

p_inflam <- build_km("data/fig5i_km_neuInflam_mp1.csv",
                      "Neu_Inflammatory",
                      "Neu_Inflammatory × MP1")
p_osm <- build_km("data/fig5i_km_neuOSM_mp1.csv",
                   "Neu_OSM_priming",
                   "Neu_OSM_priming × MP1")
p_ecm <- build_km("data/fig5i_km_neuECM_mp1.csv",
                   "Neu_ECM_remodeling",
                   "Neu_ECM_remodeling × MP1")

# Save each KM panel separately so curves aren't cramped.  Each gets a
# generous ~120 mm × 90 mm canvas — plenty of room for 4 curves + legend.
ggsave("supp/figS7h_km_inflam.pdf", p_inflam, width = 120, height = 90, units = "mm")
ggsave("supp/figS7h_km_inflam.png", p_inflam, width = 120, height = 90, units = "mm",
       dpi = 300)
ggsave("supp/figS7h_km_osm.pdf",    p_osm,    width = 120, height = 90, units = "mm")
ggsave("supp/figS7h_km_osm.png",    p_osm,    width = 120, height = 90, units = "mm",
       dpi = 300)
ggsave("supp/figS7h_km_ecm.pdf",    p_ecm,    width = 120, height = 90, units = "mm")
ggsave("supp/figS7h_km_ecm.png",    p_ecm,    width = 120, height = 90, units = "mm",
       dpi = 300)

# Combined with much more breathing room (was 250×80 → 360×130),
# AND optionally stack vertically if you'd rather use the 3 panels
# in a tall column.
p7h_horiz <- (p_inflam | p_osm | p_ecm) + plot_layout(ncol = 3)
ggsave("supp/figS7h_km_3panels.pdf", p7h_horiz,
       width = 360, height = 130, units = "mm")
ggsave("supp/figS7h_km_3panels.png", p7h_horiz,
       width = 360, height = 130, units = "mm", dpi = 300)
p7h_vert  <- p_inflam / p_osm / p_ecm
ggsave("supp/figS7h_km_3panels_vert.pdf", p7h_vert,
       width = 130, height = 270, units = "mm")
ggsave("supp/figS7h_km_3panels_vert.png", p7h_vert,
       width = 130, height = 270, units = "mm", dpi = 300)
cat("  figS7h saved (3 separate + horizontal + vertical combined)\n")


#  Summary
cat("\n==============================\n")
cat("Supp Fig S7 outputs:\n")
files <- list.files("supp", pattern = "^figS7.*\\.(pdf|png)$",
                    full.names = TRUE)
for (f in files) {
  sz <- round(file.info(f)$size / 1e6, 2)
  cat(sprintf("  %s (%.2f MB)\n", f, sz))
}
cat("==============================\n")
