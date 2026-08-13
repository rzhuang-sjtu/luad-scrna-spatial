# Figure 4: LUAD Myeloid 14-panel — strict alignment with fig1/2/3
# Panels A-N. Output: panels/fig4{a..n}_*.{pdf,png}

required_cran <- c("ggplot2", "dplyr", "tidyr", "data.table", "R.utils",
                   "patchwork", "RColorBrewer", "scales", "ggnewscale",
                   "ggrepel", "showtext", "ragg")
required_bioc <- c("ComplexHeatmap", "circlize")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg)
  }
}

library(ggplot2); library(dplyr); library(tidyr); library(data.table)
library(R.utils); library(patchwork); library(RColorBrewer); library(scales)
library(ggnewscale); library(ggrepel); library(showtext); library(ragg)
library(ComplexHeatmap); library(circlize); library(grid)

# ── Font: real Arial (Windows / WSL with copied arial.ttf) ──
.find_arial <- function() {
  cands <- c(
    "arial.ttf",                                              # Windows R search
    "C:/Windows/Fonts/arial.ttf",                             # explicit Windows
    "~/.local/share/fonts/arial.ttf",                         # WSL user font
    "/mnt/c/Windows/Fonts/arial.ttf",                         # WSL bridged
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",  # metric fallback
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
  )
  for (p in cands) {
    pp <- path.expand(p)
    if (file.exists(pp) || p == "arial.ttf") return(p)
  }
  return("sans")
}
.arial_reg <- .find_arial()
.arial_dir <- dirname(path.expand(.arial_reg))
.arial_bold <- if (file.exists(file.path(.arial_dir, "arialbd.ttf")))
  file.path(.arial_dir, "arialbd.ttf") else
  if (file.exists(file.path(.arial_dir, "LiberationSans-Bold.ttf")))
    file.path(.arial_dir, "LiberationSans-Bold.ttf") else .arial_reg
.arial_italic <- if (file.exists(file.path(.arial_dir, "ariali.ttf")))
  file.path(.arial_dir, "ariali.ttf") else
  if (file.exists(file.path(.arial_dir, "LiberationSans-Italic.ttf")))
    file.path(.arial_dir, "LiberationSans-Italic.ttf") else .arial_reg

.font_ok <- tryCatch({
  font_add("Arial", regular = .arial_reg, bold = .arial_bold,
           italic = .arial_italic)
  showtext_auto(); showtext_opts(dpi = 300)
  TRUE
}, error = function(e) FALSE)

# Register Arial as the family the grid graphics device knows ("sans" alias)
# Same path also gets used for ComplexHeatmap by aliasing my_font="Arial"
my_font <- if (.font_ok) "Arial" else "sans"
cat("Arial font registered:", .font_ok, "| using:", my_font, "\n")
cat("  regular:", .arial_reg, "\n")

.has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
.rasterise_pt <- function(p, dpi = 600) {
  if (.has_ggrastr) ggrastr::rasterise(p, dpi = dpi) else p
}
cat("ggrastr available:", .has_ggrastr, "\n")

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
  "${WORK_ROOT}/luad_figures/fig4" else
  "${WORK_ROOT}/luad_figures/fig4"
setwd(.fig_dir)
dir.create("panels", showWarnings = FALSE)


#  Theme + colors
theme_pub <- function(base_size = 10) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1))
}

# Major type 7 colors (per user spec)
major_order <- c("Macrophage", "Mono_nonclassical", "Neutrophil",
                 "cDC1", "cDC2", "cDC_LAMP3", "pDC")
major_colors <- c(
  "Macrophage"        = "#3C5488",
  "Mono_nonclassical" = "#F39B7F",
  "Neutrophil"        = "#E64B35",
  "cDC1"              = "#4DBBD5",
  "cDC2"              = "#00A087",
  "cDC_LAMP3"         = "#8491B4",
  "pDC"               = "#91D1C2"
)

# Macrophage subset 7 colors (per user spec)
macro_order <- c("Macro_C1QC", "Macro_FCN1", "Macro_FOLR2", "Macro_MARCO",
                 "Macro_SPP1", "Macro_general", "Macro_prolif")
macro_colors <- c(
  "Macro_C1QC"    = "#4DBBD5",
  "Macro_FCN1"    = "#E64B35",
  "Macro_FOLR2"   = "#00A087",
  "Macro_MARCO"   = "#3C5488",
  "Macro_SPP1"    = "#F39B7F",
  "Macro_general" = "#8491B4",
  "Macro_prolif"  = "#91D1C2"
)

# 13 subtypes (canonical order for dotplots / heatmaps)
sub_order <- c(
  "Macro_C1QC", "Macro_FCN1", "Macro_FOLR2", "Macro_MARCO", "Macro_SPP1",
  "Macro_general", "Macro_prolif",
  "Mono_nonclassical", "Neutrophil",
  "cDC1", "cDC2", "cDC_LAMP3", "pDC"
)
sub_colors <- c(macro_colors,
  "Mono_nonclassical" = "#F39B7F",
  "Neutrophil"        = "#F0A23B",
  "cDC1"              = "#8491B4",
  "cDC2"              = "#7570B3",
  "cDC_LAMP3"         = "#984EA3",
  "pDC"               = "#B09C85"
)

mp_colors <- c("MP1" = "#E64B35", "MP2" = "#4DBBD5",
               "MP3" = "#00A087", "MP4" = "#3C5488")

tissue_order <- c("Normal_Lung", "Adjacent_Normal", "Normal_LN",
                  "Precancerous", "Primary_Tumor", "LN_Metastasis",
                  "Brain_Metastasis", "Distant_Metastasis", "Pleural_Effusion")
tissue_labels <- c(
  "Normal_Lung" = "Normal", "Adjacent_Normal" = "Adjacent",
  "Normal_LN" = "Normal LN", "Precancerous" = "Pre-cancerous",
  "Primary_Tumor" = "Tumor", "LN_Metastasis" = "LN Met",
  "Brain_Metastasis" = "Brain Met", "Distant_Metastasis" = "Distant Met",
  "Pleural_Effusion" = "Pleural Effusion"
)

# Helper: stratified downsample
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

# Helper: per-row z-score for a wide matrix
row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

# Helper: long-form per-gene z-score
add_gene_zscore <- function(dt, value_col = "mean_log1p", gene_col = "gene") {
  dt <- as.data.frame(dt)
  dt$z_expr <- ave(dt[[value_col]], dt[[gene_col]],
                   FUN = function(x) {
                     s <- sd(x); if (is.na(s) || s == 0) return(rep(0, length(x)))
                     (x - mean(x)) / s
                   })
  dt
}

# Helper: stars
stars_fn <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}


#  LOAD DATA
cat("\n== Loading data ==\n")

meta_full <- as.data.frame(fread("panel_major_type_metadata.csv.gz"))
# Rename Mono_NC -> Mono_nonclassical (per spec)
meta_full$myeloid_major_type[meta_full$myeloid_major_type == "Mono_NC"] <- "Mono_nonclassical"
meta_full$major <- factor(meta_full$myeloid_major_type, levels = major_order)
meta_full$sub   <- factor(meta_full$myeloid_subtype_refined, levels = sub_order)
names(meta_full)[names(meta_full) == "UMAP1"] <- "UMAP_1"
names(meta_full)[names(meta_full) == "UMAP2"] <- "UMAP_2"
meta_full <- meta_full[!is.na(meta_full$major) & !is.na(meta_full$sub), ]
cat(sprintf("  myeloid: %d cells (%d major / %d sub)\n",
            nrow(meta_full), nlevels(meta_full$major), nlevels(meta_full$sub)))

dot_refined <- as.data.frame(fread("myeloid_dotplot_markers_refined.csv"))
prop_tissue <- as.data.frame(fread("myeloid_proportion_by_tissue_refined.csv"))
go_refined  <- as.data.frame(fread("myeloid_go_enrichment.csv"))
m12_refined <- as.data.frame(fread("myeloid_m1m2_scores_refined.csv"))
mp_refined  <- as.data.frame(fread("myeloid_mp3_association_refined.csv"))
ant_genes   <- as.data.frame(fread("panel_F_antitumor_genes.csv"))
sub_markers <- as.data.frame(fread("panel_GM_subset_markers.csv"))
gsea_n      <- as.data.frame(fread("panel_N_spp1_vs_c1qc_gsea.csv"))


#  Fig 4A — Major-type UMAP (Fig1-style: clean, shape=16, larger pts)
cat("\n-- fig4a: major UMAP --\n")
m4a <- as.data.table(meta_full)[, .(UMAP_1, UMAP_2, major)]
# rare types on top, then stratified downsample
set.seed(42)
m4a <- m4a[sample(.N)]
m4a_freq <- m4a[, .N, by = major][order(-N)]
m4a$major <- factor(m4a$major, levels = m4a_freq$major)
m4a <- m4a[order(major, decreasing = TRUE)]
m4a$major <- factor(m4a$major, levels = major_order)
m4a <- strat_sample(m4a, "major", 80000L)

p4a <- ggplot(as.data.frame(m4a), aes(x = UMAP_1, y = UMAP_2, color = major)) +
  .rasterise_pt(geom_point(size = 0.45, alpha = 0.55, stroke = 0, shape = 16),
                dpi = 600) +
  scale_color_manual(values = major_colors, breaks = major_order,
                     name = "celltype",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  umap_arrow_axes(as.data.frame(m4a), "UMAP_1", "UMAP_2")

ggsave("panels/fig4a_umap_major.pdf", p4a, width = 100, height = 85, units = "mm")
ggsave("panels/fig4a_umap_major.png", p4a, width = 100, height = 85, units = "mm",
       dpi = 300)
cat("  fig4a saved\n")


#  Fig 4B — Major-type tissue stacked bar
cat("-- fig4b: tissue major bar --\n")
prop_major <- as.data.table(meta_full)[, .N, by = .(tissue_type, major)]
prop_major[, pct := 100 * N / sum(N), by = tissue_type]
prop_major <- prop_major[tissue_type %in% tissue_order]
prop_major$tissue_type <- factor(prop_major$tissue_type, levels = tissue_order)
prop_major$major <- factor(prop_major$major, levels = rev(major_order))
prop_major <- as.data.frame(prop_major)

p4b <- ggplot(prop_major, aes(x = tissue_type, y = pct, fill = major)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = major_colors, breaks = major_order, name = NULL) +
  scale_x_discrete(labels = tissue_labels) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 100)) +
  labs(x = NULL, y = "Cell percent ratio") +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1,
                                   face = "italic"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        plot.margin = margin(4, 4, 4, 4, "pt"))

ggsave("panels/fig4b_tissue_major.pdf", p4b, width = 110, height = 65, units = "mm")
ggsave("panels/fig4b_tissue_major.png", p4b, width = 110, height = 65, units = "mm",
       dpi = 300)
cat("  fig4b saved\n")


#  Fig 4C — Macrophage subset UMAP (104k Macro cells)
cat("-- fig4c: Macro UMAP --\n")
mac <- as.data.table(meta_full)[major == "Macrophage" &
                                  myeloid_subtype_refined %in% macro_order]
mac[, sub_macro := factor(myeloid_subtype_refined, levels = macro_order)]

set.seed(42)
mac4c <- mac[sample(.N), .(UMAP_1, UMAP_2, sub_macro)]
mac_freq <- mac4c[, .N, by = sub_macro][order(-N)]
mac4c$sub_macro <- factor(mac4c$sub_macro, levels = mac_freq$sub_macro)
mac4c <- mac4c[order(sub_macro, decreasing = TRUE)]
mac4c$sub_macro <- factor(mac4c$sub_macro, levels = macro_order)
mac4c <- strat_sample(mac4c, "sub_macro", 70000L)

p4c <- ggplot(as.data.frame(mac4c), aes(x = UMAP_1, y = UMAP_2, color = sub_macro)) +
  .rasterise_pt(geom_point(size = 0.45, alpha = 0.55, stroke = 0, shape = 16),
                dpi = 600) +
  scale_color_manual(values = macro_colors, breaks = macro_order,
                     name = "Macrophage subset",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        plot.background = element_blank(),
        panel.background = element_blank()) +
  umap_arrow_axes(as.data.frame(mac4c), "UMAP_1", "UMAP_2")

ggsave("panels/fig4c_umap_macro.pdf", p4c, width = 100, height = 85, units = "mm")
ggsave("panels/fig4c_umap_macro.png", p4c, width = 100, height = 85, units = "mm",
       dpi = 300)
cat("  fig4c saved\n")


#  Fig 4D — Macrophage subset tissue stacked bar
cat("-- fig4d: tissue Macro bar --\n")
prop_mac <- mac[, .N, by = .(tissue_type, sub_macro)]
prop_mac[, pct := 100 * N / sum(N), by = tissue_type]
prop_mac <- prop_mac[tissue_type %in% tissue_order]
prop_mac$tissue_type <- factor(prop_mac$tissue_type, levels = tissue_order)
prop_mac$sub_macro <- factor(prop_mac$sub_macro, levels = rev(macro_order))
prop_mac <- as.data.frame(prop_mac)

p4d <- ggplot(prop_mac, aes(x = tissue_type, y = pct, fill = sub_macro)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = macro_colors, breaks = macro_order, name = NULL) +
  scale_x_discrete(labels = tissue_labels) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 100)) +
  labs(x = NULL, y = "Fraction within Macrophage") +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1,
                                   face = "italic"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        plot.margin = margin(4, 4, 4, 4, "pt"))

ggsave("panels/fig4d_tissue_macro.pdf", p4d, width = 110, height = 65, units = "mm")
ggsave("panels/fig4d_tissue_macro.png", p4d, width = 110, height = 65, units = "mm",
       dpi = 300)
cat("  fig4d saved\n")


#  Fig 4E — GO BP heatmap (Macro 7 subtypes × 2 terms each)
#           Single-direction warm palette, subtype row annotation
cat("-- fig4e: GO heatmap (Macro-only, 14 rows) --\n")
go_bp <- as.data.table(go_refined)[gene_set == "GO_Biological_Process_2023"]
go_bp[, term_short := gsub(" \\(GO:[0-9]+\\)$", "", term)]
go_bp[, term_short := ifelse(nchar(term_short) > 45,
                              paste0(substr(term_short, 1, 42), "..."), term_short)]

# Restrict to 7 Macrophage subtypes only — 4F shows the same set
go_bp_mac <- go_bp[subtype %in% macro_order]
# Top 2 per subtype, each term claimed by first occurring subtype
go_top_e <- go_bp_mac[order(-combined_score), head(.SD, 2), by = subtype]
term_origin_dt <- go_top_e[!duplicated(term_short),
                           .(term_short, origin = subtype)]
term_origin_dt$origin <- factor(term_origin_dt$origin, levels = macro_order)
term_origin_dt <- term_origin_dt[order(origin)]
keep_terms <- term_origin_dt$term_short

go_long <- go_bp_mac[term_short %in% keep_terms]
mat_e_dt <- dcast(go_long, term_short ~ subtype,
                   value.var = "adj_p_value",
                   fun.aggregate = function(x) min(x, na.rm = TRUE),
                   fill = 1)
sub_av_e <- intersect(macro_order, colnames(mat_e_dt))
mat_e <- as.matrix(mat_e_dt[, ..sub_av_e])
rownames(mat_e) <- mat_e_dt$term_short
mat_e[!is.finite(mat_e)] <- 1
mat_e_mlog <- pmin(-log10(mat_e), 8)

# Re-order rows by origin macro subtype
mat_e_mlog <- mat_e_mlog[keep_terms, , drop = FALSE]
row_origin <- term_origin_dt$origin

# Single-direction warm palette: white = non-sig, deep red = strong sig
col_fun_e <- colorRamp2(c(0, 1.3, 4, 8),
                         c("#F7F7F7", "#FEE0D2", "#FC9272", "#E64B35"))

# Both annotations share the same macro_colors mapping. We only need one
# legend (drawn from the row annotation) labelled "Macro subtype" — it
# explains both the top column strip and the left row strip.
ha_col_e <- HeatmapAnnotation(
  Subtype = sub_av_e,
  col = list(Subtype = macro_colors[sub_av_e]),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(2.8, "mm")
)
ha_row_e <- rowAnnotation(
  `Macro subtype` = row_origin,
  col = list(`Macro subtype` = macro_colors[macro_order]),
  show_annotation_name = FALSE, show_legend = TRUE,
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(
    `Macro subtype` = list(
      title_gp = gpar(fontsize = 7, fontfamily = my_font, fontface = "bold"),
      labels_gp = gpar(fontsize = 6, fontfamily = my_font),
      grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")
    )
  )
)

ht_e <- Heatmap(
  mat_e_mlog,
  name = "-log10(adj.P)", col = col_fun_e,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = ha_col_e, left_annotation = ha_row_e,
  row_split = row_origin,
  row_title = NULL, row_gap = unit(0.8, "mm"),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6.5, fontfamily = my_font),
  column_names_gp = gpar(fontsize = 7, fontfamily = my_font, fontface = "italic"),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.6),
  border = TRUE, border_gp = gpar(col = "black", lwd = 0.8),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontfamily = my_font, fontface = "bold"),
    labels_gp = gpar(fontsize = 6, fontfamily = my_font),
    legend_height = unit(22, "mm"), grid_width = unit(2.5, "mm")
  ),
  width = unit(ncol(mat_e_mlog) * 6, "mm"),
  height = unit(nrow(mat_e_mlog) * 5, "mm")
)
pdf("panels/fig4e_go_heatmap.pdf", width = 6.2, height = 4.8)
ComplexHeatmap::draw(ht_e, padding = unit(c(2, 4, 2, 2), "mm"),
                     merge_legend = TRUE,
                     heatmap_legend_side = "right",
                     annotation_legend_side = "right")
dev.off()
png("panels/fig4e_go_heatmap.png", width = 6.2, height = 4.8, units = "in", res = 300)
ComplexHeatmap::draw(ht_e, padding = unit(c(2, 4, 2, 2), "mm"),
                     merge_legend = TRUE,
                     heatmap_legend_side = "right",
                     annotation_legend_side = "right")
dev.off()
cat("  fig4e saved\n")


#  Fig 4F — Anti-tumor gene heatmap (per-gene z-score, Macro subsets only)
cat("-- fig4f: anti-tumor heatmap --\n")
ant_macro <- subset(ant_genes, subtype %in% macro_order)
mat_f_dt <- dcast(as.data.table(ant_macro), gene ~ subtype,
                   value.var = "mean_log1p", fun.aggregate = mean)
sub_av_f <- intersect(macro_order, colnames(mat_f_dt))
mat_f <- as.matrix(mat_f_dt[, ..sub_av_f])
rownames(mat_f) <- mat_f_dt$gene
mat_fz <- row_zscore(mat_f)
mat_fz <- pmax(pmin(mat_fz, 2), -2)

# Group genes by function
fun_group <- ifelse(grepl("^HLA-|^CD74", rownames(mat_fz)), "MHC-II",
              ifelse(grepl("^IL1|^TNF|^IL6", rownames(mat_fz)), "Pro-inflam.",
              ifelse(grepl("^IL12|^IFNG", rownames(mat_fz)), "Th1",
              ifelse(grepl("^IRF|^STAT", rownames(mat_fz)), "TF",
              ifelse(grepl("^CSF2|^CXCL|^CCL", rownames(mat_fz)), "Chemokine",
                     "Other")))))
fg_levels <- c("MHC-II", "Pro-inflam.", "Th1", "TF", "Chemokine", "Other")
fun_group <- factor(fun_group, levels = fg_levels)
fg_palette <- c("MHC-II" = "#3C5488", "Pro-inflam." = "#E64B35",
                "Th1" = "#00A087", "TF" = "#7E6148",
                "Chemokine" = "#F39B7F", "Other" = "grey80")

# Sort rows by group then alphabetic
ord <- order(fun_group, rownames(mat_fz))
mat_fz <- mat_fz[ord, , drop = FALSE]
fun_group <- fun_group[ord]

col_fun_f <- colorRamp2(c(-2, 0, 2),
                         c("#3C5488", "#F7F7F7", "#E64B35"))

ha_col_f <- HeatmapAnnotation(
  Subtype = sub_av_f,
  col = list(Subtype = macro_colors[sub_av_f]),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(2.8, "mm")
)
ha_row_f <- rowAnnotation(
  Function = fun_group,
  col = list(Function = fg_palette),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Function = list(title_gp = gpar(fontsize = 7, fontfamily = my_font),
                    labels_gp = gpar(fontsize = 6, fontfamily = my_font))
  ),
  simple_anno_size = unit(2.5, "mm")
)

ht_f <- Heatmap(
  mat_fz, name = "z-score", col = col_fun_f,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = ha_col_f, left_annotation = ha_row_f,
  row_split = fun_group, row_title = NULL, row_gap = unit(1, "mm"),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6.5, fontfamily = my_font, fontface = "italic"),
  column_names_gp = gpar(fontsize = 7, fontfamily = my_font),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.4),
  border = TRUE, border_gp = gpar(col = "black", lwd = 1),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontfamily = my_font),
    labels_gp = gpar(fontsize = 7, fontfamily = my_font),
    legend_height = unit(22, "mm"), grid_width = unit(2.5, "mm")
  ),
  width = unit(ncol(mat_fz) * 5, "mm"),
  height = unit(nrow(mat_fz) * 5, "mm")
)
pdf("panels/fig4f_antitumor_heatmap.pdf", width = 5, height = 6.5)
ComplexHeatmap::draw(ht_f, padding = unit(c(2, 12, 2, 2), "mm"),
                     merge_legend = TRUE)
dev.off()
png("panels/fig4f_antitumor_heatmap.png", width = 5, height = 6.5,
    units = "in", res = 300)
ComplexHeatmap::draw(ht_f, padding = unit(c(2, 12, 2, 2), "mm"),
                     merge_legend = TRUE)
dev.off()
cat("  fig4f saved\n")


#  Fig 4G — Macro subset markers (6 subsets combined, 2x3 dotplots)
#           Vertical italic gene labels; shared subtype y-axis on left
cat("-- fig4g: combined macro subset dotplot panel --\n")

build_subset_dotplot <- function(target_subtype, title, show_y_text = TRUE,
                                 show_legend = FALSE) {
  d <- subset(sub_markers, panel_origin == target_subtype)
  if (nrow(d) == 0) return(NULL)
  d$subtype <- factor(d$subtype, levels = rev(sub_order))
  d <- add_gene_zscore(d, "mean_log1p", "gene")
  d$gene <- factor(d$gene, levels = unique(d$gene))
  d$pct100 <- d$pct_expressing * 100

  p <- ggplot(d, aes(x = gene, y = subtype)) +
    geom_point(aes(size = pct100, fill = z_expr),
               shape = 21, color = "black", stroke = 0.2) +
    scale_size_continuous(name = "% expressing",
                          range = c(0.3, 4), breaks = c(25, 50, 75),
                          limits = c(0, 100)) +
    scale_fill_gradient2(name = "z-score",
                         low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                         midpoint = 0, limits = c(-2, 2),
                         oob = scales::squish,
                         guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                barheight = unit(12, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.3,
                                                ticks.colour = "black")) +
    labs(x = NULL, y = NULL, title = title) +
    theme_pub(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 5.5, face = "italic"),
          axis.text.y = if (show_y_text)
            element_text(size = 5.5) else element_blank(),
          axis.ticks.y = if (show_y_text)
            element_line(linewidth = 0.3) else element_blank(),
          plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          legend.position = if (show_legend) "right" else "none",
          legend.box = "vertical",
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5),
          legend.key.size = unit(2.5, "mm"),
          panel.border = element_rect(color = "black", fill = NA,
                                      linewidth = 0.4),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks.x = element_line(linewidth = 0.3),
          plot.margin = margin(2, 2, 2, 2, "pt"))
  p
}

mac_subsets <- list(
  c("Macro_MARCO",   "Macro_MARCO (alveolar-like)"),
  c("Macro_C1QC",    "Macro_C1QC (antigen-presenting)"),
  c("Macro_FCN1",    "Macro_FCN1 (mono-derived)"),
  c("Macro_FOLR2",   "Macro_FOLR2 (tissue-resident)"),
  c("Macro_SPP1",    "Macro_SPP1 (hypoxic TAM)"),
  c("Macro_general", "Macro_general (baseline)")
)

# 2 rows × 3 cols. Show y-axis labels only on leftmost column (1, 4).
# Show legend only on rightmost-bottom (6).
plots_g <- vector("list", length(mac_subsets))
for (i in seq_along(mac_subsets)) {
  show_y  <- i %in% c(1, 4)
  show_lg <- i == 6
  plots_g[[i]] <- build_subset_dotplot(mac_subsets[[i]][1],
                                       mac_subsets[[i]][2],
                                       show_y_text = show_y,
                                       show_legend = show_lg)
}
plots_g <- Filter(Negate(is.null), plots_g)

p4g <- wrap_plots(plots_g, ncol = 3, nrow = 2,
                  widths = c(1.15, 1, 1)) +
  plot_annotation(
    theme = theme(plot.margin = margin(2, 2, 2, 2, "pt"))
  )

ggsave("panels/fig4g_macro_subsets.pdf", p4g,
       width = 200, height = 110, units = "mm")
ggsave("panels/fig4g_macro_subsets.png", p4g,
       width = 200, height = 110, units = "mm", dpi = 300)
cat("  fig4g (combined macro subsets) saved\n")

# Remove obsolete per-subset panel files from any previous run
# (fig4g_marco / fig4h_c1qc ... fig4l_general — replaced by combined fig4g_macro_subsets)
old_basenames <- c("fig4g_marco", "fig4h_c1qc", "fig4i_fcn1",
                   "fig4j_folr2", "fig4k_spp1", "fig4l_general")
for (b in old_basenames) {
  for (ext in c("pdf", "png")) {
    f <- file.path("panels", paste0(b, ".", ext))
    if (file.exists(f)) file.remove(f)
  }
}


#  Fig 4M — M1/M2 score bar + MP association heatmap (combined)
cat("-- fig4m: M1/M2 + MP heatmap --\n")

# Left: M1/M2 grouped bar
m12 <- as.data.table(m12_refined)
setnames(m12, "myeloid_subtype_refined", "subtype")
m12 <- m12[subtype %in% sub_order]
m_long <- rbind(
  m12[, .(subtype, score = "M1", mean = M1_score_mean, sd = M1_score_std)],
  m12[, .(subtype, score = "M2", mean = M2_score_mean, sd = M2_score_std)]
)
ord_sub <- m12[order(-M1_M2_ratio_median), subtype]
m_long$subtype <- factor(m_long$subtype, levels = ord_sub)
m_long$score <- factor(m_long$score, levels = c("M1", "M2"))
m_long <- as.data.frame(m_long)

p4m_left <- ggplot(m_long, aes(x = subtype, y = mean, fill = score)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c("M1" = "#E64B35", "M2" = "#3C5488"),
                    name = "Polarization") +
  labs(x = NULL, y = "Polarization score (mean ± SD)",
       title = "M1/M2 polarization") +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                   face = "italic"),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.position = "top",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        plot.title = element_text(size = 9, face = "bold"))

# Right: MP × subtype Spearman heatmap (geom_tile for patchwork compat)
mp <- as.data.table(mp_refined)
mp <- mp[subtype %in% sub_order & MP %in% c("MP1", "MP2", "MP3", "MP4")]
mp$subtype <- factor(mp$subtype, levels = rev(sub_order))
mp$MP <- factor(mp$MP, levels = c("MP1", "MP2", "MP3", "MP4"))
mp$rho_clamp <- pmax(pmin(mp$spearman_rho, 0.4), -0.4)
mp$star <- stars_fn(mp$p)
mp <- as.data.frame(mp)

p4m_right <- ggplot(mp, aes(x = MP, y = subtype, fill = rho_clamp)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)),
            size = 1.9, color = "black", vjust = 1.4) +
  geom_text(aes(label = star),
            size = 2.2, color = "black", fontface = "bold", vjust = -0.6) +
  scale_fill_gradient2(name = "Spearman ρ",
                       low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                       midpoint = 0, limits = c(-0.4, 0.4),
                       oob = scales::squish,
                       guide = guide_colorbar(barwidth = unit(3, "mm"),
                                              barheight = unit(15, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.3,
                                              ticks.colour = "black")) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL, title = "MP × subtype (patient-level)") +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(size = 8, face = "bold.italic"),
        axis.text.y = element_text(size = 7),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

# Save as two separate files to allow square MP cells
ggsave("panels/fig4m_m1m2.pdf", p4m_left, width = 140, height = 75, units = "mm")
ggsave("panels/fig4m_m1m2.png", p4m_left, width = 140, height = 75, units = "mm",
       dpi = 300)
ggsave("panels/fig4m_mp_heatmap.pdf", p4m_right,
       width = 90, height = 165, units = "mm")
ggsave("panels/fig4m_mp_heatmap.png", p4m_right,
       width = 90, height = 165, units = "mm", dpi = 300)
cat("  fig4m saved\n")


#  Fig 4N — Hallmark GSEA SPP1 vs C1QC (FDR<0.05 filter, NES bar)
cat("-- fig4n: GSEA --\n")
g <- as.data.table(gsea_n)
setnames(g, c("FDR q-val", "NOM p-val"), c("FDR_q", "NOM_p"), skip_absent = TRUE)
g <- g[FDR_q < 0.05]
g[, Term := gsub("^HALLMARK_", "", Term)]
g[, Term := gsub("_", " ", Term)]
g[, direction := ifelse(NES > 0, "SPP1_enriched", "C1QC_enriched")]
g[, direction := factor(direction, levels = c("SPP1_enriched", "C1QC_enriched"))]
g <- g[order(-NES)]
g[, Term := factor(Term, levels = Term[order(NES)])]
g <- as.data.frame(g)

nes_lim <- max(abs(g$NES), na.rm = TRUE)
p4n <- ggplot(g, aes(x = Term, y = NES, fill = NES)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  scale_fill_gradient2(name = "NES",
                       low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                       midpoint = 0,
                       limits = c(-nes_lim, nes_lim),
                       oob = scales::squish,
                       guide = guide_colorbar(barwidth = unit(15, "mm"),
                                              barheight = unit(2.5, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.3,
                                              ticks.colour = "black",
                                              title.position = "top",
                                              title.hjust = 0.5)) +
  coord_flip() +
  labs(x = NULL, y = "Normalized Enrichment Score",
       title = "Hallmark GSEA: SPP1 vs C1QC  (FDR < 0.05)") +
  theme_pub(base_size = 8) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, face = "italic"),
        axis.title.x = element_text(size = 8),
        legend.position = "top",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"))

ggsave("panels/fig4n_gsea.pdf", p4n, width = 150, height = 150, units = "mm")
ggsave("panels/fig4n_gsea.png", p4n, width = 150, height = 150, units = "mm",
       dpi = 300)
cat("  fig4n saved\n")


#  SUMMARY
cat("\n==============================\n")
cat("Figure 4 (14-panel, strict-style) outputs:\n")
files <- list.files("panels", pattern = "^fig4[a-n]_.*\\.(pdf|png)$",
                     full.names = TRUE)
files <- files[order(files)]
for (f in files) {
  sz <- round(file.info(f)$size / 1e6, 2)
  cat(sprintf("  %s (%.2f MB)\n", f, sz))
}
cat("==============================\n")
