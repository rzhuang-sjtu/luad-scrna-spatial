#!/usr/bin/env Rscript
# Fig 8 overhaul:
#   8A: proportional eulerr Venn + count badges (modern, clean)
#   8B: stays (violin LUAD vs other)
#   8C: horizontal lollipop (replaces 2000-style bar)
#   8D: TvN combined 3-facet (replaces old 8D/E/F)
#   ----
#   drop: old 8G scatter (negative finding) and old 8H coexpr embedding
#   ----
#   spatial 8I-L: subtitle = sample name only
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(patchwork); library(ggrepel); library(ggbeeswarm)
  library(png); library(grid); library(jsonlite); library(scales)
  library(ComplexHeatmap)
})

DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

############################################################
# 8A: UpSet plot via ComplexHeatmap (modern, top-journal style for 3+ sets)
############################################################
cat("\n[8A] UpSet plot (ComplexHeatmap)\n")
v <- read_csv(file.path(DATA, "8A_venn_subsets.csv"), show_col_types = FALSE)
v200 <- v[v$top_N == 200, ]
# Schema: only_macro, only_mal, only_neu, macro_AND_mal_only,
#         macro_AND_neu_only, mal_AND_neu_only, all_three
n_macro <- v200$size_macro; n_mal <- v200$size_mal; n_neu <- v200$size_neu
n_mm_only  <- v200$macro_AND_mal_only
n_mn_only  <- v200$macro_AND_neu_only
n_man_only <- v200$mal_AND_neu_only
n_all      <- v200$all_three
# pre-compute combined-only set sizes for the 7 Venn regions
n_mm  <- n_mm_only  + n_all
n_mn  <- n_mn_only  + n_all
n_man <- n_man_only + n_all

# Reconstruct gene sets (synthetic since we only have counts) - use logical matrix
# build a comb_mat directly from sizes
nm <- c("Macro-SPP1->C1QC", "Mal-MP3->MP1", "Neu-OSM_p->OSM_low")
size_vec <- c(
  "100" = n_macro - n_mm - n_mn  + n_all,
  "010" = n_mal   - n_mm - n_man + n_all,
  "001" = n_neu   - n_mn - n_man + n_all,
  "110" = n_mm  - n_all,
  "101" = n_mn  - n_all,
  "011" = n_man - n_all,
  "111" = n_all
)
# Build fake set lists matching these sizes
gene_pool <- paste0("g", seq_len(sum(size_vec)))
ofs <- 0; sets <- list("Macro-SPP1->C1QC"=character(0), "Mal-MP3->MP1"=character(0),
                        "Neu-OSM_p->OSM_low"=character(0))
for (code in names(size_vec)) {
  k <- size_vec[[code]]; if (k == 0) next
  gs <- gene_pool[(ofs + 1):(ofs + k)]; ofs <- ofs + k
  if (substr(code,1,1) == "1") sets[["Macro-SPP1->C1QC"]] <- c(sets[["Macro-SPP1->C1QC"]], gs)
  if (substr(code,2,2) == "1") sets[["Mal-MP3->MP1"]]     <- c(sets[["Mal-MP3->MP1"]], gs)
  if (substr(code,3,3) == "1") sets[["Neu-OSM_p->OSM_low"]] <- c(sets[["Neu-OSM_p->OSM_low"]], gs)
}
m <- ComplexHeatmap::make_comb_mat(sets)

pal <- c("Macro-SPP1->C1QC" = "#3C5488",
         "Mal-MP3->MP1" = "#E64B35",
         "Neu-OSM_p->OSM_low" = "#00A087")

ht <- ComplexHeatmap::UpSet(m,
        set_order = names(sets),
        comb_order = order(-ComplexHeatmap::comb_size(m)),
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          "Intersection" = ComplexHeatmap::anno_barplot(
            ComplexHeatmap::comb_size(m),
            border = FALSE, gp = gpar(fill = "grey25"),
            height = unit(2.4, "cm"),
            add_numbers = TRUE, numbers_gp = gpar(fontsize = 7)),
          annotation_name_gp = gpar(fontsize = 7, fontface = "bold")),
        right_annotation = ComplexHeatmap::rowAnnotation(
          "Set size" = ComplexHeatmap::anno_barplot(
            ComplexHeatmap::set_size(m),
            border = FALSE,
            gp = gpar(fill = pal[names(sets)]),
            width = unit(2.4, "cm"),
            add_numbers = TRUE, numbers_gp = gpar(fontsize = 7)),
          annotation_name_gp = gpar(fontsize = 7, fontface = "bold")),
        row_names_gp = gpar(fontsize = 7, fontface = "bold"),
        comb_col = "grey25", bg_col = "grey94",
        pt_size = unit(3, "mm"), lwd = 1.2)

pdf(file.path(OUT, "8A_venn.pdf"), width = 4.6, height = 2.6, useDingbats = FALSE)
ComplexHeatmap::draw(ht)
dev.off()
png(file.path(OUT, "8A_venn.png"), width = 4.6, height = 2.6, units = "in", res = 300)
ComplexHeatmap::draw(ht)
dev.off()
cat("  -> 8A_venn.pdf (UpSet plot)\n")

############################################################
# 8C : horizontal lollipop (replaces 2D bar)
############################################################
cat("\n[8C] LUAD essentiality lollipop\n")
dep <- read_csv(file.path(DATA, "8B_depmap_long.csv"), show_col_types = FALSE)
dep$is_LUAD <- as.logical(dep$is_LUAD)
luad_dep <- dep %>% filter(is_LUAD) %>%
  group_by(gene) %>%
  summarise(mean_eff = mean(gene_effect, na.rm = TRUE),
            sd_eff   = sd(gene_effect, na.rm = TRUE),
            n        = sum(!is.na(gene_effect))) %>%
  mutate(se = sd_eff / sqrt(n)) %>%
  arrange(mean_eff)
luad_dep$gene <- factor(luad_dep$gene, levels = luad_dep$gene)

p_8C <- ggplot(luad_dep, aes(mean_eff, gene)) +
  geom_vline(xintercept = 0,    color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = -0.5, color = COL$ref_red, linetype = "dashed",
             linewidth = 0.35, alpha = 0.5) +
  geom_segment(aes(x = 0, xend = mean_eff, yend = gene),
               linewidth = 0.45, color = "grey55") +
  geom_errorbarh(aes(xmin = mean_eff - se, xmax = mean_eff + se),
                 height = 0.18, linewidth = 0.4, color = "grey25") +
  geom_point(size = 2.6, color = COL$LUAD, fill = COL$LUAD,
             shape = 21, stroke = 0.4) +
  geom_text(aes(label = sprintf("%.2f", mean_eff)),
            hjust = -0.5, family = FAM, size = 2.1, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Mean Chronos gene effect (LUAD lines)", y = NULL) +
  theme_pub(8) +
  theme(panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin  = margin(2, 4, 2, 2))
save_panel(p_8C, file.path(OUT, "8C_bar_luad_essentiality"), 3.0, 1.9)

############################################################
# 8D : TvN combined 3-facet (replaces 8D/E/F)
############################################################
cat("\n[8D] TvN combined 3-facet\n")
tn <- read_csv(file.path(DATA, "8D_tcga_TvN_long.csv"), show_col_types = FALSE)
tn_st <- read_csv(file.path(DATA, "8D_tcga_TvN_stats.csv"), show_col_types = FALSE)
tn$gene <- factor(tn$gene, levels = LEADS)
tn$type <- factor(tn$type, levels = c("Normal","Tumor"))
ann_d <- tn_st %>% mutate(gene = factor(gene, levels = LEADS),
                           label = sprintf("log2FC=%+.2f\n%s",
                                            log2FC_T_minus_N,
                                            fmt_p(wilcoxon_p)))
ymax_d <- max(tn$log2_TPM_p1, na.rm = TRUE)

p_8D <- ggplot(tn, aes(type, log2_TPM_p1, fill = type, color = type)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35,
              linewidth = 0.3, color = NA) +
  geom_quasirandom(width = 0.22, size = 0.32, alpha = 0.55,
                   shape = 16, stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.40,
               linewidth = 0.45, fatten = 1.6, color = "black") +
  scale_fill_manual(values = c(Normal = COL$normal, Tumor = COL$tumor),
                    guide = "none") +
  scale_color_manual(values = c(Normal = COL$normal, Tumor = COL$tumor),
                     guide = "none") +
  geom_text(data = ann_d, aes(x = 1.5, y = ymax_d * 1.04, label = label),
            inherit.aes = FALSE, family = FAM, size = 2.0) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(TPM+1)",
       subtitle = "TCGA-LUAD Tumor vs Normal") +
  theme_pub(8) +
  theme(plot.subtitle = element_text(size = 7.5, color = "grey25",
                                     margin = margin(b = 1)),
        strip.text = element_text(size = 7.5, face = "bold"),
        panel.spacing = unit(2, "mm"))
save_panel(p_8D, file.path(OUT, "8D_TvN_combined"), 6.0, 2.2)

############################################################
# Spatial subtitles -> sample name only
############################################################
spot <- read_csv(file.path(DATA, "8I_spot_long.csv"), show_col_types = FALSE)
spot$gene <- factor(spot$gene, levels = LEADS)
HE_ROOTS <- list(
  `E-MTAB-13530` = list(
    dir = "${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530",
    sub = function(s) sprintf("%s-spatial", s)),
  Okamura = list(
    dir = "${DATA_ROOT}/ST/results/step09_okamura_validation/raw",
    sub = function(s) sprintf("%s/spatial", s)))

load_he <- function(cohort, section) {
  cfg <- HE_ROOTS[[cohort]]; if (is.null(cfg)) return(NULL)
  base <- file.path(cfg$dir, cfg$sub(section))
  img_p <- file.path(base, "tissue_hires_image.png")
  sf_p  <- file.path(base, "scalefactors_json.json")
  if (!file.exists(img_p) || !file.exists(sf_p)) return(NULL)
  list(img = png::readPNG(img_p),
       w = dim(png::readPNG(img_p))[2], h = dim(png::readPNG(img_p))[1],
       sf = jsonlite::fromJSON(sf_p)$tissue_hires_scalef)
}
GENE_HIGH <- c(SEC61G = "#FFB300", SRSF9 = "#00BCD4", ANGPTL4 = "#43A047")

plot_gene_he <- function(d, he, gene, hi) {
  v <- d$expr; lo <- quantile(v, 0.02, na.rm=TRUE); hh <- quantile(v, 0.98, na.rm=TRUE)
  rng <- hh - lo
  d$alpha_pt <- 0.92 * pmin(1, pmax(0, (v - lo) / ifelse(rng > 0, rng, 1)))
  ggplot() +
    annotation_raster(he$img, xmin = 0, xmax = he$w, ymin = -he$h, ymax = 0) +
    geom_point(data = d, aes(x_he, y_he, color = expr, alpha = alpha_pt),
               size = 0.85, stroke = 0, shape = 16) +
    scale_color_gradient(low = "white", high = hi,
                         limits = c(lo, hh), oob = scales::squish,
                         name = sprintf("%s\nexpr", gene),
                         guide = guide_colorbar(barwidth = unit(2.2, "mm"),
                                                barheight = unit(16, "mm"),
                                                frame.colour = "black",
                                                frame.linewidth = 0.25,
                                                title.position = "top")) +
    scale_alpha_identity() +
    coord_fixed(xlim = c(0, he$w), ylim = c(-he$h, 0), expand = FALSE) +
    theme_pub(8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(1, 1, 1, 1),
          legend.position = "right",
          legend.margin = margin(0,0,0,1), legend.box.margin = margin(0,0,0,0),
          legend.title = element_text(size = 7, family = FAM, face = "bold"),
          legend.text  = element_text(size = 6, family = FAM),
          plot.subtitle = element_text(size = 7, color = "black",
                                       face = "bold", margin = margin(b = 1))) +
    labs(subtitle = gene)
}
draw_section <- function(cohort, section, out_stem) {
  sub <- spot %>% filter(cohort == !!cohort, sample == section)
  he <- load_he(cohort, section)
  if (is.null(he)) { cat(sprintf("  WARN: no H&E %s/%s\n", cohort, section)); return() }
  sub <- sub %>% mutate(x_he = spatial1 * he$sf, y_he = -spatial2 * he$sf)
  ps <- lapply(LEADS, function(g) plot_gene_he(sub %>% filter(gene == g), he, g,
                                                unname(GENE_HIGH[g])))
  p <- wrap_plots(ps, nrow = 1) +
        plot_annotation(
          subtitle = section,
          theme = theme(plot.subtitle = element_text(size = 8, family = FAM,
                                                     face = "bold", hjust = 0.02,
                                                     margin = margin(b = 1)),
                        plot.margin = margin(1, 1, 1, 1)))
  save_panel(p, out_stem, 6.0, 2.1)
  cat(sprintf("  saved %s\n", basename(out_stem)))
}

cat("\n[spatial] re-render with sample-name-only subtitle\n")
draw_section("E-MTAB-13530", "P10_T1",     file.path(OUT, "8I_spatial_P10_T1"))
draw_section("E-MTAB-13530", "P15_T1",     file.path(OUT, "8J_spatial_P15_T1"))
draw_section("Okamura",       "LUAD_No_4", file.path(OUT, "8K_spatial_LUAD_No_4"))
draw_section("Okamura",       "LUAD_No_1", file.path(OUT, "8L_spatial_LUAD_No_1"))

cat("\nDONE.\n")
