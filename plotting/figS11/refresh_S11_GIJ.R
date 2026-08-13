#!/usr/bin/env Rscript
# Refresh S11G (invasive front composite), S11I (cohort bubble), S11J (per-section heatmap).
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(scales); library(grid); library(patchwork)
})
DST <- "${WORK_ROOT}/luad_figures/fig_s11"
DATA <- "${DATA_ROOT}/ST/results/r_data"

############################################################
# S11G : refresh invasive_front composite
############################################################
cat("\n[S11G] refreshing invasive_front composite\n")
src_pdf <- "${DATA_ROOT}/ST/results/step_invasive_front/plots/S11H_invasive_front_composite.pdf"
src_png <- "${DATA_ROOT}/ST/results/step_invasive_front/plots/S11H_invasive_front_composite.png"
if (!file.exists(src_pdf)) {
  cat("  source missing - re-running step_invasive_front_02_plot.R\n")
  source("${PROJECT_ROOT}/plotting/fig8/step_invasive_front_02_plot.R")
}
if (file.exists(src_pdf)) {
  file.copy(src_pdf, file.path(DST, "S11G_invasive_front_composite.pdf"), overwrite = TRUE)
  file.copy(src_png, file.path(DST, "S11G_invasive_front_composite.png"), overwrite = TRUE)
  cat("  -> S11G_invasive_front_composite.pdf (refreshed)\n")
}

############################################################
# S11I + S11J : cohort bubble + per-section heatmap
############################################################
cat("\n[S11I/J] re-rendering bubble + per-section heatmap\n")
agg <- read.csv(file.path(DATA, "roi_vs_nonroi_aggregate_pvalues.csv"))
sec <- read.csv(file.path(DATA, "roi_vs_nonroi_stats_with_pvalues.csv"))
agg$cohort <- ifelse(agg$cohort == "Okamura", "Takano 2024",
                ifelse(agg$cohort == "EMTAB13530", "E-MTAB-13530", agg$cohort))
sec$cohort <- ifelse(sec$cohort == "Okamura", "Takano 2024",
                ifelse(sec$cohort == "EMTAB13530", "E-MTAB-13530", sec$cohort))
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")
agg_l <- agg[grepl("^gex_", agg$metric), ]
agg_l$gene <- sub("^gex_", "", agg_l$metric)
agg_l <- agg_l[agg_l$gene %in% LEADS, ]

if (nrow(agg_l) == 0) {
  # the lead-gene aggregate may live elsewhere - try the fig8 lead-gene tables
  agg_l <- read.csv("${WORK_ROOT}/luad_figures/fig8/v2_500/data/per_cohort_lead_gene_sig.csv")
  agg_l$cohort <- ifelse(agg_l$cohort == "Okamura", "Takano 2024",
                  ifelse(agg_l$cohort == "EMTAB13530", "E-MTAB-13530", agg_l$cohort))
}
agg_l$gene <- factor(agg_l$gene, levels = LEADS)
agg_l$cohort <- factor(agg_l$cohort, levels = c("E-MTAB-13530", "Takano 2024"))
agg_l$sig <- ifelse(is.na(agg_l$p_fdr), "",
              ifelse(agg_l$p_fdr < 1e-3, "***",
              ifelse(agg_l$p_fdr < 1e-2, "**",
              ifelse(agg_l$p_fdr < 5e-2, "*", ""))))
fill_lim <- max(abs(agg_l$delta), na.rm = TRUE)

p_bubble <- ggplot(agg_l, aes(cohort, gene)) +
  geom_point(aes(size = abs(delta), fill = delta),
             shape = 21, color = "grey25", stroke = 0.3) +
  geom_text(aes(label = sig), color = "white",
            family = FAM, size = 2.0, fontface = "bold",
            hjust = 0.5, vjust = 0.55) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, limits = c(-fill_lim, fill_lim),
                       oob = squish, name = "Delta",
                       guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                              barheight = unit(20, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.3)) +
  scale_size_continuous(range = c(2.5, 8), name = "|Delta|") +
  labs(x = NULL, y = NULL,
       subtitle = "Tumor-intrinsic ROI vs non-ROI - cohort-pooled Mann-Whitney + BH-FDR\n*** <0.001  ** <0.01  * <0.05") +
  theme_pub(8) +
  theme(axis.text.x   = element_text(face = "bold", size = 8),
        axis.text.y   = element_text(face = "italic", size = 8),
        plot.subtitle = element_text(size = 7, color = "grey25",
                                      lineheight = 1.1, margin = margin(b = 1)),
        legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text  = element_text(size = 6),
        legend.box   = "vertical")
save_panel(p_bubble, file.path(DST, "S11I_lead_gene_cohort_bubble"), 4.0, 2.0)

# S11J per-section heatmap (lead genes only)
sec_l <- read.csv("${WORK_ROOT}/luad_figures/fig8/v2_500/data/per_section_lead_gene_sig.csv")
sec_l$cohort <- ifelse(sec_l$cohort == "Okamura", "Takano 2024",
                ifelse(sec_l$cohort == "EMTAB13530", "E-MTAB-13530", sec_l$cohort))
sec_l$gene   <- factor(sec_l$gene,   levels = LEADS)
sec_l$cohort <- factor(sec_l$cohort, levels = c("E-MTAB-13530","Takano 2024"))
sample_order <- sec_l %>%
  group_by(sample) %>% summarise(s = mean(delta, na.rm = TRUE)) %>%
  arrange(-s) %>% pull(sample)
sec_l$sample <- factor(sec_l$sample, levels = sample_order)
fill_cap <- as.numeric(quantile(abs(sec_l$delta), 0.95, na.rm = TRUE))

p_hm <- ggplot(sec_l, aes(sample, gene)) +
  geom_tile(aes(fill = delta), color = "white", linewidth = 0.18) +
  geom_text(aes(label = sig), color = "black",
            family = FAM, size = 1.9, fontface = "bold",
            hjust = 0.5, vjust = 0.55) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, limits = c(-fill_cap, fill_cap),
                       oob = squish, name = "Delta\n(ROI - non-ROI)",
                       guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                              barheight = unit(20, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.3,
                                              title.position = "top")) +
  facet_grid(. ~ cohort, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       subtitle = "Per-section tumor-intrinsic ROI vs non-ROI - Mann-Whitney + per-section BH-FDR\n*** <0.001  ** <0.01  * <0.05") +
  theme_pub(8) +
  theme(panel.grid       = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                         size = 6, color = "black"),
        axis.text.y      = element_text(face = "italic", size = 7,
                                         color = "black"),
        axis.ticks       = element_blank(),
        strip.background = element_rect(fill = "grey94", color = NA),
        strip.text       = element_text(face = "bold", size = 7.5),
        plot.subtitle    = element_text(size = 7, color = "grey25",
                                         lineheight = 1.1, margin = margin(b = 1)),
        legend.position  = "right",
        legend.title     = element_text(size = 7, face = "bold"),
        legend.text      = element_text(size = 6),
        plot.margin      = margin(2, 4, 2, 2))
save_panel(p_hm, file.path(DST, "S11J_lead_gene_per_section_heatmap"), 5.5, 2.0)

cat("\nDONE.\n")
