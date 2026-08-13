#!/usr/bin/env Rscript
# Add a horizontal (single-row) variant of S11H: 6 facets in 1 row
# (3 genes x 2 cohorts), saved as S11H_tumor_intrinsic_ROI_horizontal.pdf.
# The original 2x3 (cohort rows / gene cols) version is kept untouched.
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(ggbeeswarm)
})
DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig_s11"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

ti  <- read_csv(file.path(DATA, "tumor_intrinsic_roi_long.csv"), show_col_types = FALSE)
tist<- read_csv(file.path(DATA, "tumor_intrinsic_roi_stats.csv"), show_col_types = FALSE)
ti$gene    <- factor(ti$gene, levels = LEADS)
ti$dataset <- factor(ti$dataset, levels = c("E-MTAB-13530","Okamura"),
                                  labels = c("E-MTAB-13530","Takano 2024"))
ti$grp <- factor(ifelse(ti$new_roi, "ROI", "non-ROI"),
                  levels = c("non-ROI","ROI"))

# Compose a single grouping label per facet: "gene | cohort"
ti$facet_lbl <- factor(
  paste(ti$gene, ti$dataset, sep = " | "),
  levels = unlist(lapply(LEADS, function(g)
    paste(g, c("E-MTAB-13530","Takano 2024"), sep = " | "))))

ann <- tist %>%
  mutate(gene    = factor(gene, levels = LEADS),
         dataset = factor(dataset, levels = c("E-MTAB-13530","Okamura"),
                                    labels = c("E-MTAB-13530","Takano 2024")),
         facet_lbl = factor(paste(gene, dataset, sep = " | "),
                             levels = levels(ti$facet_lbl)),
         label = sprintf("Delta=%+.2f\n%s",
                          delta_new_minus_non, fmt_p(mw_p)))
y_top <- max(ti$expr, na.rm = TRUE) * 1.04
ann$y_lab <- y_top

p <- ggplot(ti, aes(grp, expr, fill = grp, color = grp)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35,
              linewidth = 0.3, color = NA) +
  geom_quasirandom(width = 0.20, size = 0.16, alpha = 0.28,
                   shape = 16, stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.42,
               linewidth = 0.40, fatten = 1.5, color = "black") +
  scale_fill_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI),
                    guide = "none") +
  scale_color_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI),
                     guide = "none") +
  geom_text(data = ann, aes(x = 1.5, y = y_lab, label = label),
            inherit.aes = FALSE, family = FAM, size = 1.9) +
  facet_wrap(~ facet_lbl, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(TPM+1)",
       subtitle = "Tumor-intrinsic ROI vs non-ROI - 3 lead genes x 2 cohorts (horizontal)") +
  theme_pub(7.5) +
  theme(plot.subtitle = element_text(size = 7, color = "grey25",
                                      margin = margin(b = 1)),
        strip.text   = element_text(size = 7.5, face = "bold",
                                     margin = margin(t = 1, b = 1)),
        strip.background = element_rect(fill = "grey94", color = NA),
        axis.text    = element_text(size = 6.5, color = "black"),
        panel.spacing = unit(1.2, "mm"))

save_panel(p, file.path(OUT, "S11H_tumor_intrinsic_ROI_horizontal"), 9.5, 2.2)
cat("S11H horizontal variant saved (9.5 x 2.2 in, single-row 6-facet)\n")
