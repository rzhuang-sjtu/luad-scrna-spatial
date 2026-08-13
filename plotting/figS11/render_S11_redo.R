#!/usr/bin/env Rscript
# Redraw Figure S11 panels using the unified Fig 8 raincloud template
# (matching 8D / 8I square 3-facet style).
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(ggbeeswarm)
})
DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig_s11"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

############################################################
# S11D : merged GSE135222 R vs NR (single 3-facet, replaces D/E/F)
############################################################
cat("\n[S11D] merged GSE135222 R vs NR\n")
g3  <- read_csv(file.path(DATA, "S11D_GSE135222_long.csv"),  show_col_types = FALSE)
g3s <- read_csv(file.path(DATA, "S11D_GSE135222_stats.csv"), show_col_types = FALSE)
g3$gene <- factor(g3$gene, levels = LEADS)
g3$response <- factor(g3$response, levels = c("NR", "R"))
y_top <- max(g3$expr, na.rm = TRUE) * 1.04
ann_d <- g3s %>%
  mutate(gene = factor(gene, levels = LEADS),
         label = sprintf("log2FC=%+.2f\n%s",
                          log2FC_pos_minus_neg, fmt_p(p)),
         y_lab = y_top)
p_d <- ggplot(g3, aes(response, expr, fill = response, color = response)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
              color = NA) +
  geom_quasirandom(width = 0.22, size = 0.32, alpha = 0.55, shape = 16,
                   stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.40,
               linewidth = 0.45, fatten = 1.6, color = "black") +
  scale_fill_manual(values  = c(NR = COL$NR, R = COL$R), guide = "none") +
  scale_color_manual(values = c(NR = COL$NR, R = COL$R), guide = "none") +
  geom_text(data = ann_d, aes(x = 1.5, y = y_lab, label = label),
            inherit.aes = FALSE, family = FAM, size = 2.0) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(TPM+1)",
       subtitle = sprintf("GSE135222 single-cell - R (n=%d) vs NR (n=%d)",
                           g3s$n_pos[1], g3s$n_neg[1])) +
  theme_pub(8) +
  theme(plot.subtitle = element_text(size = 7.5, color = "grey25",
                                      margin = margin(b = 1)),
        strip.text = element_text(size = 8, face = "bold",
                                   margin = margin(b = 1, t = 1)),
        strip.background = element_rect(fill = "grey94", color = NA),
        panel.spacing = unit(2, "mm"))
save_panel(p_d, file.path(OUT, "S11D_GSE135222_combined"), 6.0, 2.4)

############################################################
# S11E : tumor-intrinsic ROI raincloud, 2 cohorts x 3 genes
# (was S11G previously; renamed to S11E to fill freed E/F slots)
############################################################
cat("\n[S11E] tumor-intrinsic ROI 2x3\n")
ti  <- read_csv(file.path(DATA, "tumor_intrinsic_roi_long.csv"), show_col_types = FALSE)
tist<- read_csv(file.path(DATA, "tumor_intrinsic_roi_stats.csv"), show_col_types = FALSE)
ti$gene    <- factor(ti$gene, levels = LEADS)
ti$dataset <- factor(ti$dataset, levels = c("E-MTAB-13530","Okamura"),
                                  labels = c("E-MTAB-13530","Takano 2024"))
ti$grp <- factor(ifelse(ti$new_roi, "ROI", "non-ROI"),
                  levels = c("non-ROI","ROI"))
ann_e <- tist %>%
  mutate(gene = factor(gene, levels = LEADS),
         dataset = factor(dataset, levels = c("E-MTAB-13530","Okamura"),
                           labels = c("E-MTAB-13530","Takano 2024")),
         label = sprintf("Delta=%+.2f\n%s",
                          delta_new_minus_non, fmt_p(mw_p)))
y_top_e <- max(ti$expr, na.rm = TRUE) * 1.04

p_e <- ggplot(ti, aes(grp, expr, fill = grp, color = grp)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
              color = NA) +
  geom_quasirandom(width = 0.22, size = 0.16, alpha = 0.28, shape = 16,
                   stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.42,
               linewidth = 0.40, fatten = 1.5, color = "black") +
  scale_fill_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI),
                    guide = "none") +
  scale_color_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI),
                     guide = "none") +
  geom_text(data = ann_e, aes(x = 1.5, y = y_top_e, label = label),
            inherit.aes = FALSE, family = FAM, size = 1.9) +
  facet_grid(dataset ~ gene, scales = "free_y", switch = "y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(TPM+1)",
       subtitle = "Tumor-intrinsic ROI vs non-ROI - z(Malignant)>0.5 AND z(MP3)>0.5") +
  theme_pub(8) +
  theme(plot.subtitle = element_text(size = 7.5, color = "grey25",
                                      margin = margin(b = 1)),
        strip.text   = element_text(size = 8, face = "bold",
                                     margin = margin(b = 1, t = 1)),
        strip.background = element_rect(fill = "grey94", color = NA),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        panel.spacing = unit(2, "mm"))
save_panel(p_e, file.path(OUT, "S11E_tumor_intrinsic_ROI"), 6.0, 4.0)

############################################################
# Archive the legacy S11D-G separate-panel files
############################################################
ARCH <- file.path(OUT, "_archived_legacy_panels")
dir.create(ARCH, showWarnings = FALSE, recursive = TRUE)
for (old in c("S11D_GSE135222_SEC61G", "S11E_GSE135222_SRSF9",
              "S11F_GSE135222_ANGPTL4", "S11G_tumor_intrinsic_ROI")) {
  for (ext in c("pdf","png")) {
    src <- file.path(OUT, paste0(old, ".", ext))
    if (file.exists(src)) {
      dst <- file.path(ARCH, basename(src))
      file.rename(src, dst)
      cat(sprintf("  archived %s\n", basename(src)))
    }
  }
}

############################################################
# Renumber downstream panels (H -> F, K-R -> G-N) so the
# letter scheme has no holes after the merge.
############################################################
SHIFTS <- list(
  c("S11H_invasive_front_composite",        "S11F_invasive_front_composite"),
  c("S11H_cooccurrence_curves_supp",        "S11F_cooccurrence_curves_supp"),
  c("S11K_atlas_EMTAB13530_SEC61G",         "S11G_atlas_EMTAB13530_SEC61G"),
  c("S11L_atlas_EMTAB13530_SRSF9",          "S11H_atlas_EMTAB13530_SRSF9"),
  c("S11M_atlas_EMTAB13530_ANGPTL4",        "S11I_atlas_EMTAB13530_ANGPTL4"),
  c("S11N_atlas_Takano_SEC61G",             "S11J_atlas_Takano_SEC61G"),
  c("S11O_atlas_Takano_SRSF9",              "S11K_atlas_Takano_SRSF9"),
  c("S11P_atlas_Takano_ANGPTL4",            "S11L_atlas_Takano_ANGPTL4"),
  c("S11Q_lead_gene_cohort_bubble",         "S11M_lead_gene_cohort_bubble"),
  c("S11R_lead_gene_per_section_heatmap",   "S11N_lead_gene_per_section_heatmap")
)
for (sh in SHIFTS) {
  old_stem <- sh[[1]]; new_stem <- sh[[2]]
  for (ext in c("pdf","png")) {
    src <- file.path(OUT, paste0(old_stem, ".", ext))
    dst <- file.path(OUT, paste0(new_stem, ".", ext))
    if (file.exists(src)) {
      file.rename(src, dst)
      cat(sprintf("  rename %s -> %s\n", basename(src), basename(dst)))
    }
  }
}

cat("\nDONE.\n")
