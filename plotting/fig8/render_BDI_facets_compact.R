#!/usr/bin/env Rscript
# Restore 3-facet layout (one sub-axis per gene), but compact:
# narrower width, taller height, tight facet spacing.
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(ggbeeswarm)
})
DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

facet_panel <- function(d, group_col, group_levels, fill_col, color_col,
                         ann_df, ylab, subtitle_txt, out_stem,
                         hlines = NULL) {
  d[[group_col]] <- factor(d[[group_col]], levels = group_levels)
  yvar <- intersect(c("gene_effect","log2_TPM_p1","expr"), colnames(d))[1]
  p <- ggplot(d, aes(x = .data[[group_col]], y = .data[[yvar]],
                      fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_violin(width = 0.78, trim = FALSE, alpha = 0.35,
                linewidth = 0.3, color = NA) +
    geom_quasirandom(width = 0.20, size = 0.18, alpha = 0.32,
                     shape = 16, stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.42,
                 linewidth = 0.4, fatten = 1.5, color = "black") +
    scale_fill_manual(values  = fill_col,  guide = "none") +
    scale_color_manual(values = color_col, guide = "none")
  if (!is.null(hlines)) {
    for (h in hlines) {
      p <- p + geom_hline(yintercept = h$y, linetype = h$lty,
                            color = h$col, linewidth = 0.32, alpha = 0.55)
    }
  }
  p <- p +
    geom_text(data = ann_df, aes(x = 1.5, y = y_lab, label = label),
              inherit.aes = FALSE, family = FAM, size = 2.0) +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = ylab, subtitle = subtitle_txt) +
    theme_pub(7.5) +
    theme(plot.subtitle = element_text(size = 7, color = "grey25",
                                        margin = margin(b = 1)),
          strip.text    = element_text(size = 7.5, face = "bold",
                                        margin = margin(b = 1, t = 1)),
          strip.background = element_rect(fill = "grey94", color = NA),
          axis.text.x   = element_text(size = 6.5),
          axis.text.y   = element_text(size = 6.5),
          axis.title.y  = element_text(size = 7.5),
          panel.spacing = unit(1.0, "mm"),
          plot.margin   = margin(2, 3, 2, 2))
  # Compact: narrower (3.6 in) + taller (2.6 in)
  save_panel(p, out_stem, 3.6, 2.6)
  cat(sprintf("  saved %s\n", basename(out_stem)))
}

############################################################
# 8B
############################################################
cat("\n[8B]\n")
dep <- read_csv(file.path(DATA, "8B_depmap_long.csv"), show_col_types = FALSE)
dep_st <- read_csv(file.path(DATA, "8B_depmap_stats.csv"), show_col_types = FALSE)
dep$is_LUAD <- as.logical(dep$is_LUAD)
dep$gene <- factor(dep$gene, levels = LEADS)
y_top_b <- max(dep$gene_effect, na.rm = TRUE) * 1.05
ann_b <- dep_st %>%
  mutate(gene = factor(gene, levels = LEADS),
         label = sprintf("Delta=%+.2f\n%s",
                          delta_LUAD_minus_other,
                          fmt_p(mw_p_LUAD_lt_other)),
         y_lab = y_top_b)
facet_panel(
  d = dep, group_col = "group", group_levels = c("LUAD","non-LUAD"),
  fill_col  = c("LUAD" = COL$LUAD, "non-LUAD" = COL$normal),
  color_col = c("LUAD" = COL$LUAD, "non-LUAD" = COL$normal),
  ann_df = ann_b, ylab = "CRISPR Gene Effect",
  subtitle_txt = sprintf("LUAD n=%d vs non-LUAD n=%d - DepMap 24Q2",
                          dep_st$luad_n[1], dep_st$other_n[1]),
  out_stem = file.path(OUT, "8B_violin_LUAD_vs_other"),
  hlines = list(list(y = 0,    lty = "dotted", col = "black"),
                list(y = -0.5, lty = "dashed", col = COL$ref_red))
)

############################################################
# 8D
############################################################
cat("\n[8D]\n")
tn <- read_csv(file.path(DATA, "8D_tcga_TvN_long.csv"), show_col_types = FALSE)
tn_st <- read_csv(file.path(DATA, "8D_tcga_TvN_stats.csv"), show_col_types = FALSE)
tn$gene <- factor(tn$gene, levels = LEADS)
tn$type <- factor(tn$type, levels = c("Normal","Tumor"))
y_top_d <- max(tn$log2_TPM_p1, na.rm = TRUE) * 1.05
ann_d <- tn_st %>%
  mutate(gene = factor(gene, levels = LEADS),
         label = sprintf("log2FC=%+.2f\n%s", log2FC_T_minus_N, fmt_p(wilcoxon_p)),
         y_lab = y_top_d)
facet_panel(
  d = tn, group_col = "type", group_levels = c("Normal","Tumor"),
  fill_col  = c("Normal" = COL$normal, "Tumor" = COL$tumor),
  color_col = c("Normal" = COL$normal, "Tumor" = COL$tumor),
  ann_df = ann_d, ylab = "log2(TPM+1)",
  subtitle_txt = "TCGA-LUAD Tumor vs Normal",
  out_stem = file.path(OUT, "8D_TvN_combined")
)

############################################################
# 8I
############################################################
cat("\n[8I]\n")
m  <- read_csv(file.path(DATA, "8M_GSE207422_long.csv"), show_col_types = FALSE)
ms <- read_csv(file.path(DATA, "8M_GSE207422_stats.csv"), show_col_types = FALSE)
m$gene <- factor(m$gene, levels = LEADS)
m$response <- factor(m$response, levels = c("NMPR","MPR"))
y_top_i <- max(m$expr, na.rm = TRUE) * 1.05
ann_i <- ms %>%
  mutate(gene = factor(gene, levels = LEADS),
         label = sprintf("AUC=%.2f\n%s", auc_pos_vs_neg, fmt_p(p)),
         y_lab = y_top_i)
facet_panel(
  d = m, group_col = "response", group_levels = c("NMPR","MPR"),
  fill_col  = c("NMPR" = COL$NMPR, "MPR" = COL$MPR),
  color_col = c("NMPR" = COL$NMPR, "MPR" = COL$MPR),
  ann_df = ann_i, ylab = "log2(TPM+1)",
  subtitle_txt = sprintf("GSE207422 - MPR (n=%d) vs NMPR (n=%d) - neoadjuvant chemo-IO",
                         ms$n_pos[1], ms$n_neg[1]),
  out_stem = file.path(OUT, "8I_GSE207422_3genes")
)

cat("\nDONE.\n")
