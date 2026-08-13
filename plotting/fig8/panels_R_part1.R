#!/usr/bin/env Rscript
# Fig 8 v2 (500-cell) — panels 8A, 8B, 8C, 8D-F, 8G, 8O, 8P, S11D-F.
# Loads pre-exported CSVs from ${PROJECT_ROOT}/results/fig8_plot_data/v2_500/
# Outputs to ${WORK_ROOT}/luad_figures/fig8/v2_500/

suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(VennDiagram)
  library(grid)
  library(survival)
  if (requireNamespace("survminer", quietly = TRUE)) library(survminer)
  library(ggrepel)
})

DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

cat("\n[8A] Venn + barplot\n")
v <- read_csv(file.path(DATA, "8A_venn_subsets.csv"), show_col_types = FALSE)
v200 <- v[v$top_N == 200, ]
n_macro <- v200$size_macro; n_mal <- v200$size_mal; n_neu <- v200$size_neu

# Venn — use VennDiagram::draw.triple.venn (returns grobs; render via grid alone).
# Bigger canvas + shorter labels so circle labels fit without clipping.
draw_venn_grob <- function() {
  grid.newpage()
  draw.triple.venn(
    area1 = n_macro, area2 = n_mal, area3 = n_neu,
    n12 = v200$macro_AND_mal_only + v200$all_three,
    n13 = v200$macro_AND_neu_only + v200$all_three,
    n23 = v200$mal_AND_neu_only   + v200$all_three,
    n123 = v200$all_three,
    category = c("Macro_SPP1\nto C1QC", "Mal_MP3\nto MP1", "Neu_OSM\npriming to low"),
    fill = c(COL$venn_macro, COL$venn_mal, COL$venn_neu),
    alpha = 0.55, lty = "solid", lwd = 0.8, col = "black",
    cex = 0.85, cat.cex = 0.7, cat.fontfamily = FAM, fontfamily = FAM,
    cat.col = "black", margin = 0.12,
    cat.dist = c(0.10, 0.10, 0.06)
  )
}
for (ext in c("png", "pdf")) {
  fp_e <- file.path(OUT, paste0("8A_top200_venn.", ext))
  if (ext == "png") png(fp_e, width = 3.6 * 300, height = 3.2 * 300, res = 300, bg = "white")
  else              pdf(fp_e, width = 3.6, height = 3.2, bg = "white")
  draw_venn_grob()
  dev.off()
}

# 8A right-side barplot — separate ggplot panel
bar_df <- data.frame(label = factor(c("Macro", "Mal", "Neu"), levels = c("Neu","Mal","Macro")),
                      size  = c(n_macro, n_mal, n_neu),
                      col   = c(COL$venn_macro, COL$venn_mal, COL$venn_neu))
p_8A_bar <- ggplot(bar_df, aes(size, label, fill = col)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = size), hjust = -0.2, family = FAM, size = 2.4) +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "# candidate genes", y = NULL,
       subtitle = "top-200 per transition") +
  theme_pub(8)
save_panel(p_8A_bar, file.path(OUT, "8A_pool_size_bar"), 2.4, 1.6)
cat("  wrote 8A venn + bar\n")

cat("\n[8B] DepMap violin\n")
dep <- read_csv(file.path(DATA, "8B_depmap_long.csv"), show_col_types = FALSE)
dep_st <- read_csv(file.path(DATA, "8B_depmap_stats.csv"), show_col_types = FALSE)
dep$is_LUAD <- as.logical(dep$is_LUAD)
dep$gene <- factor(dep$gene, levels = LEADS)
dep$group <- factor(dep$group, levels = c("LUAD", "non-LUAD"))
ann <- dep_st %>%
  mutate(label = sapply(mw_p_LUAD_lt_other, sig_stars),
         y = 1.2,
         gene = factor(gene, levels = LEADS))
p_8B <- ggplot(dep, aes(gene, gene_effect, fill = group)) +
  geom_violin(scale = "width", width = 0.85, position = position_dodge(0.9),
              linewidth = 0.4, color = "black", trim = FALSE) +
  geom_boxplot(width = 0.18, position = position_dodge(0.9), outlier.size = 0.3,
               linewidth = 0.3, color = "black", fill = "white") +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.4) +
  geom_hline(yintercept = -0.5, color = COL$ref_red, linetype = "dashed",
             linewidth = 0.4, alpha = 0.6) +
  geom_text(data = ann, aes(x = gene, y = y, label = label),
            inherit.aes = FALSE, family = FAM, size = 2.4) +
  scale_fill_manual(values = c("LUAD" = COL$LUAD, "non-LUAD" = COL$other),
                    name = NULL) +
  labs(x = NULL, y = "CRISPR Gene Effect",
       subtitle = sprintf("LUAD n=%d vs non-LUAD n=%d · DepMap 24Q2",
                          dep_st$luad_n[1], dep_st$other_n[1])) +
  theme_pub(8) +
  theme(legend.position = "right")
save_panel(p_8B, file.path(OUT, "8B_violin_LUAD_vs_other"), 5, 2.8)

cat("\n[8C] LUAD-only essentiality bar\n")
luad_dep <- dep %>% filter(is_LUAD) %>%
  group_by(gene) %>%
  summarise(mean_eff = mean(gene_effect, na.rm = TRUE),
            sd_eff   = sd(gene_effect, na.rm = TRUE),
            n        = sum(!is.na(gene_effect))) %>%
  mutate(se = sd_eff / sqrt(n)) %>%
  arrange(mean_eff)
luad_dep$gene <- factor(luad_dep$gene, levels = luad_dep$gene)
p_8C <- ggplot(luad_dep, aes(gene, mean_eff)) +
  geom_col(fill = COL$LUAD, color = "black", linewidth = 0.4, width = 0.65) +
  geom_errorbar(aes(ymin = mean_eff - se, ymax = mean_eff + se),
                width = 0.18, linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = -0.5, color = COL$ref_red, linetype = "dashed",
             linewidth = 0.4, alpha = 0.6) +
  geom_text(aes(label = sprintf("%.2f", mean_eff)),
            vjust = ifelse(luad_dep$mean_eff > 0, -0.4, 1.4),
            family = FAM, size = 2.2) +
  labs(x = NULL, y = "Mean CRISPR Gene Effect",
       subtitle = sprintf("LUAD lines, n=%d · DepMap 24Q2", luad_dep$n[1])) +
  theme_pub(8)
save_panel(p_8C, file.path(OUT, "8C_bar_luad_essentiality"), 2.6, 2.4)

# Half-raincloud: violin (distribution) + quasirandom dots (raw values) +
# median crossbar.  Replaces the previous flat boxplot — much richer
# visualisation of the TCGA differential expression at the same canvas size.
cat("\n[8D-F] TCGA T vs N per lead\n")
suppressPackageStartupMessages(library(ggbeeswarm))
tn <- read_csv(file.path(DATA, "8D_tcga_TvN_long.csv"), show_col_types = FALSE)
tn_st <- read_csv(file.path(DATA, "8D_tcga_TvN_stats.csv"), show_col_types = FALSE)
panel_letter <- c(SEC61G = "8D", SRSF9 = "8E", ANGPTL4 = "8F")
for (g in LEADS) {
  sub <- tn %>% filter(gene == g)
  st <- tn_st %>% filter(gene == g)
  sub$type <- factor(sub$type, levels = c("Normal", "Tumor"))
  p <- ggplot(sub, aes(type, log2_TPM_p1, fill = type, color = type)) +
    geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
                color = NA) +
    geom_quasirandom(width = 0.22, size = 0.32, alpha = 0.55, shape = 16,
                     stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.40,
                 linewidth = 0.45, fatten = 1.6, color = "black") +
    scale_fill_manual(values = c(Normal = COL$normal, Tumor = COL$tumor),
                      guide = "none") +
    scale_color_manual(values = c(Normal = COL$normal, Tumor = COL$tumor),
                       guide = "none") +
    labs(x = NULL, y = "log2(TPM+1)",
         subtitle = sprintf("%s   log2FC = %+.2f   %s",
                            g, st$log2FC_T_minus_N, fmt_p(st$wilcoxon_p))) +
    theme_pub(8)
  save_panel(p, file.path(OUT, sprintf("%s_TvN_%s", panel_letter[[g]], g)), 2.4, 2.4)
}

cat("\n[8G] expr × CRISPR scatter (LUAD only)\n")
g8 <- read_csv(file.path(DATA, "8G_expr_vs_effect.csv"), show_col_types = FALSE)
g8$is_LUAD <- as.logical(g8$is_LUAD)
g8 <- g8 %>% filter(is_LUAD)
g8$gene <- factor(g8$gene, levels = LEADS)
ann8g <- g8 %>% group_by(gene) %>%
  summarise(rho = ifelse(n() >= 5, cor(log2_TPM_p1, gene_effect, method = "spearman"), NA_real_),
            p   = ifelse(n() >= 5, suppressWarnings(cor.test(log2_TPM_p1, gene_effect, method = "spearman"))$p.value, NA_real_),
            x   = min(log2_TPM_p1, na.rm = TRUE),
            y   = max(gene_effect, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(label = sprintf("rho = %.2f\n%s", rho, fmt_p(p)))
p_8G <- ggplot(g8, aes(log2_TPM_p1, gene_effect)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.3) +
  geom_hline(yintercept = -0.5, color = COL$ref_red, linetype = "dashed",
             linewidth = 0.3, alpha = 0.5) +
  geom_point(size = 1.0, alpha = 0.65, color = COL$normal,
             shape = 21, fill = COL$normal, stroke = 0.2) +
  geom_text(data = ann8g, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, family = FAM, size = 2.0) +
  facet_wrap(~ gene, nrow = 1, scales = "free_x") +
  labs(x = "log2(TPM+1)", y = "CRISPR Gene Effect",
       subtitle = sprintf("LUAD lines, n=%d · DepMap 24Q2",
                          length(unique(g8$ModelID[g8$gene == LEADS[1]])))) +
  theme_pub(8)
save_panel(p_8G, file.path(OUT, "8G_scatter_expr_vs_effect"), 6.0, 2.2)

# Aligned with Fig 3 KM convention:
#   geom_step (no risk-table, no survminer wrapper)
#   high = #E64B35 (red),  low = #4DBBD5 (blue)
#   p-value annotated bottom-left, legend top-right (0.70, 0.90)
cat("\n[8O/8P] KM curves\n")
km <- read_csv(file.path(DATA, "8OP_km_long.csv"), show_col_types = FALSE)
km_o <- read_csv(file.path(DATA, "8O_km_SRSF9_stats.csv"), show_col_types = FALSE)
km_p <- read_csv(file.path(DATA, "8P_km_SEC61G_stats.csv"), show_col_types = FALSE)
draw_km <- function(g, stats_row, out_stem) {
  sub <- km %>% filter(gene == g)
  sub$group <- factor(sub$group, levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ group, data = sub)
  # Build a step-function frame from survfit so we can geom_step with our
  # own theme (no survminer dependency, mirrors fig3b/d/f).
  sf <- summary(fit, times = sort(unique(c(0, sub$time))))
  sf_df <- data.frame(
    time   = sf$time,
    surv   = sf$surv,
    strata = sub("group=", "", sf$strata)
  )
  lbl_high <- sprintf("%s=high (n=%d)", g, stats_row$n_high)
  lbl_low  <- sprintf("%s=low (n=%d)",  g, stats_row$n_low)
  sf_df$strata <- factor(ifelse(sf_df$strata == "High", lbl_high, lbl_low),
                         levels = c(lbl_high, lbl_low))
  pl <- ggplot(sf_df, aes(time, surv, color = strata)) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(values = setNames(c("#E64B35", "#4DBBD5"),
                                         c(lbl_high, lbl_low)),
                       name = "Strata") +
    annotate("text", x = 0, y = 0.05,
             label = paste0("Log-rank\n", fmt_p(stats_row$logrank_p)),
             hjust = 0, vjust = 0, size = 3, family = FAM) +
    ylim(0, 1) +
    labs(x = "Time (days)", y = "OS (Overall Survival)",
         title = sprintf("%s   TCGA-LUAD n=%d (events=%d)",
                          g, nrow(sub), sum(sub$event))) +
    theme_pub(9) +
    theme(legend.position = c(0.70, 0.90),
          legend.background = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 9, face = "bold"))
  save_panel(pl, out_stem, 3.0, 2.6)
}
draw_km("SRSF9",  km_o, file.path(OUT, "8O_km_SRSF9"))
draw_km("SEC61G", km_p, file.path(OUT, "8P_km_SEC61G"))

# Switched from boxplot to violin + jitter + median crossbar; cohort label
# now displays "Takano 2024" instead of the internal "Okamura" key.
cat("\n[S11G] tumor-intrinsic ROI raincloud\n")
suppressPackageStartupMessages(library(ggbeeswarm))
ti  <- read_csv(file.path(DATA, "tumor_intrinsic_roi_long.csv"), show_col_types = FALSE)
tist<- read_csv(file.path(DATA, "tumor_intrinsic_roi_stats.csv"), show_col_types = FALSE)
ti$gene    <- factor(ti$gene, levels = LEADS)
ti$dataset <- factor(ti$dataset, levels = c("E-MTAB-13530", "Okamura"),
                                  labels = c("E-MTAB-13530", "Takano 2024"))
ti$grp     <- factor(ifelse(ti$new_roi, "ROI", "non-ROI"), levels = c("non-ROI", "ROI"))
ann_s11g <- tist %>%
  mutate(gene    = factor(gene, levels = LEADS),
         dataset = factor(dataset, levels = c("E-MTAB-13530", "Okamura"),
                                    labels = c("E-MTAB-13530", "Takano 2024")),
         label   = sprintf("Δ = %+.2f\n%s",
                            delta_new_minus_non, fmt_p(mw_p)))
ymax <- max(ti$expr, na.rm = TRUE)
p_s11g <- ggplot(ti, aes(grp, expr, fill = grp, color = grp)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
              color = NA) +
  geom_quasirandom(width = 0.22, size = 0.18, alpha = 0.32, shape = 16,
                   stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               linewidth = 0.4, fatten = 1.5, color = "black") +
  scale_fill_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI), guide = "none") +
  scale_color_manual(values = c(`non-ROI` = COL$nonROI, ROI = COL$ROI), guide = "none") +
  geom_text(data = ann_s11g, aes(x = 1.5, y = ymax * 1.05, label = label),
            inherit.aes = FALSE, family = FAM, size = 2.0) +
  facet_grid(dataset ~ gene, scales = "free_y", switch = "y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Expression (log)",
       subtitle = "Tumor-intrinsic ROI = z(Malignant)>0.5 AND z(MP3_score)>0.5") +
  theme_pub(8) +
  theme(strip.placement = "outside",
        panel.spacing = unit(0.5, "lines"))
save_panel(p_s11g, file.path(OUT, "S11G_tumor_intrinsic_ROI"), 5.5, 3.3)

cat("\nALL part1 panels DONE.\n")
