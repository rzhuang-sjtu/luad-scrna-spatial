#!/usr/bin/env Rscript
# Fig 8 v2 — treatment panels: 8M (GSE207422 boxplot), 8N (GSE126044 volcano),
# S11D/E/F (GSE135222 boxplots).

suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
  library(ggrepel); library(ggbeeswarm)
})

DATA <- "${PROJECT_ROOT}/results/fig8_plot_data/v2_500"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"
LEADS <- c("SEC61G", "SRSF9", "ANGPTL4")

cat("\n[8M] GSE207422 boxplot\n")
m  <- read_csv(file.path(DATA, "8M_GSE207422_long.csv"), show_col_types = FALSE)
ms <- read_csv(file.path(DATA, "8M_GSE207422_stats.csv"), show_col_types = FALSE)
m$gene     <- factor(m$gene, levels = LEADS)
m$response <- factor(m$response, levels = c("NMPR", "MPR"))
ann_m <- ms %>% mutate(gene = factor(gene, levels = LEADS),
                       label = sprintf("AUC = %.2f\n%s",
                                       auc_pos_vs_neg, fmt_p(p)))
ymax_m <- max(m$expr, na.rm = TRUE)
p_8M <- ggplot(m, aes(response, expr, fill = response, color = response)) +
  geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
              color = NA) +
  geom_quasirandom(width = 0.22, size = 0.55, alpha = 0.75, shape = 16,
                   stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               linewidth = 0.45, fatten = 1.5, color = "black") +
  scale_fill_manual(values = c(NMPR = COL$NMPR, MPR = COL$MPR), guide = "none") +
  scale_color_manual(values = c(NMPR = COL$NMPR, MPR = COL$MPR), guide = "none") +
  geom_text(data = ann_m, aes(x = 1.5, y = ymax_m * 1.04, label = label),
            inherit.aes = FALSE, family = FAM, size = 2.0) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "log2(TPM+1)",
       subtitle = sprintf("GSE207422 — MPR (n=%d) vs NMPR (n=%d) · neoadjuvant chemo-IO",
                          ms$n_pos[1], ms$n_neg[1])) +
  theme_pub(8)
save_panel(p_8M, file.path(OUT, "8M_GSE207422_3genes"), 6.0, 2.4)

cat("\n[8N] GSE126044 volcano\n")
v <- read_csv(file.path(DATA, "8N_GSE126044_volcano.csv"), show_col_types = FALSE)
v$is_lead <- as.logical(v$is_lead)
v_sig  <- v %>% filter(p < 0.05, !is_lead)
v_lead <- v %>% filter(is_lead)
p_max <- max(v$nlog10p, na.rm = TRUE)
p_8N <- ggplot(v, aes(log2FC_R_minus_NR, nlog10p)) +
  geom_point(data = subset(v, !is_lead), color = "grey75", size = 0.4, alpha = 0.45) +
  geom_point(data = v_sig, color = "grey50", size = 0.4, alpha = 0.55) +
  geom_point(data = v_lead, aes(color = gene), size = 1.6, shape = 21,
             fill = "white", stroke = 0.7) +
  geom_text_repel(data = v_lead, aes(label = gene, color = gene),
                  family = FAM, size = 2.4, fontface = "bold",
                  box.padding = 0.4, point.padding = 0.2,
                  segment.size = 0.3, max.overlaps = Inf) +
  geom_vline(xintercept = c(-1, 0, 1), linetype = c("dashed","dotted","dashed"),
             color = c(COL$ref_red,"black",COL$ref_red), linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = COL$ref_red, linewidth = 0.3) +
  scale_color_manual(values = c(SEC61G = COL$LUAD, SRSF9 = COL$venn_neu, ANGPTL4 = COL$Macro_SPP1),
                     guide = "none") +
  labs(x = "log2 FC (R − NR)", y = "-log10(p)",
       subtitle = sprintf("GSE126044 — anti-PD-1 R vs NR · %d genes · leads in color",
                          nrow(v))) +
  theme_pub(8)
save_panel(p_8N, file.path(OUT, "8N_GSE126044_volcano"), 4.5, 3.0)

cat("\n[S11D/E/F] GSE135222 single-gene boxplots\n")
g3  <- read_csv(file.path(DATA, "S11D_GSE135222_long.csv"),  show_col_types = FALSE)
g3s <- read_csv(file.path(DATA, "S11D_GSE135222_stats.csv"), show_col_types = FALSE)
g3$response <- factor(g3$response, levels = c("NR", "R"))
panel_letter <- c(SEC61G = "S11D", SRSF9 = "S11E", ANGPTL4 = "S11F")
for (g in LEADS) {
  sub <- g3 %>% filter(gene == g)
  st  <- g3s %>% filter(gene == g)
  if (!nrow(sub)) { cat(sprintf("  WARN: no data for %s\n", g)); next }
  pl <- ggplot(sub, aes(response, expr, fill = response, color = response)) +
    geom_violin(width = 0.78, trim = FALSE, alpha = 0.35, linewidth = 0.3,
                color = NA) +
    geom_quasirandom(width = 0.22, size = 0.6, alpha = 0.8, shape = 16,
                     stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.42,
                 linewidth = 0.45, fatten = 1.5, color = "black") +
    scale_fill_manual(values = c(NR = COL$NR, R = COL$R), guide = "none") +
    scale_color_manual(values = c(NR = COL$NR, R = COL$R), guide = "none") +
    labs(x = NULL, y = "log2(TPM+1)",
         subtitle = sprintf("%s   GSE135222 R vs NR   log2FC = %+.2f   %s",
                            g, st$log2FC_pos_minus_neg, fmt_p(st$p))) +
    theme_pub(8)
  save_panel(pl, file.path(OUT, sprintf("%s_GSE135222_%s", panel_letter[[g]], g)), 2.6, 2.4)
}

cat("\nALL part3 panels DONE.\n")
