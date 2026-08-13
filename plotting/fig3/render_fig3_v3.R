#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(dplyr); library(grid)
  library(showtext); library(ggrepel); library(ggnewscale); library(patchwork)
})
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p)
  showtext_auto(); showtext_opts(dpi = 300)
  FAM <- "Arial"
} else FAM <- "sans"

setwd("${WORK_ROOT}/luad_figures/fig3")
mp_colors <- c(MP1 = "#E64B35", MP2 = "#4DBBD5", MP3 = "#00A087", MP4 = "#3C5488")
mp_lbl <- c(MP1 = "MP1: Stress/AP-1", MP2 = "MP2: Proliferative",
            MP3 = "MP3: EMT/IFN",   MP4 = "MP4: AT2-like")
theme_pub <- function(base_size = 8) {
  theme_classic(base_family = FAM, base_size = base_size) +
    theme(axis.text  = element_text(color = "black"),
          axis.line  = element_line(linewidth = 0.4, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"))
}

# Fig 3B  trajectory: middle ground arrows + corner-only axis (blank big axes)
pt_data <- as.data.frame(fread("../fig2/pseudotime_umap.csv.gz"))
seg <- as.data.frame(fread("../fig2/monocle3_graph_segments_oriented.csv"))
mark <- as.data.frame(fread("../fig2/monocle3_graph_root_tip.csv"))

set.seed(42)
n_plot <- min(nrow(pt_data), 20000)
pt_sub <- pt_data[sample(nrow(pt_data), n_plot), ]
pt_sub <- pt_sub[pt_sub$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
pt_sub <- pt_sub[sample(nrow(pt_sub)), ]

# 14 arrows: take longest segments per pt bin (14 bins)
seg$len <- sqrt((seg$x_end - seg$x_start)^2 + (seg$y_end - seg$y_start)^2)
seg$bin <- cut(seg$pt_mean, breaks = 14, include.lowest = TRUE)
arrow_seg <- seg %>% group_by(bin) %>% slice_max(len, n = 1, with_ties = FALSE) %>% ungroup()

pt_lim <- range(c(seg$pt_start, seg$pt_end), na.rm = TRUE)
arr <- arrow(length = unit(1.8, "mm"), type = "closed", ends = "last")

xr <- range(pt_sub$UMAP1, na.rm = TRUE); yr <- range(pt_sub$UMAP2, na.rm = TRUE)
inset <- 0.05; frac <- 0.16
ax0 <- xr[1] + inset*diff(xr); ay0 <- yr[1] + inset*diff(yr)
ax1 <- ax0 + frac*diff(xr);    ay1 <- ay0 + frac*diff(yr)

p3b <- ggplot() +
  geom_point(data = pt_sub,
             aes(UMAP1, UMAP2, color = dominant_MP),
             size = 0.55, alpha = 0.45, stroke = 0, shape = 16) +
  scale_color_manual(values = mp_colors, labels = mp_lbl,
                     name = "Dominant MP",
                     guide = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  ggnewscale::new_scale_color() +
  geom_segment(data = seg,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                   color = pt_mean),
               linewidth = 0.85, lineend = "round") +
  geom_segment(data = arrow_seg,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                   color = pt_mean),
               linewidth = 1.05, arrow = arr,
               lineend = "round", linejoin = "mitre") +
  scale_color_viridis_c(option = "inferno", limits = pt_lim,
                        name = "Pseudotime\n(early -> late)",
                        guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                               barheight = unit(20, "mm"),
                                               frame.colour = "black",
                                               frame.linewidth = 0.3)) +
  geom_point(data = mark, aes(x, y, shape = kind),
             size = 3.6, fill = "white", color = "black", stroke = 0.6) +
  geom_text_repel(data = mark, aes(x, y, label = label),
                  family = FAM, size = 2.6, fontface = "bold",
                  segment.size = 0.3, box.padding = 0.4, point.padding = 0.4) +
  scale_shape_manual(values = c(root = 23, tip = 21), guide = "none") +
  # Corner Component 1 / 2 small arrows (drawn last so they sit on top)
  annotate("segment", x = ax0, xend = ax1, y = ay0, yend = ay0,
           arrow = arrow(length = unit(1.2,"mm"), type = "closed"),
           linewidth = 0.3, color = "black") +
  annotate("segment", x = ax0, xend = ax0, y = ay0, yend = ay1,
           arrow = arrow(length = unit(1.2,"mm"), type = "closed"),
           linewidth = 0.3, color = "black") +
  annotate("text", x = (ax0+ax1)/2, y = ay0, label = "Component 1",
           vjust = 2.4, size = 2.4, family = FAM) +
  annotate("text", x = ax0, y = (ay0+ay1)/2, label = "Component 2",
           angle = 90, vjust = -1.4, size = 2.4, family = FAM) +
  labs(x = NULL, y = NULL,
       title = "Malignant cell trajectory (Monocle3 principal graph)") +
  coord_equal() +
  theme_void(base_family = FAM, base_size = 9) +
  theme(plot.title       = element_text(size = 9, face = "bold",
                                        margin = margin(b = 4)),
        legend.position  = "right",
        legend.title     = element_text(size = 7, face = "bold"),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"),
        legend.box       = "vertical",
        legend.spacing.y = unit(2, "mm"))

ggsave("fig3b_trajectory.pdf", p3b, width = 5.6, height = 4.5)
ggsave("fig3b_trajectory.png", p3b, width = 5.6, height = 4.5, dpi = 300)
cat("Fig 3B: 14 arrows + corner-only axes (big axes blanked)\n")

# Fig 3C  — minimal: faded grey + bold TF/Surface + halo labels
gs <- read.csv("geneswitches_results.csv", stringsAsFactors = FALSE)
tf_zscore <- read.csv("tf_activity_mp_zscore.csv", stringsAsFactors = FALSE)
known_tfs <- tf_zscore[[1]]
surface_proteins <- c("EPCAM","CD44","CD24","CEACAM5","CEACAM6","MUC1",
                      "EGFR","ERBB2","ERBB3","MET","AXL","PDGFRA",
                      "ITGA6","ITGB1","ITGB4","ICAM1","ALCAM","CXCR4",
                      "SDC1","SDC2","SDC4","TNFRSF10A","TNFRSF10B",
                      "LY6K","LY6E","NECTIN2","NECTIN4","CD274","PDCD1LG2",
                      "THY1","ENG","PECAM1","CDH1","CDH2","VIM",
                      "FN1","SPARC","LAMB3","COL1A1","COL3A1",
                      "CXCL1","CXCL2","CXCL8","CCL2","CCL5",
                      "F2","FGA","FGB","C3","CFH","VTN",
                      "APOC2","APOC3","APOB","APOH","CPS1","ARG1",
                      "SLC2A2","HNF4A")
gs$gene_type <- "Other"
gs$gene_type[gs$gene %in% known_tfs] <- "TF"
gs$gene_type[gs$gene %in% surface_proteins] <- "Surface"
gs$gene_type <- factor(gs$gene_type, levels = c("Other","TF","Surface"))
gs$switch_pct <- 100 * (rank(gs$switch_pseudotime_rank, ties.method = "average") /
                          length(gs$switch_pseudotime_rank))

# Label rule: every TF + every Surface protein + top-N "Other" by R²
n_top_other <- 15
top_other <- gs %>% filter(gene_type == "Other") %>%
  slice_max(mcfadden_R2, n = n_top_other) %>% pull(gene)
label_set <- unique(c(gs$gene[gs$gene_type != "Other"], top_other))
gs$label <- ifelse(gs$gene %in% label_set, gs$gene, "")

# Reorder so labeled / colored points draw on top of the unlabeled cloud
gs$z_order <- with(gs,
  ifelse(gene_type != "Other", 3,
  ifelse(label != "",          2, 1)))
gs <- gs[order(gs$z_order), ]

dat_other_unlab <- gs %>% filter(gene_type == "Other", label == "")
dat_other_lab   <- gs %>% filter(gene_type == "Other", label != "")
dat_hl          <- gs %>% filter(gene_type != "Other")

type_fill <- c(TF = "#00A087", Surface = "#E64B35")
type_col  <- c(TF = "#00543E", Surface = "#A4291C")

# Map every point through scale_fill_manual so all 3 tiers appear in legend
gs$plot_type <- factor(
  ifelse(gs$gene_type == "TF",      "TF",
  ifelse(gs$gene_type == "Surface", "Surface", "Other")),
  levels = c("Other","TF","Surface"))

dat_other_unlab <- gs %>% filter(plot_type == "Other", label == "")
dat_other_lab   <- gs %>% filter(plot_type == "Other", label != "")
dat_hl          <- gs %>% filter(plot_type != "Other")

tier_fill <- c(Other = "grey75", TF = "#00A087", Surface = "#E64B35")
tier_col  <- c(Other = "grey25", TF = "#00543E", Surface = "#A4291C")
tier_lbl  <- c(Other = "Other (AT2 / MHC-II / etc.)",
               TF    = "Transcription factors",
               Surface = "Surface proteins")

p3c <- ggplot() +
  # 1) unlabeled "Other" genes — light grey, still clearly visible
  geom_point(data = dat_other_unlab,
             aes(x = switch_pct, y = mcfadden_R2, fill = plot_type, color = plot_type),
             shape = 21, size = 1.4, stroke = 0.25, alpha = 0.6,
             show.legend = FALSE) +
  # 2) labeled "Other" genes — darker grey, slightly larger (drives legend entry)
  geom_point(data = dat_other_lab,
             aes(x = switch_pct, y = mcfadden_R2, fill = plot_type, color = plot_type),
             shape = 21, size = 2.0, stroke = 0.4, alpha = 0.95) +
  # 3) TF / Surface highlights
  geom_point(data = dat_hl,
             aes(x = switch_pct, y = mcfadden_R2, fill = plot_type, color = plot_type),
             shape = 21, size = 2.3, stroke = 0.5, alpha = 0.95) +
  scale_fill_manual(values = tier_fill, name = NULL, labels = tier_lbl,
                    breaks = c("TF","Surface","Other")) +
  scale_color_manual(values = tier_col, guide = "none") +
  # 4) labels with halo — italic, anchored to every colored / top-Other dot
  ggrepel::geom_text_repel(
    data = gs %>% filter(label != ""),
    aes(x = switch_pct, y = mcfadden_R2, label = label),
    size = 2.5, family = FAM, fontface = "bold.italic", color = "black",
    bg.color = "white", bg.r = 0.20,
    box.padding = 0.35, point.padding = 0.25,
    segment.size = 0.22, segment.alpha = 0.5,
    min.segment.length = 0.05, max.overlaps = Inf
  ) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%","25%","50%","75%","100%"),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Pseudo-timeline (percentile)",
       y = expression("McFadden " * R^2),
       title = "GeneSwitches analysis") +
  theme_pub(8.5) +
  theme(plot.title       = element_text(size = 9, face = "bold"),
        legend.position  = c(0.86, 0.92),
        legend.background = element_blank(),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"),
        panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3))

ggsave("fig3c_geneswitches.pdf", p3c, width = 5.4, height = 3.8)
ggsave("fig3c_geneswitches.png", p3c, width = 5.4, height = 3.8, dpi = 300)
cat("Fig 3C: 3-tier scatter (TF / Surface / Other) — no LOESS line\n")

# Fig 3A+C combined Cox forest — replace fig3a_univariate_cox + fig3c_multivariate
cox_uv <- read.csv("tcga_luad_mp_cox_univariate.csv", stringsAsFactors = FALSE)
cox_mv <- read.csv("tcga_luad_mp_cox_multivariate.csv", stringsAsFactors = FALSE)
cox_uv$model <- "Univariate"
cox_uv$row_name <- cox_uv$MP
cox_mv$model <- "Multivariate"
cox_mv$row_name <- dplyr::recode(cox_mv$covariate,
                                  "age" = "Age", "stage_num" = "Stage",
                                  "MP1_score" = "MP1", "MP2_score" = "MP2",
                                  "MP3_score" = "MP3", "MP4_score" = "MP4",
                                  .default = cox_mv$covariate)
cox_uv <- cox_uv[, c("row_name","HR","HR_lo","HR_hi","p","model")]
cox_mv <- cox_mv[, c("row_name","HR","HR_lo","HR_hi","p","model")]
cox_combined <- rbind(cox_uv, cox_mv)
cox_combined$model <- factor(cox_combined$model,
                             levels = c("Univariate","Multivariate"))
# y order: Univariate top (MP1-4 reversed), Multivariate bottom (Age, Stage, MP1-4)
uv_order <- c("MP4","MP3","MP2","MP1")
mv_order <- c("MP4","MP3","MP2","MP1","Stage","Age")
cox_combined$row_name <- factor(cox_combined$row_name,
                                levels = unique(c(uv_order, mv_order)))
cox_combined$sig <- ifelse(cox_combined$p < 0.05, "sig", "ns")
cox_combined$lbl <- sprintf("%.2f (%.2f-%.2f)",
                             cox_combined$HR, cox_combined$HR_lo, cox_combined$HR_hi)

# log-symmetric x-axis covering both models
log_span <- max(abs(log10(c(min(cox_combined$HR_lo, na.rm = TRUE),
                              max(cox_combined$HR_hi, na.rm = TRUE)))), na.rm = TRUE)
xlim_lo <- 10^(-log_span * 1.05); xlim_hi <- 10^(log_span * 1.05)
label_x <- xlim_hi * 1.15

# colour MP rows by mp_colors; non-MP rows in grey
cox_combined$mp_for_color <- ifelse(cox_combined$row_name %in% c("MP1","MP2","MP3","MP4"),
                                     as.character(cox_combined$row_name), "Other")
clr <- c(mp_colors, Other = "grey45")

p3ac <- ggplot(cox_combined, aes(x = HR, y = row_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi, color = mp_for_color),
                 height = 0.22, linewidth = 0.65) +
  geom_point(aes(fill = mp_for_color, color = mp_for_color, shape = sig),
             size = 2.4, stroke = 0.5) +
  geom_text(aes(x = label_x, label = lbl), hjust = 0,
            size = 2.0, family = FAM) +
  scale_fill_manual(values = clr, guide = "none") +
  scale_color_manual(values = clr, guide = "none") +
  scale_shape_manual(values = c(sig = 21, ns = 1), guide = "none") +
  scale_x_log10(limits = c(xlim_lo, xlim_hi),
                breaks = c(0.25, 0.5, 1, 2, 4)) +
  facet_grid(model ~ ., scales = "free_y", space = "free_y", switch = "y") +
  coord_cartesian(clip = "off") +
  labs(x = "Hazard ratio (log scale)", y = NULL,
       title = "TCGA-LUAD Cox regression",
       subtitle = "Univariate (top) | Multivariate adj. for age + stage (bottom)") +
  theme_pub(8) +
  theme(plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 7, color = "grey25"),
        strip.background = element_rect(fill = "grey94", color = NA),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 8,
                                         color = "black"),
        strip.placement  = "outside",
        plot.margin      = margin(2, 60, 2, 4),
        panel.spacing.y  = unit(2, "mm"),
        axis.text.y      = element_text(size = 7))

ggsave("fig3ac_cox_combined.pdf", p3ac, width = 5.0, height = 3.6)
ggsave("fig3ac_cox_combined.png", p3ac, width = 5.0, height = 3.6, dpi = 300)
cat("Fig 3A+C combined Cox forest saved (5.0 x 3.6 in)\n")

# Fig 3H  even more compact
chemo <- read.csv("tcga_mp3_chemokine_correlation.csv", stringsAsFactors = FALSE)
chemo$present <- tolower(as.character(chemo$present)) %in% c("true","1","t")
chemo <- chemo[chemo$present, ]
chemo$sig <- ifelse(chemo$p < 0.05, "sig", "ns")
chemo$chemokine <- factor(chemo$chemokine,
                          levels = chemo$chemokine[order(chemo$spearman_rho)])

p3h <- ggplot(chemo, aes(x = spearman_rho, y = chemokine, fill = sig)) +
  geom_col(width = 0.82, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(sig = "#00A087", ns = "grey70"), guide = "none") +
  geom_text(aes(x = spearman_rho, label = sprintf("%.2f", spearman_rho)),
            hjust = -0.18, size = 1.7, family = FAM) +
  scale_x_continuous(limits = c(0, max(chemo$spearman_rho, na.rm = TRUE) * 1.18),
                     expand = expansion(mult = c(0, 0.02)),
                     breaks = c(0, 0.2, 0.4)) +
  labs(x = "Spearman rho (TCGA-LUAD)", y = NULL,
       title = "MP3 vs chemokine correlation") +
  theme_pub(7) +
  theme(axis.text.y  = element_text(face = "italic", size = 5.5),
        axis.text.x  = element_text(size = 5.5),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title   = element_text(size = 7, face = "bold",
                                    margin = margin(b = 1)),
        plot.margin  = margin(2, 4, 2, 1))

ggsave("fig3h_chemokine_corr.pdf", p3h, width = 2.8, height = 2.0)
ggsave("fig3h_chemokine_corr.png", p3h, width = 2.8, height = 2.0, dpi = 300)
cat("Fig 3H ultra-compact saved (2.8 x 2.0 in)\n")
