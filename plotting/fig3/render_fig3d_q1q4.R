#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork)
  library(showtext); library(grid)
})
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p)
  showtext_auto(); showtext_opts(dpi = 300)
  FAM <- "Arial"
} else FAM <- "sans"

setwd("${WORK_ROOT}/luad_figures/fig3")
theme_pub <- function(base = 8) {
  theme_classic(base_family = FAM, base_size = base) +
    theme(axis.text  = element_text(color = "black"),
          axis.line  = element_line(linewidth = 0.4, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"))
}
mp_long <- c(MP1 = "MP1: Stress/AP-1", MP2 = "MP2: Proliferative",
             MP3 = "MP3: EMT/IFN",     MP4 = "MP4: AT2-like")
fmt_p <- function(p) {
  if (p < 1e-3) sprintf("p = %.1e", p)
  else if (p < 1e-2) sprintf("p = %.3f", p)
  else sprintf("p = %.3f", p)
}

curves <- as.data.frame(fread("tcga_luad_mp_q1q4_km_curves.csv.gz"))
stats  <- as.data.frame(fread("tcga_luad_mp_q1q4_km_logrank.csv"))
mvr    <- as.data.frame(fread("mp_full_panel_survival_summary.csv"))

# Panel builder
make_panel <- function(mp) {
  d  <- curves[curves$MP == mp, ]
  st <- stats[stats$MP == mp, ]
  mv <- mvr[mvr$MP == mp, ]
  d$strata <- factor(ifelse(d$group == "High",
                             paste0(mp, "=Q4 (n=", st$n_high, ")"),
                             paste0(mp, "=Q1 (n=", st$n_low,  ")")),
                     levels = c(paste0(mp, "=Q4 (n=", st$n_high, ")"),
                                paste0(mp, "=Q1 (n=", st$n_low,  ")")))
  ttl <- mp_long[[mp]]
  sub <- sprintf("Q1 vs Q4 log-rank %s | Cox MV HR=%.2f (%.2f-%.2f), %s",
                 fmt_p(st$log_rank_p),
                 mv$Cox_mv_HR, mv$Cox_mv_HR_lo, mv$Cox_mv_HR_hi,
                 fmt_p(mv$Cox_mv_p))
  ggplot(d, aes(time, surv_prob, color = strata)) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(values = c("#E64B35","#4DBBD5"), name = NULL) +
    annotate("text", x = 0, y = 0.05,
             label = paste0("log-rank ", fmt_p(st$log_rank_p)),
             hjust = 0, vjust = 0, size = 2.4, family = FAM) +
    ylim(0, 1) +
    labs(x = "Time (days)", y = "Overall survival",
         title = ttl, subtitle = sub) +
    theme_pub(8) +
    theme(plot.title       = element_text(size = 9, face = "bold"),
          plot.subtitle    = element_text(size = 6.5, color = "grey25",
                                          margin = margin(b = 2)),
          legend.position  = c(0.70, 0.92),
          legend.background = element_blank(),
          legend.text      = element_text(size = 6.5),
          legend.key.size  = unit(3, "mm"))
}

panels <- lapply(c("MP1","MP2","MP3","MP4"), make_panel)
fig3d <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  plot_annotation(
    title = "TCGA-LUAD overall survival - MP top vs bottom quartile (Q1 vs Q4)",
    subtitle = "Predefined cutoffs (no scanning); log-rank shown on curves; multivariate Cox HR adjusted for age, stage and other MPs.",
    theme = theme(plot.title = element_text(size = 10, face = "bold", family = FAM),
                  plot.subtitle = element_text(size = 7, color = "grey25",
                                                family = FAM, margin = margin(b = 4))))

ggsave("fig3d_km_q1q4.pdf", fig3d, width = 8, height = 6.4)
ggsave("fig3d_km_q1q4.png", fig3d, width = 8, height = 6.4, dpi = 300)
cat("Fig 3D Q1/Q4 KM grid saved\n")
