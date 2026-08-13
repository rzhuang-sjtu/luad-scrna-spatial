#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  source("${PROJECT_ROOT}/plotting/fig8/fig8_theme.R")
})
DATA <- "${WORK_ROOT}/luad_figures/fig8/v2_500/data"
OUT  <- "${WORK_ROOT}/luad_figures/fig8/v2_500"

km <- read_csv(file.path(DATA, "8OPQ_km_long.csv"), show_col_types = FALSE)
st <- read_csv(file.path(DATA, "8O_km_ANGPTL4_stats.csv"), show_col_types = FALSE)

sub <- km %>% filter(gene == "ANGPTL4")
lbl_hi <- sprintf("ANGPTL4=high (n=%d)", st$n_high)
lbl_lo <- sprintf("ANGPTL4=low (n=%d)",  st$n_low)
sub$strata <- factor(ifelse(sub$group == "High", lbl_hi, lbl_lo),
                     levels = c(lbl_hi, lbl_lo))

p <- ggplot(sub, aes(time, surv_prob, color = strata)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = setNames(c("#E64B35","#4DBBD5"),
                                        c(lbl_hi, lbl_lo)),
                     name = "Strata") +
  annotate("text", x = 0, y = 0.05,
           label = paste0("Log-rank ", fmt_p(st$logrank_p)),
           hjust = 0, vjust = 0, size = 3, family = FAM) +
  ylim(0, 1) +
  labs(x = "Time (days)", y = "OS (Overall Survival)",
       title = sprintf("ANGPTL4   TCGA-LUAD n=%d (events=%d)",
                        st$n_high + st$n_low,
                        st$events_high + st$events_low)) +
  theme_pub(9) +
  theme(legend.position  = c(0.70, 0.90),
        legend.background = element_blank(),
        legend.key.size   = unit(0.35, "cm"),
        legend.text       = element_text(size = 7),
        legend.title      = element_text(size = 8, face = "bold"),
        plot.title        = element_text(size = 9, face = "bold"))

save_panel(p, file.path(OUT, "8O_km_ANGPTL4"), 3.0, 2.6)
cat("8O ANGPTL4 KM (median split) saved\n")
