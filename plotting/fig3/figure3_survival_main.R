# Fig Treatment
# A: GSE207422 4 MP boxplot (, MP1
# B: GSE126044 MP1 boxplot
# C: GSE14814 MP2-high vs MP2-low × ACT/OBS KM
# D: GSE14814 MP4-low vs MP4-high × ACT/OBS KM (interaction p=0.020

setwd("${WORK_ROOT}/luad_figures/fig_treatment")

library(ggplot2)
library(dplyr)
library(survival)
library(ggbeeswarm)
library(patchwork)

my_font <- "sans"
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")

theme_elegant <- function(base_size=10) {
  theme_minimal(base_family=my_font, base_size=base_size) +
    theme(
      axis.text=element_text(color="black"),
      axis.line=element_line(color="black", linewidth=0.3),
      axis.ticks=element_line(color="black", linewidth=0.2),
      panel.grid=element_blank(),
      plot.title=element_text(face="bold", size=base_size+1),
      plot.subtitle=element_text(color="grey40", size=base_size-1),
      strip.text=element_text(face="bold", size=base_size),
      strip.background=element_blank()
    )
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

# A: GSE207422 4MP boxplot ( Fig 3H

cat("=== Panel A: GSE207422 ===\n")

bp207 <- read.csv("gse207422_boxplot_data.csv", stringsAsFactors=FALSE)
summ  <- read.csv("treatment_validation_summary.csv", stringsAsFactors=FALSE)

bp207 <- bp207 %>%
  filter((subset == "all" | is.na(subset)), score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(response = ifelse(group %in% c("MPR","R"), "Responder", "Non-responder"))

bp207$score <- factor(bp207$score, levels=c("MP1","MP2","MP3","MP4"))
bp207$response <- factor(bp207$response, levels=c("Responder","Non-responder"))

mp_facet_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
                     "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

p207 <- summ %>%
  filter(cohort=="GSE207422", score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(p_star = case_when(
    wilcoxon_p < 0.001 ~ "****",
    wilcoxon_p < 0.01 ~ "***",
    wilcoxon_p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

y_pos <- bp207 %>%
  group_by(score) %>%
  summarise(ymax=max(value, na.rm=TRUE), .groups="drop")
p207 <- merge(p207, y_pos, by="score")

# facet R NR
# fill
bp207$fill_key <- paste0(bp207$score, "_", bp207$response)

fill_vals <- c(
  "MP1_Responder"="#E64B35", "MP1_Non-responder"=alpha("#E64B35", 0.3),
  "MP2_Responder"="#4DBBD5", "MP2_Non-responder"=alpha("#4DBBD5", 0.3),
  "MP3_Responder"="#00A087", "MP3_Non-responder"=alpha("#00A087", 0.3),
  "MP4_Responder"="#3C5488", "MP4_Non-responder"=alpha("#3C5488", 0.3)
)

pA <- ggplot(bp207, aes(x=response, y=value, fill=fill_key)) +
  geom_boxplot(width=0.55, outlier.shape=NA, color="grey30", linewidth=0.3) +
  geom_beeswarm(aes(color=fill_key), size=1.5, alpha=0.8, cex=2.5) +
  facet_wrap(~score, nrow=1, scales="free_y", labeller=labeller(score=mp_facet_labels)) +
  scale_fill_manual(values=fill_vals, guide="none") +
  scale_color_manual(values=fill_vals, guide="none") +
  geom_text(data=p207, aes(x=1.5, y=ymax*1.08, label=p_star),
            inherit.aes=FALSE, size=5, family=my_font, fontface="bold") +
  labs(x=NULL, y="Score",
       title="GSE207422   Neoadjuvant Chemo-IO (n=24)",
       subtitle="MPR (filled) vs NMPR (transparent)") +
  theme_elegant(base_size=10) +
  theme(axis.text.x=element_text(size=8, angle=20, hjust=1))

ggsave("fig_treat_panelA.pdf", pA, width=9, height=3.5)
cat("   fig_treat_panelA.pdf\n")

# B: GSE126044 MP1 ( AUC

cat("=== Panel B: GSE126044 MP1 ===\n")

bp126 <- read.csv("gse126044_boxplot_data.csv", stringsAsFactors=FALSE)
bp126 <- bp126 %>%
  filter(score == "MP1") %>%
  mutate(response = ifelse(group == "R", "Responder", "Non-responder"))
bp126$response <- factor(bp126$response, levels=c("Responder","Non-responder"))

p126_row <- summ[summ$cohort=="GSE126044" & summ$score=="MP1", ]

yr <- range(bp126$value, na.rm=TRUE)

pB <- ggplot(bp126, aes(x=response, y=value, fill=response)) +
  geom_boxplot(width=0.5, outlier.shape=NA, color="grey30", linewidth=0.3,
               alpha=0.7) +
  geom_beeswarm(aes(color=response), size=2.5, alpha=0.8, cex=3) +
  scale_fill_manual(values=c("Responder"="#E64B35","Non-responder"=alpha("#E64B35",0.3)),
                    guide="none") +
  scale_color_manual(values=c("Responder"="#E64B35","Non-responder"=alpha("#E64B35",0.3)),
                     guide="none") +
  annotate("segment", x=1, xend=2, y=yr[2]+diff(yr)*0.08,
           yend=yr[2]+diff(yr)*0.08, linewidth=0.4) +
  annotate("text", x=1.5, y=yr[2]+diff(yr)*0.15,
           label=sprintf("%s\nAUC = %.2f", fmt_p(p126_row$wilcoxon_p), p126_row$auc),
           size=3.5, family=my_font, fontface="italic") +
  labs(x=NULL, y="MP1 Score",
       title="GSE126044   Anti-PD-1 (n=16)",
       subtitle="Independent IO validation") +
  theme_elegant(base_size=10)

ggsave("fig_treat_panelB.pdf", pB, width=3.5, height=3.5)
cat("   fig_treat_panelB.pdf\n")

# C/D: GSE14814 KM (MP2 & MP4

cat("=== Panel C/D: GSE14814 KM ===\n")

mp14 <- read.csv("gse14814_mp_scores.csv", stringsAsFactors=FALSE)
resp14 <- read.csv("gse14814_response_comparison.csv", stringsAsFactors=FALSE)

# treatment 0/1 ACT/OBS
if (all(mp14$treatment %in% c(0, 1, "0", "1"))) {
  mp14$treatment <- ifelse(as.numeric(mp14$treatment) == 1, "ACT", "OBS")
  cat("  Recoded treatment: 1->ACT, 0->OBS\n")
}

plot_strat_km <- function(data, mp_name, resp_data, panel_letter) {
  
  # median split
  med <- median(data[[mp_name]], na.rm=TRUE)
  data$mp_grp <- ifelse(data[[mp_name]] >= med, "High", "Low")
  data$strata <- paste0(mp_name, "-", data$mp_grp, " + ", data$treatment)
  
  # survfit
  data$stime <- data$os_time
  data$sevent <- data$os_event
  fit <- survfit(Surv(stime, sevent) ~ strata, data=data)
  
  sf <- summary(fit)
  km_df <- data.frame(
    time = sf$time,
    surv = sf$surv,
    upper = sf$upper,
    lower = sf$lower,
    strata = gsub("strata=", "", sf$strata)
  )
  
  strata_names <- unique(km_df$strata)
  starts <- data.frame(time=0, surv=1, upper=1, lower=1, strata=strata_names)
  km_df <- rbind(starts, km_df)
  km_df <- km_df[order(km_df$strata, km_df$time), ]
  
  inter_row <- resp_data[resp_data$score == mp_name, ]
  if (nrow(inter_row) > 1) inter_row <- inter_row[1, ]
  inter_p <- inter_row$interaction_p
  hr_high <- inter_row$hr_act_in_mp_high
  hr_low  <- inter_row$hr_act_in_mp_low
  
  strata_cols <- setNames(
    c("#B2182B","#EF8A62","#2166AC","#67A9CF"),
    c(paste0(mp_name,"-High + ACT"), paste0(mp_name,"-High + OBS"),
      paste0(mp_name,"-Low + ACT"),  paste0(mp_name,"-Low + OBS"))
  )
  
  strata_lty <- setNames(
    c("solid","dashed","solid","dashed"),
    names(strata_cols)
  )
  
  legend_labels <- setNames(
    c(paste0(mp_name,"-High + ACT (HR=",round(hr_high,2),")"),
      paste0(mp_name,"-High + OBS"),
      paste0(mp_name,"-Low + ACT (HR=",round(hr_low,2),")"),
      paste0(mp_name,"-Low + OBS")),
    names(strata_cols)
  )
  
  p <- ggplot(km_df, aes(x=time, y=surv, color=strata, linetype=strata)) +
    # CI ribbon ( ACT
    geom_ribbon(data=km_df[grepl("ACT", km_df$strata), ],
                aes(ymin=lower, ymax=upper, fill=strata),
                alpha=0.1, color=NA, show.legend=FALSE) +
    geom_step(linewidth=0.9) +
    scale_color_manual(values=strata_cols, labels=legend_labels, name=NULL) +
    scale_linetype_manual(values=strata_lty, labels=legend_labels, name=NULL) +
    scale_fill_manual(values=strata_cols, guide="none") +
    annotate("label", x=max(km_df$time)*0.55, y=0.15,
             label=sprintf("Interaction\n%s", fmt_p(inter_p)),
             size=3.2, family=my_font, fill="white", label.size=0.3,
             color="grey20") +
    labs(x="Time (years)", y="Overall Survival",
         title=sprintf("%s) GSE14814   %s × Adjuvant Chemotherapy", panel_letter, mp_name)) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_elegant(base_size=10) +
    theme(legend.position=c(0.65, 0.25),
          legend.background=element_rect(fill=alpha("white",0.9), color="grey80",
                                         linewidth=0.3),
          legend.key.width=unit(1.2, "cm"),
          legend.text=element_text(size=7.5),
          legend.margin=margin(3,5,3,5))
  
  return(p)
}

pC <- plot_strat_km(mp14, "MP2", resp14, "C")
pD <- plot_strat_km(mp14, "MP4", resp14, "D")

ggsave("fig_treat_panelC.pdf", pC, width=5.5, height=5)
ggsave("fig_treat_panelD.pdf", pD, width=5.5, height=5)
cat("   fig_treat_panelC.pdf (MP2)\n")
cat("   fig_treat_panelD.pdf (MP4)\n")


cat("=== Combining ===\n")

# A A A A A A A A B B B
# C C C C C D D D D D D

layout <- "
AAAAAAAABBB
AAAAAAAABBB
CCCCCDDDDDD
CCCCCDDDDDD
CCCCCDDDDDD
"

p_main <- pA + pB + pC + pD + plot_layout(design=layout) +
  plot_annotation(
    tag_levels="A",
    theme=theme(plot.tag=element_text(face="bold", size=14, family=my_font))
  )

ggsave("fig_treatment_main.pdf", p_main, width=13, height=10)
cat("   fig_treatment_main.pdf\n")

# Supplementary: ns (GSE135222 + ACT boxplot

cat("\n=== Supplementary ===\n")

# GSE135222 4 MP
bp135 <- read.csv("gse135222_boxplot_data.csv", stringsAsFactors=FALSE) %>%
  filter(score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(response = ifelse(group %in% c("R","DCB"), "Responder", "Non-responder"))

bp135$score <- factor(bp135$score, levels=c("MP1","MP2","MP3","MP4"))
bp135$response <- factor(bp135$response, levels=c("Responder","Non-responder"))

p135_vals <- summ %>%
  filter(cohort=="GSE135222", score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(p_star = ifelse(wilcoxon_p < 0.05, "*", "ns"))

y135 <- bp135 %>%
  group_by(score) %>%
  summarise(ymax=max(value,na.rm=TRUE), .groups="drop")
p135_vals <- merge(p135_vals, y135, by="score")

p_supp135 <- ggplot(bp135, aes(x=response, y=value, fill=response)) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.6, color="grey30") +
  geom_beeswarm(size=1.2, alpha=0.7, cex=2, color="grey30") +
  facet_wrap(~score, nrow=1, scales="free_y") +
  scale_fill_manual(values=c("Responder"="#ED0000","Non-responder"="#00468B"), guide="none") +
  geom_text(data=p135_vals, aes(x=1.5, y=ymax*1.08, label=p_star),
            inherit.aes=FALSE, size=5, family=my_font) +
  labs(x=NULL, y="Score",
       title="GSE135222   Anti-PD-1/PD-L1 (n=27)",
       subtitle="DCB vs non-DCB") +
  theme_elegant(base_size=10)

ggsave("fig_treat_supp_gse135222.pdf", p_supp135, width=9, height=3.5)
cat("   fig_treat_supp_gse135222.pdf\n")

cat("\n All treatment panels saved.\n")
cat("  Main: fig_treatment_main.pdf (A=207422 boxplot, B=126044 MP1, C=MP2 KM, D=MP4 KM)\n")
cat("  Supp: fig_treat_supp_gse135222.pdf\n")