# Supplementary: Treatment Validation
# Table S_treat: × MP
# Fig S_treat: boxplot + KM

setwd("${WORK_ROOT}/luad_figures/fig_treatment")

library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(ggbeeswarm)
library(patchwork)
library(data.table)

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
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  if (p < 0.01) return(sprintf("%.3f", p))
  return(sprintf("%.2f", p))
}

# Supplementary Table

cat("=== Supplementary Table ===\n")

summ <- read.csv("treatment_validation_summary.csv", stringsAsFactors=FALSE)

supp_table <- summ %>%
  filter(score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(
    p_formatted = sapply(wilcoxon_p, fmt_p),
    auc_formatted = ifelse(is.na(auc), "NA", sprintf("%.3f", auc)),
    interaction_p_formatted = ifelse(is.na(interaction_p), "NA", sapply(interaction_p, fmt_p)),
    hr_high_formatted = ifelse(is.na(hr_act_in_high), "NA", sprintf("%.2f", hr_act_in_high)),
    hr_low_formatted = ifelse(is.na(hr_act_in_low), "NA", sprintf("%.2f", hr_act_in_low)),
    significance = case_when(
      !is.na(wilcoxon_p) & wilcoxon_p < 0.05 ~ "Significant",
      !is.na(interaction_p) & interaction_p < 0.05 ~ "Significant",
      !is.na(interaction_p) & interaction_p < 0.1 ~ "Trend",
      TRUE ~ "NS"
    )
  ) %>%
  select(cohort, treatment, analysis_type, score, n,
         median_R, median_NR, delta,
         p_formatted, auc_formatted,
         interaction_p_formatted, hr_high_formatted, hr_low_formatted,
         significance)

colnames(supp_table) <- c("Cohort","Treatment","Analysis","MP","n",
                          "Median_R","Median_NR","Delta",
                          "Wilcoxon_p","AUC",
                          "Interaction_p","HR_ACT_High","HR_ACT_Low",
                          "Significance")

write.csv(supp_table, "supp_table_treatment_full.csv", row.names=FALSE)
cat("   supp_table_treatment_full.csv\n")

# Supp Fig A: IO × 4 MP ( ns GSE135222

cat("\n=== Supp Fig: IO boxplots full ===\n")

load_bp <- function(file, cohort_name, r_labels, nr_labels) {
  d <- read.csv(file, stringsAsFactors=FALSE)
  d <- d %>%
    filter(score %in% c("MP1","MP2","MP3","MP4")) %>%
    mutate(response = ifelse(group %in% r_labels, "R", "NR"),
           cohort = cohort_name)
  # subset
  if ("subset" %in% colnames(d)) d <- d %>% filter(subset == "all" | is.na(subset))
  d %>% select(score, value, response, cohort)
}

bp_207 <- load_bp("gse207422_boxplot_data.csv", 
                  "GSE207422 (chemo-IO, n=24)", c("MPR","R"), c("NMPR","NR"))
bp_126 <- load_bp("gse126044_boxplot_data.csv",
                  "GSE126044 (anti-PD-1, n=16)", "R", "NR")
bp_135 <- load_bp("gse135222_boxplot_data.csv",
                  "GSE135222 (anti-PD-1/L1, n=27)", c("R","DCB"), c("NR","non-DCB"))

io_all <- bind_rows(bp_207, bp_126, bp_135)
io_all$score <- factor(io_all$score, levels=c("MP1","MP2","MP3","MP4"))
io_all$response <- factor(io_all$response, levels=c("R","NR"))
io_all$cohort <- factor(io_all$cohort, 
                        levels=c("GSE207422 (chemo-IO, n=24)",
                                 "GSE126044 (anti-PD-1, n=16)",
                                 "GSE135222 (anti-PD-1/L1, n=27)"))

io_p <- summ %>%
  filter(analysis_type == "IO", score %in% c("MP1","MP2","MP3","MP4")) %>%
  mutate(
    p_star = case_when(
      wilcoxon_p < 0.001 ~ "****",
      wilcoxon_p < 0.01 ~ "***",
      wilcoxon_p < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    auc_label = sprintf("AUC=%.2f", auc),
    cohort = case_when(
      cohort == "GSE207422" ~ "GSE207422 (chemo-IO, n=24)",
      cohort == "GSE126044" ~ "GSE126044 (anti-PD-1, n=16)",
      cohort == "GSE135222" ~ "GSE135222 (anti-PD-1/L1, n=27)"
    )
  )

y_io <- io_all %>%
  group_by(cohort, score) %>%
  summarise(ymax = max(value, na.rm=TRUE), .groups="drop")
io_p <- merge(io_p, y_io, by=c("cohort","score"))
io_p$cohort <- factor(io_p$cohort, levels=levels(io_all$cohort))

p_supp_io <- ggplot(io_all, aes(x=response, y=value, fill=response)) +
  geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.6, color="grey30", linewidth=0.3) +
  geom_beeswarm(size=1.2, alpha=0.7, cex=2, color="grey30") +
  facet_grid(cohort ~ score, scales="free_y") +
  scale_fill_manual(values=c("R"="#E64B35","NR"="#4DBBD5"), guide="none") +
  # p + AUC
  geom_text(data=io_p, aes(x=1.5, y=ymax*1.05, label=p_star),
            inherit.aes=FALSE, size=4, family=my_font, fontface="bold") +
  geom_text(data=io_p, aes(x=1.5, y=ymax*0.97, label=auc_label),
            inherit.aes=FALSE, size=2.5, family=my_font, color="grey40") +
  labs(x=NULL, y="MP Score",
       title="Supplementary: IO Response   All Cohorts × All MPs",
       subtitle="R (red) vs NR (blue)") +
  theme_elegant(base_size=9) +
  theme(strip.text.y=element_text(angle=0, hjust=0, size=7),
        axis.text.x=element_text(size=8))

ggsave("supp_fig_io_full.pdf", p_supp_io, width=9, height=8)
cat("   supp_fig_io_full.pdf\n")

# Supp Fig B: ACT 4 MP × 2 KM

cat("\n=== Supp Fig: ACT stratified KM full ===\n")

mp14 <- read.csv("gse14814_mp_scores.csv", stringsAsFactors=FALSE)
resp14 <- read.csv("gse14814_response_comparison.csv", stringsAsFactors=FALSE)

# treatment
if (all(mp14$treatment %in% c(0, 1, "0", "1"))) {
  mp14$treatment <- ifelse(as.numeric(mp14$treatment) == 1, "ACT", "OBS")
}

plot_km_compact <- function(data, mp_name, resp_data) {
  med <- median(data[[mp_name]], na.rm=TRUE)
  data$mp_grp <- ifelse(data[[mp_name]] >= med, "High", "Low")
  data$strata <- paste0(data$mp_grp, "-", data$treatment)
  data$stime <- data$os_time
  data$sevent <- data$os_event
  
  fit <- survfit(Surv(stime, sevent) ~ strata, data=data)
  sf <- summary(fit)
  
  km_df <- data.frame(
    time = sf$time, surv = sf$surv,
    strata = gsub("strata=", "", sf$strata)
  )
  starts <- data.frame(time=0, surv=1, strata=unique(km_df$strata))
  km_df <- rbind(starts, km_df)
  km_df <- km_df[order(km_df$strata, km_df$time), ]
  
  inter_row <- resp_data[resp_data$score == mp_name, ]
  if (nrow(inter_row) > 1) inter_row <- inter_row[1, ]
  inter_p <- if(nrow(inter_row)>0) inter_row$interaction_p else NA
  hr_h <- if(nrow(inter_row)>0) inter_row$hr_act_in_mp_high else NA
  hr_l <- if(nrow(inter_row)>0) inter_row$hr_act_in_mp_low else NA
  
  km_cols <- c("High-ACT"="#B2182B","High-OBS"="#EF8A62",
               "Low-ACT"="#2166AC","Low-OBS"="#67A9CF")
  km_lty <- c("High-ACT"="solid","High-OBS"="dashed",
              "Low-ACT"="solid","Low-OBS"="dashed")
  
  sig_mark <- if(!is.na(inter_p) && inter_p < 0.05) "*" else 
    if(!is.na(inter_p) && inter_p < 0.1) "†" else ""
  
  ggplot(km_df, aes(x=time, y=surv, color=strata, linetype=strata)) +
    geom_step(linewidth=0.7) +
    scale_color_manual(values=km_cols, name=NULL) +
    scale_linetype_manual(values=km_lty, name=NULL) +
    annotate("text", x=0.2, y=0.08,
             label=sprintf("inter %s%s\nHR(H)=%.2f\nHR(L)=%.2f",
                           fmt_p(inter_p), sig_mark, hr_h, hr_l),
             hjust=0, size=2.3, family=my_font, color="grey30") +
    coord_cartesian(ylim=c(0,1)) +
    labs(x="Time (yr)", y="OS", title=mp_name) +
    theme_elegant(base_size=8) +
    theme(legend.position="none",
          plot.title=element_text(size=9, color=mp_colors[mp_name]))
}

# GSE14814 4 MP
km14_list <- lapply(c("MP1","MP2","MP3","MP4"), function(m) {
  plot_km_compact(mp14, m, resp14)
})

p_km14 <- wrap_plots(km14_list, nrow=1) +
  plot_annotation(title="GSE14814 (JBR.10, n=133)   MP × ACT Stratified KM",
                  theme=theme(plot.title=element_text(face="bold", size=11, family=my_font)))

ggsave("supp_fig_act_km_gse14814.pdf", p_km14, width=12, height=3.5)
cat("   supp_fig_act_km_gse14814.pdf\n")

# GSE42127
mp42 <- read.csv("gse42127_mp_scores.csv", stringsAsFactors=FALSE)
resp42_file <- "gse42127_response_comparison.csv"

cat("  GSE42127 columns:", paste(colnames(mp42), collapse=", "), "\n")

has_surv_42 <- all(c("os_time","os_event") %in% colnames(mp42)) ||
  all(c("time","event") %in% colnames(mp42))

if (has_surv_42) {
  if (!"os_time" %in% colnames(mp42)) {
    tc <- intersect(c("time","survival_time"), colnames(mp42))[1]
    ec <- intersect(c("event","status"), colnames(mp42))[1]
    if (!is.na(tc)) mp42$os_time <- mp42[[tc]]
    if (!is.na(ec)) mp42$os_event <- mp42[[ec]]
  }
  
  if (all(mp42$treatment %in% c(0, 1, "0", "1"))) {
    mp42$treatment <- ifelse(as.numeric(mp42$treatment) == 1, "ACT", "OBS")
  }
  
  resp42 <- if(file.exists(resp42_file)) read.csv(resp42_file, stringsAsFactors=FALSE) else NULL
  
  if (!is.null(resp42)) {
    km42_list <- lapply(c("MP1","MP2","MP3","MP4"), function(m) {
      plot_km_compact(mp42, m, resp42)
    })
    
    p_km42 <- wrap_plots(km42_list, nrow=1) +
      plot_annotation(title="GSE42127 (UT SPORE, n=176)   MP × ACT Stratified KM",
                      theme=theme(plot.title=element_text(face="bold", size=11, family=my_font)))
    
    ggsave("supp_fig_act_km_gse42127.pdf", p_km42, width=12, height=3.5)
    cat("   supp_fig_act_km_gse42127.pdf\n")
    
    # ACT
    p_act_all <- p_km14 / p_km42
    ggsave("supp_fig_act_km_combined.pdf", p_act_all, width=12, height=7)
    cat("   supp_fig_act_km_combined.pdf\n")
  }
} else {
  cat("  GSE42127: no survival columns, skipping KM\n")
}

# Supp Fig C: AUC + Interaction p

cat("\n=== Supp Fig: Heatmap summary ===\n")

# IO: AUC
io_heat <- summ %>%
  filter(analysis_type == "IO", score %in% c("MP1","MP2","MP3","MP4")) %>%
  select(cohort, score, auc, wilcoxon_p) %>%
  mutate(label = sprintf("%.2f\n(%s)", auc, sapply(wilcoxon_p, fmt_p)))

io_heat$score <- factor(io_heat$score, levels=c("MP1","MP2","MP3","MP4"))
io_heat$cohort <- factor(io_heat$cohort)

p_heat_io <- ggplot(io_heat, aes(x=score, y=cohort, fill=auc)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=label), size=3, family=my_font, color="black") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0.5, limits=c(0,1), name="AUC") +
  labs(x=NULL, y=NULL, title="IO Response   AUC Summary") +
  theme_elegant(base_size=10) +
  theme(axis.text.x=element_text(face="bold"))

# ACT: interaction p
act_heat <- summ %>%
  filter(analysis_type == "ACT", score %in% c("MP1","MP2","MP3","MP4"),
         !is.na(interaction_p)) %>%
  select(cohort, score, interaction_p, interaction_hr) %>%
  mutate(
    neg_log_p = -log10(interaction_p),
    label = sprintf("HR=%.2f\n(%s)", interaction_hr, sapply(interaction_p, fmt_p))
  )

act_heat$score <- factor(act_heat$score, levels=c("MP1","MP2","MP3","MP4"))

p_heat_act <- ggplot(act_heat, aes(x=score, y=cohort, fill=neg_log_p)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=label), size=3, family=my_font) +
  scale_fill_gradient(low="white", high="#E64B35", name="-log10(p)") +
  labs(x=NULL, y=NULL, title="ACT Benefit   Interaction Test Summary") +
  theme_elegant(base_size=10) +
  theme(axis.text.x=element_text(face="bold"))

p_heat_combined <- p_heat_io / p_heat_act + 
  plot_layout(heights=c(1, 1)) +
  plot_annotation(tag_levels="A")

ggsave("supp_fig_heatmap_summary.pdf", p_heat_combined, width=7, height=6)
cat("   supp_fig_heatmap_summary.pdf\n")


cat("\n All supplementary outputs saved:\n")
cat("  Table: supp_table_treatment_full.csv\n")
cat("  Fig:   supp_fig_io_full.pdf (3 cohorts × 4 MPs boxplots)\n")
cat("  Fig:   supp_fig_act_km_gse14814.pdf (4 MPs stratified KM)\n")
cat("  Fig:   supp_fig_heatmap_summary.pdf (AUC + interaction p heatmap)\n")