## Figure S10 — Supplementary spatial validation
## Layout (10 panels):
##   S10A : E-MTAB section #3 (different patient) cell2location 2×3 grid
##   S10B : E-MTAB section #4 (different patient) cell2location 2×3 grid
##   S10C : Okamura LUAD_No_1 cell2location 2×3 grid (validation cohort)
##   S10D : Okamura LUAD_No_1 composite cell-type map
##   S10E : Okamura LUAD_No_1 COMMOT (OSM vectors + OSM + IL1B expr)
##   S10F : Okamura LUAD_No_5 COMMOT
##   S10G : Okamura LUAD_No_17 COMMOT (third strongest)
##   S10H : PROGENy pathway activity ranking, dual cohort, mean ± SEM
##   S10I : Spot-level Neutrophil-total vs MP3 (EMT/IFN) scatter, dual cohort
##   S10J : Dual-cohort ROI vs non-ROI delta barplot (top 20 metrics)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(scales); library(dplyr); library(tidyr)
  library(viridis); library(grid); library(png); library(jsonlite)
  library(ggnewscale); library(showtext); library(sysfonts); library(ggrepel)
})

R_DATA   <- "${DATA_ROOT}/ST/results/r_data"
PER      <- file.path(R_DATA, "per_section")
HE_DIR   <- file.path(R_DATA, "he")
FIG_DIR  <- "${WORK_ROOT}/luad_figures/fig_s10"
PANEL_DIR <- file.path(FIG_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

## fonts
arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p,
    bold = "~/.local/share/fonts/arialbd.ttf",
    italic = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300); my_font <- "Arial"
} else { my_font <- "sans" }

SEC_S10A <- list(stem = "EMTAB13530__P15_T1",  label = "P15_T1 · E-MTAB-13530 · LUAD")
SEC_S10B <- list(stem = "EMTAB13530__P24_T1",  label = "P24_T1 · E-MTAB-13530 · LUAD")
SEC_S10C <- list(stem = "Okamura__LUAD_No_1",  label = "LUAD_No_1 · Takano 2024 · LUAD")
SEC_S10E <- list(stem = "Okamura__LUAD_No_1",  label = "LUAD_No_1 · Takano 2024 · LUAD")
SEC_S10F <- list(stem = "Okamura__LUAD_No_5",  label = "LUAD_No_5 · Takano 2024 · LUAD")
SEC_S10G <- list(stem = "Okamura__LUAD_No_17", label = "LUAD_No_17 · Takano 2024 · LUAD")

neu_colors <- c(
  Neu_Inflammatory = "#E64B35", Neu_Angiogenic = "#F39B7F",
  Neu_Metastatic = "#3C5488", Neu_ECM_remodeling = "#4DBBD5",
  Neu_OSM_priming = "#00A087", Neu_OSM_low = "#8491B4",
  Neu_IFN_response = "#91D1C2")
macro_colors <- c(Macro_C1QC = "#4DBBD5", Macro_FCN1 = "#E64B35",
  Macro_FOLR2 = "#00A087", Macro_MARCO = "#3C5488",
  Macro_SPP1 = "#F39B7F", Macro_general = "#8491B4", Macro_prolif = "#91D1C2")
ct_colors <- c(Fibroblast = "#B09C85", Endothelial = "#651FFF",
  T_NK = "#FFCB5C", B = "#9C27B0", Plasma = "#D9D9D9", Mast = "#F0A23B",
  Epithelial = "#D2B48C", Mono_nonclassical = "#999999",
  cDC1 = "#80CBC4", cDC2 = "#7570B3", cDC_LAMP3 = "#984EA3", pDC = "#B2B2B2",
  Malignant = "#C73E2A")

cohort_colors <- c(EMTAB13530 = "#0072B2", Okamura = "#D55E00")

theme_panel <- function(base_size = 8) {
  theme_void(base_family = my_font, base_size = base_size) +
    theme(plot.subtitle = element_text(face = "bold", size = rel(1),
                                        margin = margin(b = 1.5)),
          plot.caption  = element_text(size = rel(0.75), color = "grey25",
                                        hjust = 0.5, margin = margin(t = 1)),
          legend.title = element_text(size = rel(0.8)),
          legend.text  = element_text(size = rel(0.75)),
          legend.key.width  = unit(0.18, "cm"),
          legend.key.height = unit(0.42, "cm"),
          plot.margin = margin(2,2,2,2))
}

theme_xy <- function(base_size = 8) {
  theme_classic(base_family = my_font, base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3),
          plot.subtitle = element_text(face = "bold", size = rel(1)),
          plot.title = element_blank(),
          legend.text = element_text(size = rel(0.85)),
          legend.title = element_text(size = rel(0.85)),
          plot.margin = margin(3,3,3,3))
}

load_section <- function(stem) fread(file.path(PER, paste0(stem, ".csv")))
load_he <- function(stem) {
  d <- file.path(HE_DIR, stem)
  img <- png::readPNG(file.path(d, "tissue_hires_image.png"))
  sf  <- jsonlite::fromJSON(file.path(d, "scalefactors_json.json"))
  list(img = img, sf = sf$tissue_hires_scalef,
       w = dim(img)[2], h = dim(img)[1])
}
prep_xy <- function(d, he) { d$x_he <- d$spatial1*he$sf; d$y_he <- -d$spatial2*he$sf; d }

plot_spatial_he <- function(d, color_col, he, panel_name,
                             high_color = "#B2182B", low_color = "white",
                             value_clip = c(0.02, 0.98), point_size = 0.85,
                             alpha_max = 0.92, diverging = FALSE) {
  v <- d[[color_col]]
  if (diverging) {
    a <- max(abs(quantile(v, value_clip[1], na.rm=TRUE)),
             abs(quantile(v, value_clip[2], na.rm=TRUE)))
    lim <- c(-a, a); a_norm <- pmin(1, abs(v) / lim[2])
  } else {
    lim <- c(quantile(v, value_clip[1], na.rm=TRUE),
             quantile(v, value_clip[2], na.rm=TRUE))
    rng <- diff(lim); a_norm <- pmin(1, pmax(0, (v - lim[1]) / ifelse(rng>0, rng, 1)))
  }
  d$alpha_pt <- alpha_max * a_norm
  p <- ggplot() +
    annotation_raster(he$img, xmin=0, xmax=he$w, ymin=-he$h, ymax=0) +
    geom_point(data=d, aes(x=x_he, y=y_he, color=.data[[color_col]], alpha=alpha_pt),
               size=point_size, stroke=0)
  if (diverging) {
    p <- p + scale_color_gradient2(low="#2166AC", mid="white", high="#B2182B",
                                     midpoint=0, limits=lim, oob=scales::squish, name="")
  } else {
    p <- p + scale_color_gradient(low=low_color, high=high_color,
                                    limits=lim, oob=scales::squish, name="")
  }
  p + scale_alpha_identity() +
    coord_fixed(xlim=c(0,he$w), ylim=c(-he$h,0), expand=FALSE) +
    theme_panel() + labs(subtitle = panel_name)
}

cell_panels_main <- list(
  list(col="ct_Neu_Inflammatory", name="Neu_Inflammatory", color=unname(neu_colors["Neu_Inflammatory"])),
  list(col="ct_Neu_OSM_priming",  name="Neu_OSM_priming",  color=unname(neu_colors["Neu_OSM_priming"])),
  list(col="ct_Macro_SPP1",       name="Macro_SPP1",       color=unname(macro_colors["Macro_SPP1"])),
  list(col="ct_Fibroblast",       name="Fibroblast",       color="#5C4830"),
  list(col="ct_Endothelial",      name="Endothelial",      color=unname(ct_colors["Endothelial"])),
  list(col="ct_Malignant",        name="Malignant",        color=unname(ct_colors["Malignant"]))
)

make_grid <- function(d, he, sample_label, panels_list) {
  ps <- lapply(panels_list, function(spec) {
    plot_spatial_he(d, spec$col, he, spec$name,
                    high_color = spec$color, point_size = 0.85, alpha_max = 0.95)
  })
  wrap_plots(ps, nrow=2, ncol=3) +
    plot_annotation(title = sample_label,
      theme = theme(plot.title = element_text(face="bold", size=9,
                                                family=my_font, hjust=0.02)))
}

DA1 <- prep_xy(load_section(SEC_S10A$stem), HEA1 <- load_he(SEC_S10A$stem))
DB1 <- prep_xy(load_section(SEC_S10B$stem), HEB1 <- load_he(SEC_S10B$stem))
cat(sprintf("[S10A] %s spots=%d\n", SEC_S10A$stem, nrow(DA1)))
cat(sprintf("[S10B] %s spots=%d\n", SEC_S10B$stem, nrow(DB1)))

P_S10A <- make_grid(DA1, HEA1, SEC_S10A$label, cell_panels_main)
ggsave(file.path(PANEL_DIR, "S10A.pdf"), P_S10A, width=130, height=90, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10A.png"), P_S10A, width=130, height=90, units="mm", dpi=300)

P_S10B <- make_grid(DB1, HEB1, SEC_S10B$label, cell_panels_main)
ggsave(file.path(PANEL_DIR, "S10B.pdf"), P_S10B, width=130, height=90, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10B.png"), P_S10B, width=130, height=90, units="mm", dpi=300)

DC1 <- prep_xy(load_section(SEC_S10C$stem), HEC1 <- load_he(SEC_S10C$stem))
cat(sprintf("[S10C] %s spots=%d\n", SEC_S10C$stem, nrow(DC1)))
P_S10C <- make_grid(DC1, HEC1, SEC_S10C$label, cell_panels_main)
ggsave(file.path(PANEL_DIR, "S10C.pdf"), P_S10C, width=130, height=90, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10C.png"), P_S10C, width=130, height=90, units="mm", dpi=300)

plot_composite_he <- function(d, he, sample_label, focus_cts, palette,
                               top_q = 0.5, point_size = 1.1, alpha_pt = 0.9) {
  ab_cols <- paste0("ct_", focus_cts)
  ab_cols <- ab_cols[ab_cols %in% names(d)]
  m <- as.matrix(d[, ..ab_cols]); colnames(m) <- sub("^ct_","", ab_cols)
  best_idx <- apply(m, 1, which.max)
  d2 <- copy(d)
  d2$dom_ct  <- factor(colnames(m)[best_idx], levels=colnames(m))
  d2$dom_val <- m[cbind(seq_len(nrow(m)), best_idx)]
  high_mask <- rep(FALSE, nrow(d2))
  for (ct in colnames(m)) {
    idx <- d2$dom_ct == ct; if (sum(idx)<5) next
    th <- quantile(d2$dom_val[idx], top_q, na.rm=TRUE)
    high_mask <- high_mask | (idx & d2$dom_val >= th)
  }
  d2$is_high <- high_mask
  ggplot() +
    annotation_raster(he$img, xmin=0, xmax=he$w, ymin=-he$h, ymax=0) +
    geom_point(data=d2[is_high == TRUE], aes(x=x_he, y=y_he, color=dom_ct),
               size=point_size, stroke=0, alpha=alpha_pt) +
    scale_color_manual(values=palette, name=NULL,
                       guide=guide_legend(override.aes=list(size=2.4, alpha=1))) +
    coord_fixed(xlim=c(0,he$w), ylim=c(-he$h,0), expand=FALSE) +
    theme_panel() +
    theme(legend.position="right", legend.key.size=unit(0.32,"cm"),
          legend.text=element_text(size=7)) +
    labs(subtitle="Composite (top-50% per cell type)", caption=sample_label)
}

focus_cts <- c("Neu_Inflammatory","Neu_OSM_priming","Macro_SPP1",
               "Fibroblast","Endothelial","Malignant")
focus_palette <- c(
  Neu_Inflammatory = unname(neu_colors["Neu_Inflammatory"]),
  Neu_OSM_priming  = unname(neu_colors["Neu_OSM_priming"]),
  Macro_SPP1       = unname(macro_colors["Macro_SPP1"]),
  Fibroblast       = "#5C4830",
  Endothelial      = unname(ct_colors["Endothelial"]),
  Malignant        = unname(ct_colors["Malignant"]))
P_S10D <- plot_composite_he(DC1, HEC1, SEC_S10C$label, focus_cts, focus_palette)
ggsave(file.path(PANEL_DIR, "S10D.pdf"), P_S10D, width=100, height=95, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10D.png"), P_S10D, width=100, height=95, units="mm", dpi=300)

plot_commot_field_he <- function(d, he, sample_label, pathway = "OSM",
                                  arrow_topN = 60, target_arrow_frac = 0.10) {
  send_col <- sprintf("commot_s_%s", pathway)
  tot_col  <- sprintf("commot_total_%s", pathway)
  vfx <- sprintf("vf_s_%s_dx", pathway); vfy <- sprintf("vf_s_%s_dy", pathway)
  d <- d[!is.na(get(vfx)) & !is.na(get(vfy))]
  d$tot <- d[[tot_col]]; d$tot[is.na(d$tot)] <- 0
  span <- max(diff(range(d$x_he)), diff(range(d$y_he)))
  vmag <- sqrt(d[[vfx]]^2 + d[[vfy]]^2)
  med <- median(vmag[vmag > 0], na.rm = TRUE)
  scale_factor <- if (is.na(med) || med == 0) 0 else span * target_arrow_frac / med
  ord <- order(d[[send_col]], decreasing = TRUE)
  arr <- d[ord[seq_len(min(arrow_topN, nrow(d)))]]
  arr$ax <- arr$x_he; arr$ay <- arr$y_he
  arr$bx <- arr$x_he + scale_factor * arr[[vfx]]
  arr$by <- arr$y_he - scale_factor * arr[[vfy]]
  d$alpha_pt <- pmin(1, d$tot / quantile(d$tot[d$tot>0], 0.95, na.rm=TRUE))*0.95
  d$alpha_pt[is.na(d$alpha_pt)] <- 0
  ggplot() +
    annotation_raster(he$img, xmin=0, xmax=he$w, ymin=-he$h, ymax=0) +
    geom_point(data=d, aes(x=x_he, y=y_he, color=tot, alpha=alpha_pt),
               size=0.8, stroke=0) +
    scale_color_gradient(low="#FFFAF0", high="#B2182B",
                          limits=c(0, max(quantile(d$tot[d$tot>0], 0.95, na.rm=TRUE), 1e-6)),
                          oob=scales::squish, name=sprintf("%s total", pathway)) +
    scale_alpha_identity() +
    geom_segment(data=arr, aes(x=ax, y=ay, xend=bx, yend=by),
                 arrow=arrow(length=unit(0.10,"cm"), type="closed"),
                 color="#003C30", linewidth=0.32) +
    coord_fixed(xlim=c(0,he$w), ylim=c(-he$h,0), expand=FALSE) +
    theme_panel() +
    labs(subtitle=sprintf("COMMOT %s sender vectors", pathway), caption=sample_label)
}

make_commot_set <- function(stem, sample_label, h_or_v = "h") {
  d <- prep_xy(load_section(stem), he <- load_he(stem))
  vec <- plot_commot_field_he(d, he, sample_label, "OSM")
  osm <- plot_spatial_he(d, "gex_OSM",  he, "OSM expression",
                          high_color = "#7B0033", point_size = 1.0)
  il1 <- plot_spatial_he(d, "gex_IL1B", he, "IL1B expression",
                          high_color = "#005C5A", point_size = 1.0)
  osm <- osm + labs(caption = sample_label)
  il1 <- il1 + labs(caption = sample_label)
  vec + osm + il1 + plot_layout(nrow = 1)
}

P_S10E <- make_commot_set(SEC_S10E$stem, SEC_S10E$label)
ggsave(file.path(PANEL_DIR, "S10E.pdf"), P_S10E, width=200, height=80, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10E.png"), P_S10E, width=200, height=80, units="mm", dpi=300)
P_S10F <- make_commot_set(SEC_S10F$stem, SEC_S10F$label)
ggsave(file.path(PANEL_DIR, "S10F.pdf"), P_S10F, width=200, height=80, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10F.png"), P_S10F, width=200, height=80, units="mm", dpi=300)
P_S10G <- make_commot_set(SEC_S10G$stem, SEC_S10G$label)
ggsave(file.path(PANEL_DIR, "S10G.pdf"), P_S10G, width=200, height=80, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10G.png"), P_S10G, width=200, height=80, units="mm", dpi=300)

## Aggregate per-section mean PROGENy score, then plot by cohort with error bars.
all_csvs <- list.files(PER, pattern = "\\.csv$", full.names = TRUE)
get_progeny_means <- function(fp) {
  d <- fread(fp)
  pcols <- grep("^progeny_", names(d), value = TRUE)
  out <- as.data.table(t(colMeans(d[, ..pcols], na.rm=TRUE)))
  out$sample <- d$sample[1]; out$cohort <- d$cohort[1]
  out
}
prog_long <- rbindlist(lapply(all_csvs, get_progeny_means), fill = TRUE)
prog_long <- melt(prog_long, id.vars = c("sample","cohort"),
                  variable.name = "pathway", value.name = "score")
prog_long[, pathway := sub("^progeny_", "", pathway)]
prog_long[, pathway := gsub("\\.", "-", pathway)]
prog_summary <- prog_long[, .(mean = mean(score, na.rm=TRUE),
                               sem  = sd(score, na.rm=TRUE)/sqrt(.N),
                               n    = .N),
                          by = .(pathway, cohort)]
## order by overall mean across cohorts
order_pw <- prog_summary[, .(overall = mean(mean)), by = pathway][order(-overall), pathway]
prog_summary[, pathway := factor(pathway, levels = order_pw)]

P_S10H <- ggplot(prog_summary,
                 aes(x = pathway, y = mean, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           color = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.78),
                width = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = cohort_colors,
                    labels = c(EMTAB13530 = "E-MTAB-13530",
                               Okamura    = "Takano 2024"),
                    name = "Cohort") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  labs(subtitle = "PROGENy pathway activity (mean ± SEM across sections)",
       x = "PROGENy pathway", y = "MLM activity score") +
  theme_xy() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7))
ggsave(file.path(PANEL_DIR, "S10H.pdf"), P_S10H, width=150, height=70, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10H.png"), P_S10H, width=150, height=70, units="mm", dpi=300)

get_scatter_data <- function(fp) {
  d <- fread(fp)
  d[, .(neu_total, MP3_score, sample, cohort)]
}
scatter_dt <- rbindlist(lapply(all_csvs, get_scatter_data), fill = TRUE)
scatter_dt <- scatter_dt[is.finite(neu_total) & is.finite(MP3_score)]

## per-cohort Spearman on FULL data (not subsample)
sp_em <- cor.test(scatter_dt[cohort=="EMTAB13530", neu_total],
                  scatter_dt[cohort=="EMTAB13530", MP3_score], method = "spearman")
sp_ok <- cor.test(scatter_dt[cohort=="Okamura",    neu_total],
                  scatter_dt[cohort=="Okamura",    MP3_score], method = "spearman")

## p-value formatter: threshold-text (avoid scientific notation in figures)
fmt_p <- function(p) {
  if (!is.finite(p) || p < .Machine$double.xmin) return("p < 0.0001 ****")
  if (p < 1e-4) return("p < 0.0001 ****")
  if (p < 1e-3) return("p < 0.001 ***")
  if (p < 1e-2) return("p < 0.01 **")
  if (p < 0.05) return(sprintf("p = %.3f *", p))
  return(sprintf("p = %.3f ns", p))
}

## subsample for rendering only (correlation already computed on full data)
set.seed(0)
plot_dt <- scatter_dt[sample(.N, min(20000, .N))]

## clip y to robust range to suppress outliers
y_lo <- quantile(plot_dt$MP3_score, 0.005, na.rm=TRUE)
y_hi <- quantile(plot_dt$MP3_score, 0.995, na.rm=TRUE)
plot_dt <- plot_dt[MP3_score >= y_lo & MP3_score <= y_hi]

## per-cohort linear-fit lines drawn in transformed space
ann_text <- data.frame(
  cohort = c("EMTAB13530","Okamura"),
  label  = c(sprintf("EMTAB13530 (n=%s)\nρ = %.3f, %s",
                     formatC(scatter_dt[cohort=="EMTAB13530", .N], format="d", big.mark=","),
                     sp_em$estimate, fmt_p(sp_em$p.value)),
             sprintf("Takano (n=%s)\nρ = %.3f, %s",
                     formatC(scatter_dt[cohort=="Okamura", .N], format="d", big.mark=","),
                     sp_ok$estimate, fmt_p(sp_ok$p.value))))

## use log1p x-axis to spread out 0-5 pile-up; faceted to avoid overlap
P_S10I <- ggplot(plot_dt, aes(x = neu_total, y = MP3_score)) +
  geom_point(aes(color = cohort), size = 0.35, stroke = 0, alpha = 0.28) +
  geom_smooth(aes(color = cohort), method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = cohort_colors,
                     labels = c(EMTAB13530 = "E-MTAB-13530",
                                Okamura    = "Takano 2024"),
                     name = "Cohort", guide = "none") +
  scale_x_continuous(transform = "log1p",
                     breaks = c(0, 1, 3, 10, 30, 100),
                     labels = c("0", "1", "3", "10", "30", "100")) +
  geom_label(data = ann_text, aes(label = label),
             x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             size = 2.3, family = my_font, label.size = 0.2,
             fill = alpha("white", 0.92), label.padding = unit(1.2, "mm")) +
  facet_wrap(~ cohort, nrow = 1, scales = "free_y") +
  labs(subtitle = "Spot-level Neutrophil-total vs MP3 (EMT/IFN) score (log1p x-axis)",
       x = "Neutrophil total abundance (sum of 7 subtypes, log1p)",
       y = "MP3 (EMT/IFN) score") +
  theme_xy() +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(size = 8, face = "bold"))

ggsave(file.path(PANEL_DIR, "S10I.pdf"), P_S10I, width=160, height=85, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10I.png"), P_S10I, width=160, height=85, units="mm", dpi=300)

## Now reads the file with Mann-Whitney p-values + BH-FDR (per-spot test
## pooled within cohort across sections) and overlays significance stars.
roi_pv_path <- file.path(R_DATA, "roi_vs_nonroi_aggregate_pvalues.csv")
if (file.exists(roi_pv_path)) {
  roi_agg <- fread(roi_pv_path)
} else {
  roi_agg <- fread(file.path(R_DATA, "roi_vs_nonroi_aggregate.csv"))
  roi_agg[, `:=`(p_fdr = NA_real_, sig = "")]
}
roi_agg[, type := factor(type, levels = c("celltype","obs","gene"),
                          labels = c("Cell type","Pathway / MP","Gene"))]
metric_in_both <- roi_agg[, .(n_coh = uniqueN(cohort)), by = metric][n_coh == 2, metric]
roi_pair <- roi_agg[metric %in% metric_in_both]
mean_delta <- roi_pair[, .(absd = mean(abs(delta), na.rm=TRUE)), by = metric]
top_metrics <- mean_delta[order(-absd)][1:24, metric]
roi_pair <- roi_pair[metric %in% top_metrics]
roi_pair[, metric := factor(metric, levels = rev(top_metrics))]
## position significance label just past the bar tip (sign-aware)
roi_pair[, sig_x := ifelse(delta >= 0,
                            delta + 0.04 * max(abs(delta), na.rm = TRUE),
                            delta - 0.04 * max(abs(delta), na.rm = TRUE))]
roi_pair[is.na(sig) | sig == "", sig := ""]

## Bubble heatmap (cohort × metric) with sig stars overlaid in white.
roi_pair[, cohort := factor(cohort, levels = c("EMTAB13530","Okamura"))]
lim_d <- max(abs(roi_pair$delta), na.rm = TRUE)

P_S10J <- ggplot(roi_pair, aes(x = cohort, y = metric)) +
  geom_point(aes(size = abs(delta), fill = delta),
             shape = 21, color = "grey25", stroke = 0.25) +
  geom_text(aes(label = sig), vjust = 0.55, hjust = 0.5,
            size = 1.8, family = FAM, fontface = "bold", color = "white") +
  scale_x_discrete(labels = c(EMTAB13530 = "E-MTAB-13530",
                              Okamura    = "Takano 2024"),
                   position = "top",
                   expand = expansion(add = c(0.7, 0.7))) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, limits = c(-lim_d, lim_d),
                       oob = scales::squish, name = "Delta",
                       guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                              barheight = unit(20, "mm"),
                                              frame.colour = "black",
                                              frame.linewidth = 0.25)) +
  scale_size_continuous(range = c(1.6, 7), name = "|Delta|",
                        guide = guide_legend(override.aes = list(fill = "grey55"))) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(subtitle = "ROI vs non-ROI · top-24 |delta| · per-spot Mann-Whitney + BH-FDR\n*** <0.001  ** <0.01  * <0.05",
       x = NULL, y = NULL) +
  theme_minimal(base_family = FAM, base_size = 7) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
        panel.background   = element_rect(fill = "white", color = NA),
        axis.text          = element_text(color = "black"),
        axis.text.x        = element_text(size = 7, face = "bold", margin = margin(t = 1)),
        axis.ticks         = element_blank(),
        axis.line          = element_blank(),
        strip.placement    = "outside",
        strip.background   = element_rect(fill = "grey94", color = NA),
        strip.text.y.left  = element_text(angle = 0, face = "bold",
                                          size = 6.5, color = "black"),
        plot.subtitle = element_text(size = 6, lineheight = 1.1, margin = margin(b = 2)),
        legend.position  = "right",
        legend.box       = "vertical",
        legend.key.size  = unit(3, "mm"),
        legend.text      = element_text(size = 6),
        legend.title     = element_text(size = 6, face = "bold"),
        plot.margin = margin(2, 4, 2, 4))
ggsave(file.path(PANEL_DIR, "S10J.pdf"), P_S10J, width=130, height=215, units="mm", device=cairo_pdf)
ggsave(file.path(PANEL_DIR, "S10J.png"), P_S10J, width=130, height=215, units="mm", dpi=300)

cat("[done] S10 panels saved -> ", PANEL_DIR, "\n")
