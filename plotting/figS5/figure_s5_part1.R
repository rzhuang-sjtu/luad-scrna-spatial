## Supplementary Figure 5  |  GSEA Hallmark (4 panels)
# NES diverging bar chart MP panel

library(data.table)
library(ggplot2)
suppressPackageStartupMessages({
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig2" else
  "${WORK_ROOT}/luad_figures/fig2")

.find_arial <- function() {
  for (p in c("arial.ttf", "C:/Windows/Fonts/arial.ttf",
              "~/.local/share/fonts/arial.ttf",
              "/mnt/c/Windows/Fonts/arial.ttf")) {
    pp <- path.expand(p)
    if (file.exists(pp) || p == "arial.ttf") return(p)
  }
  NA_character_
}
my_font <- "sans"
if (requireNamespace("showtext", quietly = TRUE)) {
  .ar <- .find_arial()
  if (!is.na(.ar)) {
    .dir <- dirname(path.expand(.ar))
    .bd  <- file.path(.dir, "arialbd.ttf"); if (!file.exists(.bd)) .bd <- .ar
    .it  <- file.path(.dir, "ariali.ttf"); if (!file.exists(.it)) .it <- .ar
    tryCatch({
      sysfonts::font_add("Arial", regular = .ar, bold = .bd, italic = .it)
      showtext::showtext_auto(); showtext::showtext_opts(dpi = 300)
      my_font <- "Arial"
    }, error = function(e) {})
  }
}

mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")
mp_titles <- c(
  "MP1"="MP1: Stress/AP-1",
  "MP2"="MP2: Proliferative",
  "MP3"="MP3: EMT/IFN",
  "MP4"="MP4: AT2-like"
)

theme_pub <- function(base_size = 8) {
  theme_classic(base_family = my_font, base_size = base_size) +
    theme(axis.text       = element_text(color = "black"),
          axis.line       = element_line(linewidth = 0.4, color = "black"),
          axis.ticks      = element_line(linewidth = 0.3, color = "black"),
          axis.ticks.length = unit(1.5, "pt"),
          legend.title    = element_text(size = rel(0.95), face = "bold"),
          legend.text     = element_text(size = rel(0.85)),
          legend.key.size = unit(3, "mm"),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text      = element_text(face = "bold", size = rel(1)),
          plot.title      = element_text(face = "bold", size = rel(1.1)))
}

dir.create("../fig_s5", showWarnings=FALSE)

df <- fread("hallmark_gsea_top20_per_mp.csv")
cat("columns:", paste(names(df), collapse=", "), "\n")
cat("shape:", nrow(df), "x", ncol(df), "\n")
cat("MP counts:\n")
print(table(df$MP))

# Term +
df[, Term_clean := gsub("^HALLMARK_", "", Term)]
df[, Term_clean := gsub("_", " ", Term_clean)]
if (requireNamespace("stringr", quietly = TRUE)) {
  df[, Term_clean := stringr::str_wrap(Term_clean, width = 32)]
}

# direction NES
if(!"direction" %in% names(df)) {
  df[, direction := ifelse(NES > 0, "Up", "Down")]
}

# FDR
df[, sig_label := ifelse(FDR < 0.001, "***",
                         ifelse(FDR < 0.01, "**",
                                ifelse(FDR < 0.05, "*",
                                       ifelse(FDR < 0.25, "\u2020", ""))))]

panel_map <- c("MP1"="a","MP2"="b","MP3"="c","MP4"="d")

# panel : A=MP4(AT2), B=MP2(Prolif), C=MP3(EMT), D=MP1(Stress

for(mp in c("MP1","MP2","MP3","MP4")) {
  
  d <- df[MP == mp]
  d <- d[order(NES)]
  d[, Term_clean := factor(Term_clean, levels=Term_clean)]
  
  p <- ggplot(d, aes(x=NES, y=Term_clean, fill=direction)) +
    geom_col(width=0.7) +
    geom_vline(xintercept=0, linewidth=0.4, color="black") +
    # NES
    geom_text(aes(label=sprintf("%.2f", NES),
                  hjust=ifelse(NES > 0, -0.1, 1.1)),
              size=2.2, family=my_font) +
    # FDR
    geom_text(aes(label=sig_label, x=0,
                  hjust=ifelse(NES > 0, 1.3, -0.3)),
              size=2.5, family=my_font, color="grey40") +
    scale_fill_manual(values=c("Up"="#D73027","Down"="#4575B4"),
                      labels=c("Up"="Activated","Down"="Suppressed"),
                      name=NULL) +
    labs(x="NES", y=NULL, title=mp_titles[mp]) +
    theme_pub(base_size=7) +
    theme(axis.text.y=element_text(size=6, lineheight=0.85),
          plot.title=element_text(size=9, face="bold", color=mp_colors[mp]),
          legend.position="top",
          legend.text=element_text(size=7),
          legend.key.size=unit(3,"mm"),
          plot.margin=margin(4, 8, 4, 4, "pt"))

  letter <- panel_map[mp]
  ggsave(sprintf("../fig_s5/fig_s5%s_%s.pdf", letter, mp),
         p, width=180, height=120, units="mm")
  ggsave(sprintf("../fig_s5/fig_s5%s_%s.png", letter, mp),
         p, width=180, height=120, units="mm", dpi=300)
  cat("  ", mp, "saved\n")
}

# ── 4 panel (2×2) ──
library(patchwork)

plots <- list()
for(mp in c("MP1","MP2","MP3","MP4")) {
  d <- df[MP == mp][order(NES)]
  d[, Term_clean := factor(Term_clean, levels=Term_clean)]
  
  plots[[mp]] <- ggplot(d, aes(x=NES, y=Term_clean, fill=direction)) +
    geom_col(width=0.7) +
    geom_vline(xintercept=0, linewidth=0.3, color="black") +
    geom_text(aes(label=sprintf("%.1f", NES),
                  hjust=ifelse(NES > 0, -0.1, 1.1)),
              size=2, family=my_font) +
    scale_fill_manual(values=c("Up"="#D73027","Down"="#4575B4"),
                      labels=c("Up"="Activated","Down"="Suppressed"),
                      name=NULL) +
    labs(x="NES", y=NULL, title=mp_titles[mp]) +
    theme_pub(base_size=6) +
    theme(axis.text.y=element_text(size=5.5, lineheight=0.85),
          plot.title=element_text(size=8, face="bold", color=mp_colors[mp]),
          legend.position="none")
}

p_combined <- patchwork::wrap_plots(plots, ncol=2) +
  patchwork::plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(size=11, face="bold", family=my_font))

ggsave("../fig_s5/FigureS5_NES.pdf", p_combined,
       width=380, height=240, units="mm")
ggsave("../fig_s5/FigureS5_NES.png", p_combined,
       width=380, height=240, units="mm", dpi=300)
# Backwards compat: keep old name too
ggsave("../fig_s5/fig_s5_combined.pdf", p_combined,
       width=380, height=240, units="mm")
ggsave("../fig_s5/fig_s5_combined.png", p_combined,
       width=380, height=240, units="mm", dpi=300)
cat("  combined saved\n")

message("=== Supp Fig 5 done ===")
