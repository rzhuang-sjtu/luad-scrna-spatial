# Supplementary Figure 5 v3 | NES

library(data.table)
library(ggplot2)
library(ggtext)    # install.packages("ggtext")
library(R.utils)
library(patchwork)
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

mp_titles <- c("MP1"="MP1","MP2"="MP2","MP3"="MP3","MP4"="MP4")

dir.create("../fig_s5", showWarnings=FALSE)

es   <- fread("gsea_running_es.csv.gz")
hits <- fread("gsea_hit_positions.csv")
top20 <- fread("hallmark_gsea_top20_per_mp.csv")

clean_term <- function(x) gsub("^HALLMARK_", "", x)

plot_mp <- function(mp_name, n_terms=13) {
  
  mp_top <- top20[MP == mp_name][order(-abs(NES))][1:min(n_terms, .N)]
  mp_top <- mp_top[order(-NES)]
  terms <- mp_top$Term
  
  es_sub <- es[MP == mp_name & Term %in% terms]
  hit_sub <- hits[MP == mp_name & Term %in% terms]
  
  nes_map <- setNames(mp_top$NES, mp_top$Term)
  
  # label NES HTML NES
  # Wrap long term names with HTML <br> so names stay legible in narrow strips.
  wrap_html <- function(s, width = 28) {
    if (nchar(s) <= width) return(s)
    words <- strsplit(s, " ")[[1]]
    out <- ""; line <- ""
    for (w in words) {
      if (nchar(line) + nchar(w) + 1 > width) {
        out <- paste0(out, line, "<br>"); line <- w
      } else {
        line <- if (line == "") w else paste(line, w)
      }
    }
    paste0(out, line)
  }

  label_map <- sapply(terms, function(t) {
    nes <- nes_map[t]
    name <- wrap_html(gsub("_", " ", clean_term(t)))
    if(nes > 0) {
      sprintf("<b style='color:#D73027'>%s</b>&nbsp;&nbsp;NES: %.2f", name, nes)
    } else {
      sprintf("<b>%s</b>&nbsp;&nbsp;NES: %.2f", name, nes)
    }
  })
  names(label_map) <- terms
  
  es_sub[, Term := factor(Term, levels=terms)]
  hit_sub[, Term := factor(Term, levels=terms)]
  
  p <- ggplot(es_sub, aes(x=rank, y=running_es)) +
    geom_hline(yintercept=0, linewidth=0.15, color="grey60", linetype="dashed") +
    geom_line(linewidth=0.3, color="black") +
    geom_segment(data=hit_sub,
                 aes(x=hit_rank, xend=hit_rank, y=-0.05, yend=0),
                 inherit.aes=FALSE, linewidth=0.08, color="grey30") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(breaks=NULL) +
    facet_grid(Term ~ ., scales="free_y", switch="y",
               labeller=labeller(Term=label_map)) +
    labs(x=NULL, y=NULL) +
    theme_void(base_family=my_font) +
    theme(
      # element_markdown HTML
      strip.text.y.left=element_markdown(angle=0, hjust=1, size=5,
                                         margin=margin(0,2,0,0,"mm")),
      strip.placement="outside",
      panel.spacing=unit(0.3, "mm"),
      plot.margin=margin(5,8,5,2,"mm"),
      plot.title=element_text(size=12, face="bold", family=my_font, hjust=0.5)
    ) +
    ggtitle(mp_titles[mp_name])
  
  p
}

for(mp in c("MP1","MP2","MP3","MP4")) {
  cat("── ", mp, " ──\n")
  p <- plot_mp(mp, n_terms=13)
  letter <- c("MP1"="a","MP2"="b","MP3"="c","MP4"="d")[mp]
  
  # Different filename pattern from S5-1 so they don't overwrite each other
  ggsave(sprintf("../fig_s5/fig_s5_es_%s_%s.pdf", letter, mp),
         p, width=160, height=170, units="mm")
  ggsave(sprintf("../fig_s5/fig_s5_es_%s_%s.png", letter, mp),
         p, width=160, height=170, units="mm", dpi=300)
  cat("  saved\n")
}

p_all <- patchwork::wrap_plots(
  list(plot_mp("MP1",13), plot_mp("MP2",13),
       plot_mp("MP3",13), plot_mp("MP4",13)),
  ncol = 2
) + patchwork::plot_annotation(tag_levels="a") &
    theme(plot.tag = element_text(size=11, face="bold", family=my_font))

ggsave("../fig_s5/FigureS5_ES.pdf", p_all,
       width=360, height=320, units="mm")
ggsave("../fig_s5/FigureS5_ES.png", p_all,
       width=360, height=320, units="mm", dpi=300)
cat("  combined saved (FigureS5_ES)\n")

message("=== Supp Fig 5 v3 done ===")
