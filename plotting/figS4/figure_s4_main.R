library(data.table)
setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

d <- fread("mp_marker_gene_order.csv")
cat("shape:", nrow(d), "x", ncol(d), "\n")
cat("columns:", paste(names(d), collapse=", "), "\n")
print(table(d$MP))
print(head(d, 8))

# top100
if(file.exists("mp_signatures_top100.csv")) {
  d2 <- fread("mp_signatures_top100.csv")
  cat("\ntop100 shape:", nrow(d2), "x", ncol(d2), "\n")
  print(table(d2$MP))
}

print(list.files(".", pattern="enrich|metascape|GO|gsea", ignore.case=TRUE))
## Supplementary Figure 4  |  Functional enrichment (4 panels)
# enrichR Metascape

# install.packages("enrichR") #
library(data.table)
library(ggplot2)
# enrichR CSV API
.use_cache <- file.exists("../fig_s4/enrichment_results_all_MPs.csv")
if (!.use_cache && requireNamespace("enrichR", quietly = TRUE)) library(enrichR)
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

# top100
sig <- fread("mp_signatures_top100.csv")
cat("columns:", paste(names(sig), collapse=", "), "\n")

# enrichR Metascape
dbs <- c("GO_Biological_Process_2023",
         "GO_Molecular_Function_2023",
         "KEGG_2021_Human",
         "Reactome_2022",
         "WikiPathway_2023_Human",
         "MSigDB_Hallmark_2020")

if (!.use_cache) setEnrichrSite("Enrichr")

dir.create("../fig_s4", showWarnings=FALSE)

all_results <- list()
panels <- list()

.cached <- if (.use_cache) {
  tmp <- as.data.table(fread("../fig_s4/enrichment_results_all_MPs.csv"))
  # CSV may contain duplicate columns from earlier runs (e.g. two `MP` cols)
  if (any(duplicated(names(tmp)))) {
    tmp <- tmp[, !duplicated(names(tmp)), with = FALSE]
  }
  tmp
} else NULL

build_s4_panel <- function(mp) {
  if (.use_cache) {
    res <- copy(.cached[MP == mp])
    if (nrow(res) == 0) return(NULL)
    # cache already has neg_log10_p + Term_clean drop them so we recompute
    drop_cols <- intersect(c("neg_log10_p", "Term_clean", "p_label", "db_short"),
                           names(res))
    if (length(drop_cols) > 0) res[, (drop_cols) := NULL]
  } else {
    genes <- sig[MP == mp]$gene
    enriched <- enrichr(genes, dbs)
    res_list <- lapply(names(enriched), function(db) {
      d <- enriched[[db]]
      if (nrow(d) == 0) return(NULL)
      d$Database <- db
      d
    })
    res <- do.call(rbind, res_list[!sapply(res_list, is.null)])
    if (nrow(res) == 0) return(NULL)
    res <- as.data.table(res)
  }
  res[, neg_log10_p := -log10(P.value)]
  res <- res[order(Adjusted.P.value)][1:min(15, .N)]
  res[, Term_clean := gsub(" \\(GO:\\d+\\)$", "", Term)]
  res[, Term_clean := gsub(" \\(R-HSA-\\d+\\)$", "", Term_clean)]
  res[, Term_clean := gsub(" WP\\d+$", "", Term_clean)]
  # Wrap term names at ~38 chars onto 2 lines instead of truncating
  # full content stays visible even for long pathway names.
  if (requireNamespace("stringr", quietly = TRUE)) {
    res[, Term_clean := stringr::str_wrap(Term_clean, width = 38)]
  } else {
    res[, Term_clean := ifelse(nchar(Term_clean) > 80,
                               paste0(substr(Term_clean, 1, 77), "..."),
                               Term_clean)]
  }
  res <- res[order(neg_log10_p)]
  res[, Term_clean := factor(Term_clean, levels = Term_clean)]
  all_results[[mp]] <<- res

  # Original gradient style (light → MP colour) user prefers
  # Wider canvas (170 mm) gives 60-char term names full room.
  ggplot(res, aes(x = neg_log10_p, y = Term_clean, fill = neg_log10_p)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    scale_fill_gradient(low  = adjustcolor(mp_colors[mp], alpha.f = 0.4),
                        high = mp_colors[mp]) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = expression(-log[10]~italic(P)), y = NULL,
         title = mp_titles[mp]) +
    theme_pub(base_size = 7) +
    theme(axis.text.y = element_text(size = 6, color = "black",
                                     lineheight = 0.85),
          axis.text.x = element_text(size = 6.5),
          plot.title  = element_text(size = 9, face = "bold",
                                     color = mp_colors[mp]),
          plot.margin = margin(4, 8, 4, 4, "pt"))
}

for (mp in c("MP1", "MP2", "MP3", "MP4")) {
  cat("── ", mp, " ──\n")
  p <- build_s4_panel(mp)
  if (is.null(p)) next
  panels[[mp]] <- p
  letter <- c("MP1"="a","MP2"="b","MP3"="c","MP4"="d")[mp]
  ggsave(sprintf("../fig_s4/fig_s4%s_%s.pdf", letter, mp), p,
         width = 195, height = 110, units = "mm")
  ggsave(sprintf("../fig_s4/fig_s4%s_%s.png", letter, mp), p,
         width = 170, height = 90, units = "mm", dpi = 300)
  cat("  saved\n")
}

# Combined 2×2 wide enough that left-side term names aren't cramped
if (length(panels) == 4) {
  combo <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 11, face = "bold",
                                  family = my_font))
  ggsave("../fig_s4/FigureS4.pdf", combo,
         width = 380, height = 220, units = "mm")
  ggsave("../fig_s4/FigureS4.png", combo,
         width = 380, height = 220, units = "mm", dpi = 300)
  cat("   FigureS4.pdf (combined 2x2)\n")
}

if (!.use_cache) {
  all_res_dt <- rbindlist(all_results, idcol = "MP")
  fwrite(all_res_dt, "../fig_s4/enrichment_results_all_MPs.csv")
}
cat("full results table saved\n")

message("=== Supp Fig 4 done ===")
