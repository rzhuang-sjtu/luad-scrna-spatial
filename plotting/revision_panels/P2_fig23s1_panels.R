##############################################################################
## Panels added to or replacing parts of Fig. 2, Fig. 3 and Fig. S1.
##
## Fonts and theme (theme_pub, Arial via showtext, mp_colors) are taken
## unchanged from the Fig. 3 rendering script so the new panels match the
## published ones. Only the input data differ: they come from the analyses
## recomputed during revision.
##
## Outputs, each as .pdf and .png at 400 dpi:
##   Fig2_new__pseudotime_composition   dominant-programme composition along
##                                      pseudotime (Fig. 2K)
##   Fig3_new__mp3_model_ladder         MP3 hazard ratio across four nested
##                                      models (Fig. 3K)
##   FigS1_alt__kbet_by_celltype        kBET under its original definition,
##                                      by cell type (Fig. S1E)
##   FigS1_alt__lisi                    iLISI and cLISI before and after
##                                      integration (Fig. S1F)
##
## Usage: Rscript P2_fig23s1_panels.R
##############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
suppressPackageStartupMessages({
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

# Arial via showtext, as in the original plotting scripts
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

theme_pub <- function(base_size = 8) {
  theme_classic(base_family = my_font, base_size = base_size) +
    theme(axis.text       = element_text(color = "black"),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          legend.key.size = unit(3, "mm"),
          legend.title    = element_text(size = rel(0.95), face = "bold"),
          legend.text     = element_text(size = rel(0.85)),
          legend.background = element_blank(),
          plot.margin     = margin(2, 2, 2, 2),
          strip.text      = element_text(face = "bold", size = rel(1)),
          plot.title      = element_text(face = "bold", size = rel(1.1)))
}

RES <- "${PROJECT_ROOT}/results"
OUT <- "${WORK_ROOT}/revision_panels"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PNG_DPI <- 400

save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = grDevices::cairo_pdf)
  ## showtext draws glyphs at the dpi it was told, while ggsave rasterises at
  ## its own. Left at 300 against a 400-dpi PNG the text came out about a tenth
  ## too small, so the PNG and the PDF of the same panel did not match. The
  ## vector output is unaffected (cairo_pdf writes glyph outlines at true size),
  ## so only the raster save needs the dpi handed to it.
  if (requireNamespace("showtext", quietly = TRUE) && my_font != "sans") {
    showtext::showtext_opts(dpi = PNG_DPI)
    on.exit(showtext::showtext_opts(dpi = 300), add = TRUE)
  }
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = PNG_DPI)
  cat("  ✓", name, "\n")
}

## Fail loudly on an empty selection. ggplot renders an empty data frame
## without complaint, so a stale filter produces a blank panel that still
## reports success. This caught a filter left behind when the result tables
## were renamed to English.
need_rows <- function(d, what) {
  if (nrow(d) == 0) stop("no rows selected for ", what,
                         " - check the filter against the source table")
  invisible(d)
}

## 1. New Fig. 2K: dominant-programme composition along pseudotime.
## The published (H) plots smoothed per-programme scores, from which the
## identity of the trajectory endpoint cannot be read.
CP <- read.csv(file.path(RES, "dpt_root_sensitivity/composition_by_ventile.csv"),
               check.names = FALSE)
names(CP)[1] <- "ventile"
long <- CP %>%
  select(ventile, MP1, MP2, MP3, MP4) %>%
  pivot_longer(-ventile, names_to = "MP", values_to = "pct") %>%
  mutate(MP = factor(MP, levels = c("MP1", "MP2", "MP3", "MP4")))

need_rows(long, "pseudotime composition")
p1 <- ggplot(long, aes(ventile, pct, fill = MP)) +
  geom_area(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = mp_colors, breaks = c("MP1","MP2","MP3","MP4")) +
  scale_x_continuous(expand = c(0, 0)) +
  ## Leave a hairline of headroom: in a 100% stacked area the top band's
  ## white stroke is clipped by the panel edge and reads as truncated.
  ## Expand the upper limit only; the lower one stays flush with zero.
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02)),
                     breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Pseudotime ventile", y = "Cells by dominant program (%)",
       fill = NULL) +
  theme_pub() +
  theme(legend.position = "right")
save_panel(p1, "Fig2_new__pseudotime_composition", 3.0, 2.0)

## 2. New Fig. 3K: MP3 across nested Cox models
L <- read.csv(file.path(RES, "results23_canonical/mp3_ladder.csv"),
              check.names = FALSE)
keep <- c("MP3 alone", "MP3 + age, stage, sex", "MP3 + the other three MPs",
          "MP3 + three MPs + clinical (published model)")
lab  <- c("MP3 alone", "+ clinical", "+ other MPs", "+ MPs + clinical")
D <- L[match(keep, L$model), ]
D$lab <- factor(lab, levels = lab)
D$sig <- ifelse(D$p < 0.05, "yes", "no")

need_rows(D, "MP3 nested models")
p2 <- ggplot(D, aes(lab, HR, color = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3,
             color = "grey40") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_text(aes(y = hi, label = sprintf("%.2f", HR)), vjust = -0.9,
            size = 2.1, family = my_font, show.legend = FALSE) +
  scale_color_manual(values = c(yes = "#00A087", no = "grey55"),
                     guide = "none") +
  scale_y_continuous(limits = c(0.85, 2.15)) +
  ## Two lines: on one line this title is as long as the panel is tall, and
  ## ggplot neither shrinks nor wraps it, so the closing "I)" was being clipped
  ## off the top of the vector output.
  labs(x = NULL, y = "MP3 hazard ratio\nper 1 SD (95% CI)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_panel(p2, "Fig3_new__mp3_model_ladder", 2.6, 2.1)

## 3. Replacement Fig. S1E: kBET under its original definition, by cell type
KB <- read.csv(file.path(RES, "batch_metrics/kbet.csv"), check.names = FALSE)
names(KB) <- c("level", "celltype", "embedding", "kbet_after", "n",
               "kbet_before", "n_batch")
ct <- need_rows(KB[KB$level == "by_cell_type", ], "kBET by cell type")
ct <- ct[order(ct$kbet_after), ]
ct$celltype <- factor(ct$celltype, levels = ct$celltype)
disp <- c(T_NK = "T/NK", B = "B cell", Plasma = "Plasma cell")
levels(ct$celltype) <- ifelse(levels(ct$celltype) %in% names(disp),
                              disp[levels(ct$celltype)],
                              levels(ct$celltype))

p3 <- ggplot(ct) +
  geom_segment(aes(x = celltype, xend = celltype,
                   y = kbet_before, yend = kbet_after),
               color = "grey80", linewidth = 0.5) +
  geom_point(aes(celltype, kbet_before), color = "grey55", size = 1.4) +
  geom_point(aes(celltype, kbet_after), color = "#3C5488", size = 1.8) +
  scale_y_continuous(limits = c(0.88, 1.005)) +
  labs(x = NULL, y = "kBET rejection rate") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_panel(p3, "FigS1_alt__kbet_by_celltype", 2.8, 2.0)

## 4. Replacement Fig. S1F: iLISI and cLISI
LI <- read.csv(file.path(RES, "batch_metrics/lisi.csv"), check.names = FALSE)
names(LI) <- c("embedding", "iLISI", "iLISI_max", "cLISI", "cLISI_max")
M <- data.frame(
  metric = factor(c("iLISI\n(batch mixing)", "cLISI\n(cell-type purity)"),
                  levels = c("cLISI\n(cell-type purity)",
                             "iLISI\n(batch mixing)")),
  before = c(LI$iLISI[1], LI$cLISI[1]),
  after  = c(LI$iLISI[2], LI$cLISI[2]))

p4 <- ggplot(M) +
  geom_segment(aes(y = metric, yend = metric, x = before, xend = after),
               color = "grey80", linewidth = 0.5) +
  geom_point(aes(before, metric), color = "grey55", size = 1.5) +
  geom_point(aes(after, metric), color = "#3C5488", size = 2.0) +
  geom_text(aes(after, metric, label = sprintf("%.2f", after)),
            hjust = -0.35, size = 2.1, family = my_font, color = "#3C5488") +
  scale_x_continuous(limits = c(0.7, 2.6)) +
  labs(x = "LISI", y = NULL) +
  theme_pub()
save_panel(p4, "FigS1_alt__lisi", 2.4, 1.3)

cat("\nwritten to", OUT, "\n")
