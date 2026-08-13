#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(dplyr); library(grid)
  library(showtext); library(ggrepel); library(ggnewscale)
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

theme_pub <- function(base_size = 9) {
  theme_classic(base_family = FAM, base_size = base_size) +
    theme(axis.text  = element_text(color = "black"),
          axis.line  = element_line(linewidth = 0.4, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"))
}
umap_arrow_axes <- function(data, x_col, y_col, frac = 0.16,
                            label_x = "Component 1", label_y = "Component 2",
                            text_size = 2.4, line_size = 0.3, arrow_mm = 1.2,
                            inset_frac = 0.05) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr); y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr);          y1 <- y0 + frac * diff(yr)
  arr <- arrow(length = unit(arrow_mm, "mm"), ends = "last", type = "closed")
  list(
    annotate("segment", x = x0, xend = x1, y = y0, yend = y0,
             arrow = arr, linewidth = line_size, color = "black"),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y1,
             arrow = arr, linewidth = line_size, color = "black"),
    annotate("text", x = (x0 + x1)/2, y = y0, label = label_x,
             vjust = 2.4, size = text_size, family = FAM),
    annotate("text", x = x0, y = (y0 + y1)/2, label = label_y,
             angle = 90, vjust = -1.4, size = text_size, family = FAM),
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(),
          panel.grid = element_blank())
  )
}

# Fig 3B  trajectory: subset arrows + corner Component arrows back
pt_data <- as.data.frame(fread("../fig2/pseudotime_umap.csv.gz"))
seg <- as.data.frame(fread("../fig2/monocle3_graph_segments_oriented.csv"))
mark <- as.data.frame(fread("../fig2/monocle3_graph_root_tip.csv"))

set.seed(42)
n_plot <- min(nrow(pt_data), 20000)
pt_sub <- pt_data[sample(nrow(pt_data), n_plot), ]
pt_sub <- pt_sub[pt_sub$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
pt_sub <- pt_sub[sample(nrow(pt_sub)), ]

# Pick 7 representative arrows: longest segments stratified by pt_mean bins,
# so the arrowheads scatter across the pt range without crowding.
seg$len <- sqrt((seg$x_end - seg$x_start)^2 + (seg$y_end - seg$y_start)^2)
seg$bin <- cut(seg$pt_mean, breaks = 7, include.lowest = TRUE)
arrow_seg <- seg %>% group_by(bin) %>% slice_max(len, n = 1, with_ties = FALSE) %>% ungroup()

pt_lim <- range(c(seg$pt_start, seg$pt_end), na.rm = TRUE)
arr <- arrow(length = unit(2.0, "mm"), type = "closed", ends = "last")

p3b <- ggplot() +
  geom_point(data = pt_sub,
             aes(UMAP1, UMAP2, color = dominant_MP),
             size = 0.55, alpha = 0.45, stroke = 0, shape = 16) +
  scale_color_manual(values = mp_colors, labels = mp_lbl,
                     name = "Dominant MP",
                     guide = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  ggnewscale::new_scale_color() +
  # All 59 segments as plain lines (no arrows) coloured by pseudotime
  geom_segment(data = seg,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                   color = pt_mean),
               linewidth = 0.85, lineend = "round") +
  # Only 7 representative arrows scattered across pt bins
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
                  segment.size = 0.3, box.padding = 0.4, point.padding = 0.4,
                  max.overlaps = Inf) +
  scale_shape_manual(values = c(root = 23, tip = 21), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Malignant cell trajectory (Monocle3 principal graph)") +
  coord_equal() +
  umap_arrow_axes(pt_sub, "UMAP1", "UMAP2",
                  label_x = "Component 1", label_y = "Component 2") +
  theme_pub(9) +
  theme(plot.title       = element_text(size = 9, face = "bold"),
        legend.position  = "right",
        legend.title     = element_text(size = 7, face = "bold"),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"),
        legend.box       = "vertical",
        legend.spacing.y = unit(2, "mm"))

ggsave("fig3b_trajectory.pdf", p3b, width = 5.6, height = 4.5)
ggsave("fig3b_trajectory.png", p3b, width = 5.6, height = 4.5, dpi = 300)
cat("Fig 3B: 7 arrows + corner Component arrows restored\n")

# Fig 3C  GeneSwitches scatter — redesigned
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
                      "F2","FGA","FGB","C3","CFH","VTN","AMBP","TF",
                      "APOC2","APOC3","APOB","APOH","ORM1",
                      "KNG1","AGX1","SERPINC1","HRG","A1BG","HPX",
                      "CPS1","ARG1","SLC2A2","SULT2A1","CYP2D6","HNF4A","NR1I3")

gs$gene_type <- "Other"
gs$gene_type[gs$gene %in% known_tfs] <- "TF"
gs$gene_type[gs$gene %in% surface_proteins] <- "Surface"
gs$gene_type <- factor(gs$gene_type, levels = c("Other","TF","Surface"))

# percentile x
gs$switch_pct <- 100 * (rank(gs$switch_pseudotime_rank, ties.method = "average") /
                          length(gs$switch_pseudotime_rank))

# label set: high R^2 + curated key genes
key_genes <- c("EPCAM","NAPSA","HOPX","AQP1","SFTPC","SFTPB","NKX2-1",
               "VIM","FN1","CDH2","CD44","MET","AXL","CXCR4",
               "KLF5","KLF6","SOX2","SOX4","IRF9","TEAD1",
               "TOP2A","MKI67","STMN1","ATF3","FOS","FOSB","JUN","JUNB",
               "SRGN","MCL1","NFKBIA")
top_r2 <- gs %>% slice_max(mcfadden_R2, n = 20) %>% pull(gene)
label_set <- unique(c(top_r2, intersect(gs$gene, key_genes)))
# cap at 35 to keep the plot readable
if (length(label_set) > 35) {
  label_set <- gs %>% filter(gene %in% label_set) %>%
    slice_max(mcfadden_R2, n = 35) %>% pull(gene)
}
gs$label <- ifelse(gs$gene %in% label_set, gs$gene, "")

# Draw order: background first, highlights last
bg <- gs %>% filter(gene_type == "Other")
hl <- gs %>% filter(gene_type != "Other")

# Custom palette: keep TF / Surface vivid; shrink "Other" to a desaturated grey
type_fill <- c(Other = "grey80", TF = "#00A087", Surface = "#E64B35")
type_col  <- c(Other = "grey60", TF = "#00744F", Surface = "#A4291C")

p3c <- ggplot() +
  # 1) faded background: tiny grey points for the bulk of "Other" genes
  geom_point(data = bg, aes(x = switch_pct, y = mcfadden_R2),
             size = 0.45, alpha = 0.18, color = "grey55", shape = 16) +
  # 2) density contours over Other-gene cloud (extra structure cue)
  geom_density_2d(data = bg, aes(x = switch_pct, y = mcfadden_R2),
                   color = "grey45", linewidth = 0.22, alpha = 0.5,
                   bins = 5) +
  # 3) overall trend line for the entire gene set
  geom_smooth(data = gs, aes(x = switch_pct, y = mcfadden_R2),
              method = "loess", se = FALSE,
              color = "grey25", linewidth = 0.5, linetype = "dashed",
              span = 0.35) +
  # 4) highlight TF / Surface points: filled circles with dark stroke
  geom_point(data = hl,
             aes(x = switch_pct, y = mcfadden_R2,
                 fill = gene_type, color = gene_type),
             shape = 21, size = 2.0, stroke = 0.4, alpha = 0.95) +
  scale_fill_manual(values = type_fill[c("TF","Surface")],
                    name = NULL,
                    labels = c(TF = "Transcription factors",
                               Surface = "Surface proteins")) +
  scale_color_manual(values = type_col[c("TF","Surface")], guide = "none") +
  # 4) gene labels — bold italic, slightly larger, with halo via bg.color
  ggrepel::geom_text_repel(
    data = gs %>% filter(label != ""),
    aes(x = switch_pct, y = mcfadden_R2, label = label),
    size = 2.4, family = FAM, fontface = "bold.italic", color = "black",
    bg.color = "white", bg.r = 0.18,
    box.padding = 0.30, point.padding = 0.20,
    segment.size = 0.22, segment.alpha = 0.55,
    min.segment.length = 0.05, max.overlaps = Inf
  ) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%","25%","50%","75%","100%"),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Pseudo-timeline (percentile)",
       y = expression("Quality of fitting (McFadden " * R^2 * ")"),
       title = "GeneSwitches analysis") +
  theme_pub(9) +
  theme(plot.title       = element_text(size = 9, face = "bold"),
        legend.position  = "right",
        legend.title     = element_text(size = 7, face = "bold"),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"),
        legend.box       = "vertical",
        legend.spacing.y = unit(2, "mm"),
        panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3))

ggsave("fig3c_geneswitches.pdf", p3c, width = 5.4, height = 4.0)
ggsave("fig3c_geneswitches.png", p3c, width = 5.4, height = 4.0, dpi = 300)
cat("Fig 3C redesigned: hex density background + highlighted TF/Surface points\n")
