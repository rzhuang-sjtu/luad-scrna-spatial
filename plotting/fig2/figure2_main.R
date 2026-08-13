## Figure 2A–E  |  LUAD Meta-Program  (v4)

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)
library(patchwork)
suppressPackageStartupMessages({
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
  if (requireNamespace("ggrepel",  quietly = TRUE)) library(ggrepel)
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

# ── UMAP corner-arrow helper ──
umap_arrow_axes <- function(data, x_col, y_col,
                            frac = 0.16, label_x = "UMAP 1", label_y = "UMAP 2",
                            text_size = 2.4, line_size = 0.3, arrow_mm = 1.2,
                            inset_frac = 0.02) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr); y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr);          y1 <- y0 + frac * diff(yr)
  arr <- grid::arrow(length = grid::unit(arrow_mm, "mm"), ends = "last", type = "closed")
  fam <- if (exists("my_font")) my_font else "sans"
  list(
    annotate("segment", x=x0, xend=x1, y=y0, yend=y0, arrow=arr, linewidth=line_size, color="black"),
    annotate("segment", x=x0, xend=x0, y=y0, yend=y1, arrow=arr, linewidth=line_size, color="black"),
    annotate("text", x=(x0+x1)/2, y=y0, label=label_x, vjust=2.4, size=text_size, family=fam),
    annotate("text", x=x0, y=(y0+y1)/2, label=label_y, angle=90, vjust=-1.4, size=text_size, family=fam),
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          panel.grid=element_blank())
  )
}

mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F","Unassigned"="grey80")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

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


# 2A GEP 770 GEP

cat("── Fig 2A ──\n")

corr <- fread("gep_spearman_corr.csv")
gep_ids <- corr[[1]]
corr_mat <- as.matrix(corr[, -1, with=FALSE])
rownames(corr_mat) <- gep_ids
colnames(corr_mat) <- gep_ids

dist_mat <- as.dist(1 - corr_mat)
hc <- hclust(dist_mat, method="average")

anno_df <- fread("gep_mp_annotation.csv")
anno_df <- anno_df[match(gep_ids, anno_df$gep_id)]

# MP block
dend_order <- hc$order
mp_ordered <- anno_df$MP[dend_order]

get_mp_blocks <- function(mp_vec) {
  rle_res <- rle(mp_vec)
  ends <- cumsum(rle_res$lengths)
  starts <- c(1, ends[-length(ends)] + 1)
  data.frame(MP=rle_res$values, start=starts, end=ends, stringsAsFactors=FALSE)
}
mp_blocks <- get_mp_blocks(mp_ordered)
mp_blocks <- mp_blocks[mp_blocks$MP %in% c("MP1","MP2","MP3","MP4"), ]

# + MP annotation
ha_left <- rowAnnotation(
  MetaProgram = anno_df$MP,
  col = list(MetaProgram = mp_colors),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize=8, fontfamily=my_font),
  annotation_legend_param = list(title_gp = gpar(fontsize=9, fontfamily=my_font),
                                 labels_gp = gpar(fontsize=8, fontfamily=my_font)),
  width = unit(4, "mm")
)

col_fun <- colorRamp2(c(-0.5, 0, 0.5, 1), c("#3C5488","white","#E64B35","#8B0000"))

ht_name <- "Spearman\nCorrelation"

draw_2a <- function() {
  ht <- Heatmap(corr_mat,
                name = ht_name,
                col = col_fun,
                cluster_rows = hc,
                cluster_columns = hc,
                show_row_names = FALSE,       # 770 
                show_column_names = FALSE,     # 770 
                left_annotation = ha_left,
                row_dend_width = unit(15, "mm"),
                column_dend_height = unit(15, "mm"),
                use_raster = TRUE,             # 
                raster_quality = 3,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=9, fontfamily=my_font),
                  labels_gp = gpar(fontsize=8, fontfamily=my_font),
                  legend_height = unit(25, "mm")
                ))
  draw(ht, merge_legend=TRUE)
  
  for(i in seq_len(nrow(mp_blocks))) {
    s <- mp_blocks$start[i]; e <- mp_blocks$end[i]; mp <- mp_blocks$MP[i]
    decorate_heatmap_body(ht_name, {
      n <- ncol(corr_mat)
      x1 <- (s-1)/n; x2 <- e/n
      y1 <- 1 - e/n; y2 <- 1 - (s-1)/n
      grid.rect(x=unit(x1,"npc"), y=unit(y1,"npc"),
                width=unit(x2-x1,"npc"), height=unit(y2-y1,"npc"),
                just=c("left","bottom"),
                gp=gpar(col=mp_colors[mp], fill=NA, lwd=2))
      grid.text(mp, x=unit((x1+x2)/2,"npc"), y=unit((y1+y2)/2,"npc"),
                gp=gpar(col=mp_colors[mp], fontsize=10, fontfamily=my_font, fontface="bold"))
    })
  }
}

pdf("fig2a_gep_heatmap.pdf", width=7, height=6)
draw_2a()
dev.off()

png("fig2a_gep_heatmap.png", width=7, height=6, units="in", res=300)
draw_2a()
dev.off()
cat("  2A saved\n")


## 2B  MP score UMAP × 4

cat("── Fig 2B ──\n")

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)
umap <- fread("malignant_umap_metadata.csv.gz")

plot_mp_umap <- function(dt, score_col, title, color_high) {
  dt_sorted <- dt[order(get(score_col))]
  v <- dt_sorted[[score_col]]
  # Map score to per-point alpha so low-score cells fade into the grey
  # background and high-score cells stay vivid. Clamp at the 2/98 % tails
  # to avoid a single outlier flattening the dynamic range.
  lo <- as.numeric(quantile(v, 0.02, na.rm=TRUE))
  hi <- as.numeric(quantile(v, 0.98, na.rm=TRUE))
  rng <- max(hi - lo, 1e-9)
  a_norm <- pmin(1, pmax(0, (v - lo) / rng))
  dt_sorted[, alpha_pt := 0.15 + 0.80 * a_norm]
  p <- ggplot(dt_sorted, aes(x=UMAP1, y=UMAP2, color=.data[[score_col]], alpha=alpha_pt))
  if(has_ggrastr) {
    p <- p + ggrastr::rasterise(geom_point(size=0.18, stroke=0, shape=16), dpi=600)
  } else {
    p <- p + geom_point(size=0.18, stroke=0, shape=16)
  }
  # Non-linear ramp: bottom ~35 % of the score range stays grey, then
  # quickly transitions to the highlight colour. This stops low cells
  # from looking falsely "expressed".
  p + scale_color_gradientn(
        colors = c("grey92", "grey88", color_high),
        values = scales::rescale(c(0, 0.35, 1)),
        name = "Score",
        guide = guide_colorbar(barwidth = unit(2.5,"mm"),
                               barheight = unit(15,"mm"),
                               frame.colour = "black",
                               frame.linewidth = 0.3)) +
    scale_alpha_identity() +
    labs(title=title) +
    coord_equal(clip = "off") +
    theme_void(base_family=my_font) +
    theme(plot.title=element_text(size=8, face="bold", hjust=0),
          legend.position="right",
          legend.title=element_text(size=6),
          legend.text=element_text(size=5),
          plot.margin=margin(1,1,1,1)) +
    umap_arrow_axes(dt_sorted, "UMAP1", "UMAP2", inset_frac = 0.05)
}

p2b_1 <- plot_mp_umap(umap, "MP1_score", "MP1", color_high="#E64B35")
p2b_2 <- plot_mp_umap(umap, "MP2_score", "MP2", color_high="#4DBBD5")
p2b_3 <- plot_mp_umap(umap, "MP3_score", "MP3", color_high="#00A087")
p2b_4 <- plot_mp_umap(umap, "MP4_score", "MP4", color_high="#3C5488")

# Tight 2x2: kill inter-panel spacing via patchwork.
p2b <- ((p2b_1 | p2b_2) / (p2b_3 | p2b_4)) +
       plot_layout(widths = c(1, 1), heights = c(1, 1)) &
       theme(plot.margin = margin(1, 1, 1, 1))

ggsave("fig2b_mp_score_umap.pdf", p2b, width=125, height=115, units="mm")
ggsave("fig2b_mp_score_umap.png", p2b, width=125, height=115, units="mm", dpi=600)
cat("  2B saved\n")


# 2C Dominant MP UMAP high-confidence cells only
##   Keep cells where dominant MP score is meaningfully higher than
##   the runner-up (gap > delta). Cells with mixed program activity
##   (small gap) are plotted as light grey "ambiguous" so they don't
##   bleed colour into other MP territories.

cat("── Fig 2C ──\n")

umap_mp <- umap[dominant_MP %in% c("MP1","MP2","MP3","MP4")]
umap_mp[, dominant_MP := factor(dominant_MP, levels=c("MP1","MP2","MP3","MP4"))]

# Compute second-highest MP score per cell, then confidence gap
score_cols <- c("MP1_score","MP2_score","MP3_score","MP4_score")
score_mat <- as.matrix(umap_mp[, ..score_cols])
sorted    <- t(apply(score_mat, 1, sort, decreasing=TRUE))
umap_mp[, top1_score := sorted[, 1]]
umap_mp[, top2_score := sorted[, 2]]
umap_mp[, mp_gap := top1_score - top2_score]

# Aggressive filter: keep top 50% by mp_gap (drops the more ambiguous 50%).
# Drop ambiguous completely (no grey background) so visible cells are
# unambiguously a single MP. Add 2D density contours per MP so each
# territory is outlined regardless of cell density.
delta <- as.numeric(quantile(umap_mp$mp_gap, 0.50, na.rm = TRUE))
umap_mp[, confident := mp_gap >= delta]
n_total <- nrow(umap_mp); n_conf <- sum(umap_mp$confident)
cat(sprintf("  gap delta (50th pct) = %.3f, kept %d / %d cells (%.1f%%) as high-confidence\n",
            delta, n_conf, n_total, 100 * n_conf / n_total))

umap_conf <- umap_mp[confident == TRUE]

# Plot order: MPs with fewer cells drawn last (on top) so they aren't buried
mp_freq <- umap_conf[, .N, by = dominant_MP][order(-N)]
umap_conf[, dominant_MP := factor(dominant_MP, levels = mp_freq$dominant_MP)]
umap_conf <- umap_conf[order(dominant_MP, decreasing = TRUE)]
umap_conf[, dominant_MP := factor(dominant_MP, levels = c("MP1","MP2","MP3","MP4"))]

p2c <- ggplot(umap_conf, aes(x = UMAP1, y = UMAP2, color = dominant_MP)) +
  ggrastr::rasterise(
    geom_point(size = 0.55, stroke = 0, alpha = 0.5, shape = 16),
    dpi = 600
  ) +
  scale_color_manual(values=mp_colors, labels=mp_labels, name="Dominant MP") +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  coord_equal() +
  theme_pub(base_size=9) +
  theme(legend.position="right") +
  labs(x = NULL, y = NULL,
       caption = sprintf("%d high-confidence cells (top 50%% by MP-score gap).",
                         n_conf)) +
  umap_arrow_axes(umap_conf, "UMAP1", "UMAP2")

ggsave("fig2c_dominant_mp_umap.pdf", p2c, width=130, height=90, units="mm")
ggsave("fig2c_dominant_mp_umap.png", p2c, width=130, height=90, units="mm", dpi=600)
cat("  2C saved\n")


# 2C v2 UMAP recomputed in 4-D MP score space
##   Runs `compute_mp_space_umap.R` upstream to produce
##   malignant_mp_umap_metadata.csv.gz with MP_UMAP1/2 columns.
##   In this space, MPs separate by construction: each MP "axis" of
##   the 4-D score simplex projects to its own region of 2-D UMAP.
mp_umap_path <- "malignant_mp_umap_metadata.csv.gz"
if (file.exists(mp_umap_path)) {
  cat("── Fig 2C (MP-space UMAP variant) ──\n")
  umap2 <- fread(mp_umap_path)
  umap2 <- umap2[dominant_MP %in% c("MP1","MP2","MP3","MP4")]
  umap2[, dominant_MP := factor(dominant_MP, levels = c("MP1","MP2","MP3","MP4"))]
  set.seed(42); umap2 <- umap2[sample(.N)]

  p2c_mp <- ggplot(umap2, aes(MP_UMAP1, MP_UMAP2, color = dominant_MP)) +
    ggrastr::rasterise(
      geom_point(size = 0.45, stroke = 0, alpha = 0.55, shape = 16),
      dpi = 600
    ) +
    scale_color_manual(values = mp_colors, labels = mp_labels,
                       name = "Dominant MP") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal() +
    theme_pub(base_size = 9) +
    theme(legend.position = "right") +
    labs(x = NULL, y = NULL) +
    umap_arrow_axes(umap2, "MP_UMAP1", "MP_UMAP2")

  ggsave("fig2c_mp_umap.pdf", p2c_mp, width = 130, height = 90, units = "mm")
  ggsave("fig2c_mp_umap.png", p2c_mp, width = 130, height = 90, units = "mm",
         dpi = 600)
  cat("  2C MP-space saved\n")
} else {
  cat("  [skip] fig2c_mp_umap   run compute_mp_space_umap.R first\n")
}


# 2D MP proportion by tissue scale +

cat("── Fig 2D ──\n")

prop <- fread("mp_proportion_by_tissue.csv")
prop_long <- melt(prop, id.vars="tissue_type", variable.name="MP", value.name="pct")

max_pct <- max(prop_long$pct, na.rm=TRUE)
if(max_pct <= 1.01) {
  prop_long[, pct := pct * 100]
}

prop_long[, MP := factor(MP, levels=c("MP4","MP3","MP2","MP1"))]

# tissue name
tissue_abbrev <- c(
  "Primary_Tumor"="Tumor", "Primary"="Tumor",
  "LN_Metastasis"="LN_Met", "LN_Met"="LN_Met",
  "Distant_Metastasis"="Dist_Met", "Distant_Met"="Dist_Met",
  "Brain_Metastasis"="Brain_Met",
  "Pleural_Effusion"="Pleural", "Pleural"="Pleural",
  "Precancerous"="Precanc",
  "Normal"="Normal", "Adjacent"="Adjacent"
)

prop_long[, tissue_label := ifelse(tissue_type %in% names(tissue_abbrev),
                                   tissue_abbrev[tissue_type],
                                   tissue_type)]

tissue_rank <- c("Normal","Adjacent","Precanc","Tumor","LN_Met",
                 "Dist_Met","Brain_Met","Pleural")
tissue_present <- intersect(tissue_rank, unique(prop_long$tissue_label))
if(length(tissue_present) == 0) tissue_present <- unique(prop_long$tissue_label)
prop_long[, tissue_label := factor(tissue_label, levels=tissue_present)]

p2d <- ggplot(prop_long, aes(x=tissue_label, y=pct, fill=MP)) +
  geom_bar(stat="identity", width=0.75, color="white", linewidth=0.1) +
  scale_fill_manual(values=mp_colors, labels=mp_labels, name="Cell Type") +
  scale_y_continuous(labels=function(x) paste0(x,"%"),
                     expand=expansion(mult=c(0, 0.02)),
                     limits=c(0, 100)) +
  labs(x=NULL, y="Cell percent ratio") +
  theme_pub(base_size=8) +
  theme(axis.text.x=element_text(size=5.5, angle=45, hjust=1, face="italic"),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=7),
        legend.text=element_text(size=6),
        legend.title=element_blank(),
        legend.key.size=unit(3,"mm"))

ggsave("fig2d_mp_by_tissue.pdf", p2d, width=65, height=65, units="mm")
ggsave("fig2d_mp_by_tissue.png", p2d, width=65, height=65, units="mm", dpi=300)
cat("  2D saved\n")


# 2E Dot plot + z-score + top 12/MP + Fig1B

cat("── Fig 2E ──\n")

dot <- fread("mp_dotplot_markers.csv")
gene_order_df <- fread("mp_marker_gene_order.csv")

# MP top 12
gene_order_df <- gene_order_df[order(factor(MP, levels=c("MP1","MP2","MP3","MP4")), rank)]
gene_top <- gene_order_df[rank <= 5]
dot_sub <- dot[gene %in% gene_top$gene]

# z-score mean_expression per gene ( 4 MP
dot_sub[, z_expr := scale(mean_expression), by=gene]
dot_sub[is.nan(z_expr), z_expr := 0]

dot_sub[, gene := factor(gene, levels=gene_top$gene)]
dot_sub[, MP := factor(MP, levels=rev(c("MP1","MP2","MP3","MP4")))]

# frac_expressing
max_frac <- max(dot_sub$frac_expressing, na.rm=TRUE)
if(max_frac <= 1.01) {
  dot_sub[, pct_expr := frac_expressing * 100]
} else {
  dot_sub[, pct_expr := frac_expressing]
}

p2e <- ggplot(dot_sub, aes(x=gene, y=MP)) +
  geom_point(aes(size=pct_expr, fill=z_expr),
             shape=21, color="black", stroke=0.25) +
  scale_size_continuous(
    name="Fraction of cells\nin group (%)",
    range=c(0.5, 6),
    breaks=c(25, 50, 75),
    limits=c(0, 100)
  ) +
  scale_fill_gradient2(
    name="Mean expression\nin group",
    low="#3C5488", mid="#F7F7F7", high="#E64B35",
    midpoint=0,
    limits=c(-2, 2),
    oob=squish,
    guide=guide_colorbar(barwidth=unit(3.5,"mm"), barheight=unit(20,"mm"),
                         frame.colour="black", frame.linewidth=0.3,
                         ticks.colour="black")
  ) +
  labs(x=NULL, y=NULL) +
  theme_pub(base_size=7) +
  theme(
    axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=5.5, face="italic"),
    axis.text.y=element_text(size=7),
    legend.position="right",
    legend.box="vertical",
    legend.spacing.y=unit(4,"mm"),
    panel.border=element_rect(color="black", fill=NA, linewidth=0.4),
    panel.grid=element_blank(),
    axis.line=element_blank(),
    axis.ticks=element_line(linewidth=0.3)
  )

ggsave("fig2e_dotplot_markers.pdf", p2e, width=180, height=80, units="mm")
ggsave("fig2e_dotplot_markers.png", p2e, width=180, height=80, units="mm", dpi=300)
cat("  2E saved\n")

message("=== Fig 2A-E v4 done ===")

# Figure 2A |

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)

mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F","Unassigned"="grey80")

corr <- fread("gep_spearman_corr.csv")
gep_ids <- corr[[1]]
corr_mat <- as.matrix(corr[, -1, with=FALSE])
rownames(corr_mat) <- gep_ids
colnames(corr_mat) <- gep_ids

dist_mat <- as.dist(1 - corr_mat)
hc <- hclust(dist_mat, method="average")

anno_df <- fread("gep_mp_annotation.csv")
anno_df <- anno_df[match(gep_ids, anno_df$gep_id)]

dend_order <- hc$order
mp_ordered <- anno_df$MP[dend_order]

get_mp_blocks <- function(mp_vec) {
  rle_res <- rle(mp_vec)
  ends <- cumsum(rle_res$lengths)
  starts <- c(1, ends[-length(ends)] + 1)
  data.frame(MP=rle_res$values, start=starts, end=ends, stringsAsFactors=FALSE)
}
mp_blocks <- get_mp_blocks(mp_ordered)
mp_blocks <- mp_blocks[mp_blocks$MP %in% c("MP1","MP2","MP3","MP4"), ]

ha_left <- rowAnnotation(
  MetaProgram = anno_df$MP,
  col = list(MetaProgram = mp_colors),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize=8, fontfamily=my_font),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontfamily=my_font),
    labels_gp = gpar(fontsize=8, fontfamily=my_font)
  ),
  width = unit(4, "mm")
)

col_fun <- colorRamp2(c(-0.5, 0, 0.5, 1), c("#3C5488","white","#E64B35","#8B0000"))
ht_name <- "Spearman\nCorrelation"

draw_mp_borders <- function() {
  for(i in seq_len(nrow(mp_blocks))) {
    s <- mp_blocks$start[i]; e <- mp_blocks$end[i]; mp <- mp_blocks$MP[i]
    decorate_heatmap_body(ht_name, {
      n <- ncol(corr_mat)
      x1 <- (s-1)/n; x2 <- e/n
      y1 <- 1 - e/n; y2 <- 1 - (s-1)/n
      grid.rect(x=unit(x1,"npc"), y=unit(y1,"npc"),
                width=unit(x2-x1,"npc"), height=unit(y2-y1,"npc"),
                just=c("left","bottom"),
                gp=gpar(col=mp_colors[mp], fill=NA, lwd=2))
      grid.text(mp, x=unit((x1+x2)/2,"npc"), y=unit((y1+y2)/2,"npc"),
                gp=gpar(col=mp_colors[mp], fontsize=10,
                        fontfamily=my_font, fontface="bold"))
    })
  }
}


# 1 dendrogram

draw_2a_v1 <- function() {
  ht <- Heatmap(corr_mat,
                name = ht_name,
                col = col_fun,
                cluster_rows = hc,
                cluster_columns = hc,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                left_annotation = ha_left,
                use_raster = TRUE,
                raster_quality = 3,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=9, fontfamily=my_font),
                  labels_gp = gpar(fontsize=8, fontfamily=my_font),
                  legend_height = unit(25, "mm")
                ))
  draw(ht, merge_legend=TRUE)
  draw_mp_borders()
}

pdf("fig2a_v1_no_dend.pdf", width=6.5, height=6)
draw_2a_v1()
dev.off()

png("fig2a_v1_no_dend.png", width=6.5, height=6, units="in", res=300)
draw_2a_v1()
dev.off()

cat("1 saved: fig2a_v1_no_dend\n")


# 2 dendrogram +

draw_2a_v2 <- function() {
  ht <- Heatmap(corr_mat,
                name = ht_name,
                col = col_fun,
                cluster_rows = hc,
                cluster_columns = hc,
                show_row_dend = TRUE,
                show_column_dend = TRUE,
                row_dend_width = unit(8, "mm"),
                column_dend_height = unit(8, "mm"),
                row_dend_gp = gpar(lwd=0.15),
                column_dend_gp = gpar(lwd=0.15),
                show_row_names = FALSE,
                show_column_names = FALSE,
                left_annotation = ha_left,
                use_raster = TRUE,
                raster_quality = 3,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=9, fontfamily=my_font),
                  labels_gp = gpar(fontsize=8, fontfamily=my_font),
                  legend_height = unit(25, "mm")
                ))
  draw(ht, merge_legend=TRUE)
  draw_mp_borders()
}

pdf("fig2a_v2_thin_dend.pdf", width=7, height=6.5)
draw_2a_v2()
dev.off()

png("fig2a_v2_thin_dend.png", width=7, height=6.5, units="in", res=300)
draw_2a_v2()
dev.off()

cat("2 saved: fig2a_v2_thin_dend\n")

message("=== Fig 2A  done ===")

## Figure 2F–J  |  LUAD Meta-Program

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)
library(patchwork)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)

mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F","Unassigned"="grey80")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}


# 2F Hallmark NES

cat("── Fig 2F ──\n")

nes <- fread("hallmark_nes_heatmap.csv")
fdr <- fread("hallmark_fdr_heatmap.csv")

# term
clean_term <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

terms <- nes$Term
nes_mat <- as.matrix(nes[, .(MP1, MP2, MP3, MP4)])
rownames(nes_mat) <- clean_term(terms)

fdr_mat <- as.matrix(fdr[, .(MP1, MP2, MP3, MP4)])

# NES MP
max_mp <- apply(nes_mat, 1, which.max)
row_ord <- order(max_mp, -apply(nes_mat, 1, max))
nes_mat <- nes_mat[row_ord, ]
fdr_mat <- fdr_mat[row_ord, ]

sig_mat <- matrix("", nrow(fdr_mat), ncol(fdr_mat))
sig_mat[fdr_mat < 0.05]  <- "*"
sig_mat[fdr_mat < 0.01]  <- "**"
sig_mat[fdr_mat < 0.001] <- "***"

col_nes <- colorRamp2(c(-3, 0, 3), c("#3C5488","white","#E64B35"))

draw_2f <- function() {
  ht <- Heatmap(nes_mat,
                name="NES",
                col=col_nes,
                cluster_rows=FALSE,
                cluster_columns=FALSE,
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_side="left",
                row_names_gp=gpar(fontsize=6, fontfamily=my_font),
                column_names_gp=gpar(fontsize=10, fontfamily=my_font),
                column_names_rot=0,
                rect_gp=gpar(col="white", lwd=0.5),
                cell_fun=function(j, i, x, y, w, h, fill) {
                  if(sig_mat[i, j] != "")
                    grid.text(sig_mat[i, j], x, y, gp=gpar(fontsize=5, fontfamily=my_font))
                },
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=9, fontfamily=my_font),
                  labels_gp=gpar(fontsize=8, fontfamily=my_font),
                  legend_height=unit(25,"mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig2f_hallmark_nes.pdf", width=5, height=7)
draw_2f()
dev.off()

png("fig2f_hallmark_nes.png", width=5, height=7, units="in", res=300)
draw_2f()
dev.off()
cat("  2F saved\n")


## 2G  Pseudotime UMAP

cat("── Fig 2G ──\n")

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)

pt <- fread("pseudotime_umap_winsorized.csv.gz")

p2g <- ggplot(pt[order(pt_winsorized)],
              aes(x=UMAP1, y=UMAP2, color=pt_winsorized))

if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}

p2g <- p2g +
  scale_color_viridis_c(option="inferno", name="Pseudotime",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(18,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  coord_equal() +
  theme_void(base_family=my_font) +
  theme(legend.position="right",
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.margin=margin(5,5,5,5))

ggsave("fig2g_pseudotime_umap.pdf", p2g, width=90, height=80, units="mm")
ggsave("fig2g_pseudotime_umap.png", p2g, width=90, height=80, units="mm", dpi=600)
cat("  2G saved\n")


# 2H MP score vs pseudotime

cat("── Fig 2H ──\n")

curves <- fread("pseudotime_mp_score_curves.csv")

# LOWESS + SE ribbon
# mean
mean_cols <- c("MP1_mean","MP2_mean","MP3_mean","MP4_mean")
lowess_cols <- c("MP1_lowess","MP2_lowess","MP3_lowess","MP4_lowess")
se_cols <- c("MP1_se","MP2_se","MP3_se","MP4_se")

build_long <- function(curves) {
  dfs <- list()
  for(i in 1:4) {
    mp <- paste0("MP", i)
    d <- data.table(
      pt_center = curves$pt_center,
      n_cells   = curves$n_cells,
      mean_val  = curves[[mean_cols[i]]],
      lowess    = curves[[lowess_cols[i]]],
      se        = curves[[se_cols[i]]],
      MP        = mp
    )
    dfs[[i]] <- d
  }
  rbindlist(dfs)
}

cl <- build_long(curves)
cl[, MP := factor(MP, levels=c("MP1","MP2","MP3","MP4"))]

p2h <- ggplot(cl, aes(x=pt_center, color=MP, fill=MP)) +
  geom_ribbon(aes(ymin=mean_val - se, ymax=mean_val + se), alpha=0.12, color=NA) +
  geom_line(aes(y=lowess), linewidth=0.8) +
  scale_color_manual(values=mp_colors, labels=mp_labels) +
  scale_fill_manual(values=mp_colors, labels=mp_labels) +
  labs(x="Pseudotime order", y="MP Score",
       caption="MP score gradient along the inferred pseudotime ranking. Cells are not assumed to transition between MPs.") +
  theme_pub(base_size=9) +
  theme(legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size=7),
        legend.key.size=unit(3,"mm"),
        plot.caption=element_text(size=5.5, color="grey30",
                                  hjust=0, lineheight=1.0))

ggsave("fig2h_mp_vs_pseudotime.pdf", p2h, width=110, height=65, units="mm")
ggsave("fig2h_mp_vs_pseudotime.png", p2h, width=110, height=65, units="mm", dpi=300)
cat("  2H saved\n")


# 2I Gavish (overlap coefficient

cat("── Fig 2I ──\n")

gavish <- fread("gavish_overlap.csv")
luad_mps <- gavish[[1]]  # LUAD_MP 
gavish_mat <- as.matrix(gavish[, -1, with=FALSE])
rownames(gavish_mat) <- luad_mps

# MP1-MP4
keep <- grep("^MP[1-4]$", rownames(gavish_mat))
gavish_mat <- gavish_mat[keep, , drop=FALSE]

gavish_mat <- gavish_mat[, colSums(gavish_mat > 0) > 0, drop=FALSE]

colnames(gavish_mat) <- sub("^MP[0-9]+ ", "", colnames(gavish_mat))

col_overlap <- colorRamp2(c(0, 0.15, 0.3), c("white","#FCBBA1","#E64B35"))

draw_2i <- function() {
  ht <- Heatmap(gavish_mat,
                name="Overlap\nCoefficient",
                col=col_overlap,
                cluster_rows=FALSE,
                cluster_columns=TRUE,
                clustering_distance_columns="euclidean",
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_gp=gpar(fontsize=9, fontfamily=my_font),
                column_names_gp=gpar(fontsize=5.5, fontfamily=my_font),
                column_names_rot=60,
                column_dend_height=unit(8,"mm"),
                column_dend_gp=gpar(lwd=0.3),
                rect_gp=gpar(col="grey90", lwd=0.3),
                cell_fun=function(j, i, x, y, w, h, fill) {
                  v <- gavish_mat[i, j]
                  if(v >= 0.1)
                    grid.text(sprintf("%.2f", v), x, y,
                              gp=gpar(fontsize=4, fontfamily=my_font))
                },
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=8, fontfamily=my_font),
                  labels_gp=gpar(fontsize=7, fontfamily=my_font)
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig2i_gavish_overlap.pdf", width=11, height=3)
draw_2i()
dev.off()

png("fig2i_gavish_overlap.png", width=11, height=3, units="in", res=300)
draw_2i()
dev.off()
cat("  2I saved\n")


# 2J Wilkerson × MP ( +

cat("── Fig 2J ──\n")

wilk_pct <- fread("wilkerson_MP_crosstab_pct.csv")
wilk_cnt <- fread("wilkerson_MP_crosstab_count.csv")

mp_names <- wilk_pct$dominant_MP
mat_pct <- as.matrix(wilk_pct[, .(TRU, PP, PI)])
rownames(mat_pct) <- mp_names

mat_cnt <- as.matrix(wilk_cnt[, .(TRU, PP, PI)])
rownames(mat_cnt) <- mp_names

if(max(mat_pct, na.rm=TRUE) <= 1.01) mat_pct <- mat_pct * 100

col_pct <- colorRamp2(c(0, 50, 100), c("white","#FCBBA1","#E64B35"))

draw_2j <- function() {
  ht <- Heatmap(mat_pct,
                name="% of cells",
                col=col_pct,
                cluster_rows=FALSE,
                cluster_columns=FALSE,
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_gp=gpar(fontsize=8, fontfamily=my_font),
                column_names_gp=gpar(fontsize=7, fontfamily=my_font),
                column_names_rot=45,
                rect_gp=gpar(col="grey50", lwd=0.8),
                cell_fun=function(j, i, x, y, w, h, fill) {
                  grid.text(sprintf("%.0f%%\n(%d)", mat_pct[i,j], mat_cnt[i,j]),
                            x, y, gp=gpar(fontsize=7, fontfamily=my_font))
                },
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=9, fontfamily=my_font),
                  labels_gp=gpar(fontsize=8, fontfamily=my_font)
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig2j_wilkerson_mp.pdf", width=3.5, height=3)
draw_2j()
dev.off()

png("fig2j_wilkerson_mp.png", width=3.5, height=3, units="in", res=300)
draw_2j()
dev.off()
cat("  2J saved\n")


message("=== Fig 2F-J done ===")

# Figure 2F / 2G / 2I FJ panel

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F","Unassigned"="grey80")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)


# 2F Hallmark NES =MP, =term, dendrogram

cat("── Fig 2F ──\n")

nes <- fread("hallmark_nes_heatmap.csv")
fdr <- fread("hallmark_fdr_heatmap.csv")

terms <- nes$Term
clean_term <- function(x) sub("^HALLMARK_", "", x)

nes_mat <- as.matrix(nes[, .(MP1, MP2, MP3, MP4)])
rownames(nes_mat) <- clean_term(terms)

fdr_mat <- as.matrix(fdr[, .(MP1, MP2, MP3, MP4)])
rownames(fdr_mat) <- clean_term(terms)

# MP =term
nes_t <- t(nes_mat)    # 4 × 50
fdr_t <- t(fdr_mat)

sig_t <- matrix("", nrow(fdr_t), ncol(fdr_t))
sig_t[fdr_t < 0.05]  <- "*"
sig_t[fdr_t < 0.01]  <- "**"
sig_t[fdr_t < 0.001] <- "***"

ha_left <- rowAnnotation(
  CellType = rownames(nes_t),
  col = list(CellType = mp_colors),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize=7, fontfamily=my_font),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=8, fontfamily=my_font),
    labels_gp = gpar(fontsize=7, fontfamily=my_font)
  ),
  width = unit(4, "mm")
)

col_nes <- colorRamp2(c(-3, 0, 3), c("#3C5488","white","#E64B35"))

draw_2f <- function() {
  ht <- Heatmap(nes_t,
                name="NES",
                col=col_nes,
                cluster_rows=FALSE,
                cluster_columns=TRUE,
                clustering_distance_columns="euclidean",
                column_dend_height=unit(5,"mm"),
                column_dend_gp=gpar(lwd=0.4),
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_side="left",
                row_names_gp=gpar(fontsize=8, fontfamily=my_font),
                column_names_gp=gpar(fontsize=5.5, fontfamily=my_font),
                column_names_rot=45,
                column_names_max_height=unit(35, "mm"),
                # Cells are ~4 mm wide × 4.5 mm tall clearly rectangular
                # (slightly taller than wide), and the canvas itself is
                # narrower so the strip is no longer over-elongated.
                height=unit(18, "mm"),
                left_annotation=ha_left,
                rect_gp=gpar(col="white", lwd=0.3),
                cell_fun=function(j, i, x, y, w, h, fill) {
                  if(sig_t[i, j] != "")
                    grid.text(sig_t[i, j], x, y,
                              gp=gpar(fontsize=4, fontfamily=my_font))
                },
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=8, fontfamily=my_font),
                  labels_gp=gpar(fontsize=7, fontfamily=my_font),
                  legend_height=unit(15,"mm")
                ))
  # Extra left padding so the leftmost (rotated) column label is not
  # clipped by the page edge.
  ComplexHeatmap::draw(ht, merge_legend=TRUE,
                       padding=unit(c(2, 6, 2, 2), "mm"))
}

pdf("fig2f_hallmark_nes.pdf", width=8.5, height=2.6)
draw_2f()
dev.off()

png("fig2f_hallmark_nes.png", width=8.5, height=2.6, units="in", res=300)
draw_2f()
dev.off()
cat("  2F saved\n")


# 2G Pseudotime UMAP +

cat("── Fig 2G ──\n")

pt <- fread("pseudotime_umap_winsorized.csv.gz")

p2g <- ggplot(pt[order(pt_winsorized)],
              aes(x=UMAP1, y=UMAP2, color=pt_winsorized))

if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}

p2g <- p2g +
  scale_color_viridis_c(option="inferno", name="Pseudotime",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(15,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  coord_equal() +
  labs(x="UMAP_1", y="UMAP_2") +
  theme_pub(base_size=8) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.4, color="black"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))

ggsave("fig2g_pseudotime_umap.pdf", p2g, width=90, height=80, units="mm")
ggsave("fig2g_pseudotime_umap.png", p2g, width=90, height=80, units="mm", dpi=600)
cat("  2G saved\n")


# 2I Gavish

cat("── Fig 2I ──\n")

gavish <- fread("gavish_overlap.csv")
luad_mps <- gavish[[1]]
gavish_mat <- as.matrix(gavish[, -1, with=FALSE])
rownames(gavish_mat) <- luad_mps

# MP1-MP4
keep <- grep("^MP[1-4]$", rownames(gavish_mat))
gavish_mat <- gavish_mat[keep, , drop=FALSE]

gavish_mat <- gavish_mat[, colSums(gavish_mat > 0) > 0, drop=FALSE]

cn <- colnames(gavish_mat)
cn <- sub("^MP[0-9]+ ", "", cn)         #  "MP1 " 
cn <- gsub("_", " ", cn)                 # 
cn <- gsub(" - ", "/", cn)               # 
colnames(gavish_mat) <- cn

col_overlap <- colorRamp2(c(0, 0.1, 0.25, 0.4),
                          c("#F7F7F7","#FDD49E","#EF6548","#990000"))

draw_2i <- function() {
  ht <- Heatmap(gavish_mat,
                name="Overlap\nCoefficient",
                col=col_overlap,
                cluster_rows=FALSE,
                cluster_columns=FALSE,        #  
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_side="left",
                row_names_gp=gpar(fontsize=9, fontfamily=my_font, fontface="bold"),
                column_names_gp=gpar(fontsize=5.5, fontfamily=my_font),
                column_names_rot=60,
                rect_gp=gpar(col="grey70", lwd=0.5),   # 
                width=unit(ncol(gavish_mat) * 4, "mm"),  #  4mm 
                height=unit(nrow(gavish_mat) * 6, "mm"), #  6mm 
                cell_fun=function(j, i, x, y, w, h, fill) {
                  v <- gavish_mat[i, j]
                  if(v >= 0.1)
                    grid.text(sprintf("%.2f", v), x, y,
                              gp=gpar(fontsize=4, fontfamily=my_font,
                                      col=ifelse(v > 0.25, "white", "black")))
                },
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=8, fontfamily=my_font),
                  labels_gp=gpar(fontsize=7, fontfamily=my_font),
                  legend_height=unit(20,"mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig2i_gavish_overlap.pdf", width=10, height=2.8)
draw_2i()
dev.off()

png("fig2i_gavish_overlap.png", width=10, height=2.8, units="in", res=300)
draw_2i()
dev.off()
cat("  2I saved\n")


message("=== Fig 2F / 2G / 2I  done ===")

# Figure 2A / 2G / 2I v3 +

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F","Unassigned"="grey80")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)


cat("── Fig 2A ──\n")

corr <- fread("cnmf_consensus_corr.csv")
prog_ids <- corr[[1]]
corr_mat <- as.matrix(corr[, -1, with=FALSE])
rownames(corr_mat) <- prog_ids
colnames(corr_mat) <- prog_ids

anno <- fread("cnmf_consensus_mp_annotation.csv")
anno <- anno[match(prog_ids, anno$program_id)]

# MP MP hclust
mp_order_levels <- c("MP1","MP2","MP3","MP4","MP5")
ordered_idx <- c()
for(mp in mp_order_levels) {
  idx <- which(anno$MP == mp)
  if(length(idx) == 0) next
  if(length(idx) == 1) {
    ordered_idx <- c(ordered_idx, idx)
  } else {
    sub_mat <- corr_mat[idx, idx, drop=FALSE]
    sub_hc <- hclust(as.dist(1 - sub_mat), method="average")
    ordered_idx <- c(ordered_idx, idx[sub_hc$order])
  }
}
corr_mat <- corr_mat[ordered_idx, ordered_idx]
anno <- anno[ordered_idx]

# MP block
mp_vec <- anno$MP
mp_blocks <- data.frame(MP=character(), start=integer(), end=integer(),
                        stringsAsFactors=FALSE)
for(mp in mp_order_levels) {
  idx <- which(mp_vec == mp)
  if(length(idx) == 0) next
  mp_blocks <- rbind(mp_blocks, data.frame(MP=mp, start=min(idx), end=max(idx)))
}

ha_left <- rowAnnotation(
  MetaProgram = anno$MP,
  col = list(MetaProgram = mp_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,   # hide name to avoid overlap with cNMF labels
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontfamily=my_font, fontface="bold"),
    labels_gp = gpar(fontsize=8, fontfamily=my_font)
  ),
  width = unit(4, "mm")
)

ha_top <- HeatmapAnnotation(
  MetaProgram = anno$MP,
  col = list(MetaProgram = mp_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  height = unit(4, "mm")
)

cat("  corr range:", range(corr_mat), "\n")
diag(corr_mat) <- NA

# off-diagonal
off_diag <- corr_mat[!is.na(corr_mat)]
cat("  off-diagonal range:", round(range(off_diag), 3), "\n")
cat("  off-diagonal quantiles:", round(quantile(off_diag, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)), 3), "\n")

val_min <- min(off_diag)
val_mid <- median(off_diag)
val_max <- max(off_diag)

col_fun <- colorRamp2(
  c(val_min, val_mid, val_max),
  c("#4575B4", "#F7F7F7", "#D73027")  # →→ (RdYlBu reversed)
)

ht_name <- "Spearman\nCorrelation"

draw_2a <- function() {
  n <- ncol(corr_mat)
  cellsize <- 7  # mm per cell -> square panel
  ht <- Heatmap(corr_mat,
                name = ht_name,
                col = col_fun,
                na_col = "#7F0000",  #  =1 
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize=8, fontfamily=my_font),
                column_names_gp = gpar(fontsize=8, fontfamily=my_font),
                column_names_rot = 45,
                left_annotation = ha_left,
                top_annotation = ha_top,
                rect_gp = gpar(col="white", lwd=0.5),
                width  = unit(n * cellsize, "mm"),
                height = unit(n * cellsize, "mm"),
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=9, fontfamily=my_font, fontface="bold"),
                  labels_gp = gpar(fontsize=8, fontfamily=my_font),
                  legend_height = unit(25, "mm")
                ))
  ComplexHeatmap::draw(ht, merge_legend=TRUE)

  # Bold MP labels with white halo + colored fill rect on diagonal blocks
  # so the 5 MP markers stand out clearly against the heatmap.
  for(i in seq_len(nrow(mp_blocks))) {
    s <- mp_blocks$start[i]; e <- mp_blocks$end[i]; mp <- mp_blocks$MP[i]
    decorate_heatmap_body(ht_name, {
      x1 <- (s-1)/n; x2 <- e/n
      y1 <- 1-e/n;   y2 <- 1-(s-1)/n
      cx <- (x1+x2)/2; cy <- (y1+y2)/2
      # Bordered diagonal block (thicker)
      grid.rect(x=unit(x1,"npc"), y=unit(y1,"npc"),
                width=unit(x2-x1,"npc"), height=unit(y2-y1,"npc"),
                just=c("left","bottom"),
                gp=gpar(col=mp_colors[mp], fill=NA, lwd=3))
      # Colored pill behind label so text doesn't fight the heatmap fill
      lab_w <- unit(7, "mm"); lab_h <- unit(4.5, "mm")
      grid.roundrect(x=unit(cx,"npc"), y=unit(cy,"npc"),
                     width=lab_w, height=lab_h,
                     r = unit(1, "mm"),
                     gp=gpar(fill=mp_colors[mp], col="white", lwd=0.8))
      grid.text(mp, x=unit(cx,"npc"), y=unit(cy,"npc"),
                gp=gpar(col="white", fontsize=10,
                        fontfamily=my_font, fontface="bold"))
    })
  }
}

# Square heatmap panel: 15 × 15 cells × 7mm = 105×105mm + margins for
# annotations / labels / legend. Use 7×7 in canvas to fit comfortably.
pdf("fig2a_gep_heatmap.pdf", width=7, height=7)
draw_2a()
dev.off()

png("fig2a_gep_heatmap.png", width=7, height=7, units="in", res=300)
draw_2a()
dev.off()
cat("  2A saved\n")


# 2G Trajectory coord_equal

cat("── Fig 2G ──\n")

traj_dt <- fread("monocle_trajectory.csv.gz")
traj <- as.data.frame(traj_dt)
traj$pseudotime <- as.numeric(traj$pseudotime)

coord_cols <- names(traj)
x_col <- coord_cols[grep("^(DC1|Component_1|diffmap_1|DC_1)", coord_cols)][1]
y_col <- coord_cols[grep("^(DC2|Component_2|diffmap_2|DC_2)", coord_cols)][1]
cat("  using x=", x_col, " y=", y_col, "\n")

traj$x <- as.numeric(traj[[x_col]])
traj$y <- as.numeric(traj[[y_col]])

cat("  x range:", range(traj$x, na.rm=TRUE), "\n")
cat("  y range:", range(traj$y, na.rm=TRUE), "\n")

# 1-99 percentile
clip <- function(v, lo=0.01, hi=0.99) {
  q <- quantile(v, c(lo, hi), na.rm=TRUE)
  pmax(pmin(v, q[2]), q[1])
}
traj$x <- clip(traj$x)
traj$y <- clip(traj$y)

# pseudotime
traj <- traj[order(traj$pseudotime), ]

# ── 1: pseudotime ──
p2g <- ggplot(traj, aes(x=x, y=y, color=pseudotime))
if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g <- p2g +
  scale_color_viridis_c(option="inferno", name="Relative\norder",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(18,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  labs(x="Component 1", y="Component 2") +
  theme_pub(base_size=9) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.4, color="black"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))

ggsave("fig2g_pseudotime_trajectory.pdf", p2g, width=80, height=70, units="mm")
ggsave("fig2g_pseudotime_trajectory.png", p2g, width=80, height=70, units="mm", dpi=600)
cat("  2G pseudotime saved\n")

traj_mp <- traj[traj$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
traj_mp$dominant_MP <- factor(traj_mp$dominant_MP,
                              levels=c("MP1","MP2","MP3","MP4"))
set.seed(42)
traj_mp <- traj_mp[sample(nrow(traj_mp)), ]

p2g_mp <- ggplot(traj_mp, aes(x=x, y=y, color=dominant_MP))
if(has_ggrastr) {
  p2g_mp <- p2g_mp + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g_mp <- p2g_mp + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g_mp <- p2g_mp +
  scale_color_manual(values=mp_colors, labels=mp_labels, name="celltype") +
  guides(color=guide_legend(override.aes=list(size=2.5, alpha=1))) +
  labs(x="Component 1", y="Component 2") +
  theme_pub(base_size=9) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.4, color="black"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))

ggsave("fig2g_trajectory_by_mp.pdf", p2g_mp, width=95, height=70, units="mm")
ggsave("fig2g_trajectory_by_mp.png", p2g_mp, width=95, height=70, units="mm", dpi=600)
cat("  2G MP saved\n")


# 2I Gavish cosine binary cosine ≥0 →

cat("── Fig 2I ──\n")

cos_df <- fread("gavish_cosine_similarity_cnmf.csv")
prog_ids_i <- cos_df[[1]]
cos_mat <- as.matrix(cos_df[, -1, with=FALSE])
rownames(cos_mat) <- prog_ids_i

cn <- colnames(cos_mat)
cn <- sub("^MP[0-9]+ ", "", cn)
cn <- gsub("_", " ", cn)
colnames(cos_mat) <- cn

anno_i <- fread("cnmf_consensus_mp_annotation.csv")
anno_i <- anno_i[match(prog_ids_i, anno_i$program_id)]
row_order <- order(factor(anno_i$MP, levels=c("MP1","MP2","MP3","MP4","MP5")),
                   prog_ids_i)
cos_mat <- cos_mat[row_order, ]
anno_i <- anno_i[row_order]

# Drop reference programs that don't meaningfully match ANY cnmf program
# their columns are flat-white and just dilute the figure.
# Cosine values are small overall (max ≈ 0.24); 0.10 corresponds to
# roughly the 95th percentile of all entries, i.e. a real signal.
ref_thresh <- 0.10
col_max  <- apply(cos_mat, 2, max, na.rm = TRUE)
keep_col <- col_max >= ref_thresh
cat(sprintf("  reference programs kept: %d / %d (max cosine ≥ %.2f)\n",
            sum(keep_col), length(col_max), ref_thresh))
# Order kept columns by their best match strength (strongest hit first)
ord_col  <- order(-col_max[keep_col])
cos_mat  <- cos_mat[, which(keep_col)[ord_col]]

cat("  cosine range (filtered):", range(cos_mat), "\n")

# Sequential palette: cosine ∈ [0, 1] is non-negative diverging colours
# would flood the many near-zero cells with strong colour. Map zero → white
# so only meaningful (high cosine) matches stand out. Quantile-aware
# breakpoints keep the dynamic range readable even when most values are low.
val_max <- max(cos_mat, na.rm = TRUE)
qs <- quantile(cos_mat[cos_mat > 0], probs = c(0.5, 0.85, 0.98), na.rm = TRUE)
col_cos <- colorRamp2(
  c(0, as.numeric(qs[1]), as.numeric(qs[2]), as.numeric(qs[3]), val_max),
  c("white", "#FEE0D2", "#FC9272", "#E64B35", "#67000D")
)

draw_2i <- function() {
  ht <- Heatmap(cos_mat,
                name="Cosine\nsimilarity",
                col=col_cos,
                cluster_rows=FALSE,
                cluster_columns=FALSE,
                show_row_names=TRUE,
                show_column_names=TRUE,
                row_names_side="left",
                row_names_gp=gpar(fontsize=7, fontfamily=my_font),
                column_names_gp=gpar(fontsize=5, fontfamily=my_font),
                column_names_rot=60,
                row_title="Query Programs",
                row_title_gp=gpar(fontsize=9, fontfamily=my_font, fontface="bold"),
                column_title="Reference Programs",
                column_title_side="bottom",
                column_title_gp=gpar(fontsize=9, fontfamily=my_font, fontface="bold"),
                rect_gp=gpar(col="grey80", lwd=0.3),
                width=unit(ncol(cos_mat) * 4.5, "mm"),
                height=unit(nrow(cos_mat) * 4.5, "mm"),
                heatmap_legend_param=list(
                  title_gp=gpar(fontsize=8, fontfamily=my_font),
                  labels_gp=gpar(fontsize=7, fontfamily=my_font),
                  legend_height=unit(20,"mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig2i_gavish_cosine.pdf", width=12, height=4.5)
draw_2i()
dev.off()

png("fig2i_gavish_cosine.png", width=12, height=4.5, units="in", res=300)
draw_2i()
dev.off()
cat("  2I saved\n")


message("=== Fig 2A / 2G / 2I v3 done ===")

# Figure 2G UMAP + pseudotime_winsorized

library(data.table)
library(ggplot2)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)

traj_dt <- fread("monocle_trajectory.csv.gz")
traj <- as.data.frame(traj_dt)
# Column is "pseudotime" in monocle_trajectory.csv.gz (no _winsorized suffix
# in this file; earlier patches to other blocks used a different file).
pt_col <- if ("pseudotime_winsorized" %in% names(traj))
  "pseudotime_winsorized" else "pseudotime"
traj$pt <- as.numeric(traj[[pt_col]])
# Cap at the 99th percentile so a few outliers don't compress the colour range
.cap <- quantile(traj$pt, 0.99, na.rm = TRUE)
traj$pt <- pmin(traj$pt, .cap)
traj$x <- as.numeric(traj$Component_1)
traj$y <- as.numeric(traj$Component_2)

# ── 1: pseudotime ──
traj_pt <- traj[order(traj$pt), ]

p2g <- ggplot(traj_pt, aes(x=x, y=y, color=pt))
if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g <- p2g +
  scale_color_viridis_c(option="inferno", name="Relative\norder",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(18,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  coord_equal() +
  labs(x="Component 1", y="Component 2") +
  theme_pub(base_size=9) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.4, color="black"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))

ggsave("fig2g_pseudotime_trajectory.pdf", p2g, width=80, height=70, units="mm")
ggsave("fig2g_pseudotime_trajectory.png", p2g, width=80, height=70, units="mm", dpi=600)
cat("  2G pseudotime saved\n")

traj_mp <- traj[traj$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
traj_mp$dominant_MP <- factor(traj_mp$dominant_MP, levels=c("MP1","MP2","MP3","MP4"))
set.seed(42)
traj_mp <- traj_mp[sample(nrow(traj_mp)), ]

p2g_mp <- ggplot(traj_mp, aes(x=x, y=y, color=dominant_MP))
if(has_ggrastr) {
  p2g_mp <- p2g_mp + ggrastr::rasterise(geom_point(size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g_mp <- p2g_mp + geom_point(size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g_mp <- p2g_mp +
  scale_color_manual(values=mp_colors, labels=mp_labels, name="celltype") +
  guides(color=guide_legend(override.aes=list(size=2.5, alpha=1))) +
  coord_equal() +
  labs(x="Component 1", y="Component 2") +
  theme_pub(base_size=9) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.4, color="black"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=6))

ggsave("fig2g_trajectory_by_mp.pdf", p2g_mp, width=95, height=70, units="mm")
ggsave("fig2g_trajectory_by_mp.png", p2g_mp, width=95, height=70, units="mm", dpi=600)
cat("  2G MP saved\n")

message("=== Fig 2G done ===")

# Figure 2G FA layout +

library(data.table)
library(ggplot2)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)

cells <- as.data.frame(fread("trajectory_cells.csv.gz"))
skel  <- as.data.frame(fread("trajectory_graph.csv"))

cells$pt <- as.numeric(cells$pseudotime_winsorized)
cells$x  <- as.numeric(cells$FA1)
cells$y  <- as.numeric(cells$FA2)
skel$x   <- as.numeric(skel$x)
skel$y   <- as.numeric(skel$y)
skel$pt  <- as.numeric(skel$pseudotime)
skel <- skel[order(skel$pt), ]

# NOTE (Option A framing): we are now treating the Monocle3 ordering as a
# *pseudotime continuum*, not a directional trajectory. Direction arrows
# and "Less / More differentiated" anchor labels were removed because
#   1) snapshot scRNA-seq cannot establish causal directionality between MPs;
#   2) MPs are co-existing programs, not discrete cell states.
# Caption notes this explicitly.

# ── 1: pseudotime ──
cells_pt <- cells[order(cells$pt), ]

p2g <- ggplot(cells_pt, aes(x=x, y=y))
if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(
    geom_point(aes(color=pt), size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(aes(color=pt), size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g <- p2g +
  scale_color_viridis_c(option="inferno", name="Pseudotime\norder",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(18,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  labs(x=NULL, y=NULL,
       caption="Pseudotime order = inferred ranking, not a temporal sequence of state transitions.") +
  theme_pub(base_size=9) +
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.caption=element_text(size=5.5, color="grey30",
                                  hjust=0, lineheight=1.0)) +
  umap_arrow_axes(cells_pt, "x", "y",
                  label_x = "Component 1", label_y = "Component 2")

ggsave("fig2g_pseudotime_continuum.pdf", p2g, width=85, height=75, units="mm")
ggsave("fig2g_pseudotime_continuum.png", p2g, width=85, height=75, units="mm", dpi=600)
# remove the trajectory-named outputs from any earlier run
for (ext in c("pdf","png"))
  if (file.exists(paste0("fig2g_pseudotime_trajectory.", ext)))
    file.remove(paste0("fig2g_pseudotime_trajectory.", ext))
cat("  2G pseudotime continuum saved\n")


cells_mp <- cells[cells$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
cells_mp$dominant_MP <- factor(cells_mp$dominant_MP, levels=c("MP1","MP2","MP3","MP4"))
set.seed(42)
cells_mp <- cells_mp[sample(nrow(cells_mp)), ]

p2g_mp <- ggplot(cells_mp, aes(x=x, y=y))
if(has_ggrastr) {
  p2g_mp <- p2g_mp + ggrastr::rasterise(
    geom_point(aes(color=dominant_MP), size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g_mp <- p2g_mp + geom_point(aes(color=dominant_MP), size=0.18, stroke=0, alpha=0.45, shape=16)
}
p2g_mp <- p2g_mp +
  scale_color_manual(values=mp_colors, labels=mp_labels, name="Dominant MP") +
  guides(color=guide_legend(override.aes=list(size=2.5, alpha=1))) +
  labs(x=NULL, y=NULL,
       caption="MP-coloured pseudotime continuum. Cells are not assumed to transition between MPs.") +
  theme_pub(base_size=9) +
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.caption=element_text(size=5.5, color="grey30",
                                  hjust=0, lineheight=1.0)) +
  umap_arrow_axes(cells_mp, "x", "y",
                  label_x = "Component 1", label_y = "Component 2")

ggsave("fig2g_continuum_by_mp.pdf", p2g_mp, width=100, height=75, units="mm")
ggsave("fig2g_continuum_by_mp.png", p2g_mp, width=100, height=75, units="mm", dpi=600)
# remove the trajectory-named outputs from any earlier run
for (ext in c("pdf","png"))
  if (file.exists(paste0("fig2g_trajectory_by_mp.", ext)))
    file.remove(paste0("fig2g_trajectory_by_mp.", ext))
cat("  2G MP continuum saved\n")

message("=== Fig 2G done ===")


# Figure 2G/H Monocle3 principal graph on UMAP

library(data.table)
library(ggplot2)
library(R.utils)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig2" else "${WORK_ROOT}/luad_figures/fig2")

# (override removed: keep Arial from top)
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")
mp_labels <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
               "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

has_ggrastr <- requireNamespace("ggrastr", quietly=TRUE)

cells <- as.data.frame(fread("monocle3_cells.csv.gz"))
seg   <- as.data.frame(fread("monocle3_graph_segments.csv"))

cells$pt <- as.numeric(cells$monocle3_pseudotime)

cat("  UMAP1 range:", range(cells$UMAP1, na.rm=TRUE), "\n")
cat("  UMAP2 range:", range(cells$UMAP2, na.rm=TRUE), "\n")
cat("  pseudotime range:", range(cells$pt, na.rm=TRUE), "\n")
cat("  segments:", nrow(seg), "\n")


# ── 1: pseudotime + ──

cells_pt <- cells[order(cells$pt), ]

p2g <- ggplot(cells_pt, aes(x=UMAP1, y=UMAP2))

if(has_ggrastr) {
  p2g <- p2g + ggrastr::rasterise(
    geom_point(aes(color=pt), size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2g <- p2g + geom_point(aes(color=pt), size=0.18, stroke=0, alpha=0.45, shape=16)
}

p2g <- p2g +
  # Monocle3 principal graph skeleton (kept as a non-directional ordering scaffold)
  geom_segment(data=seg,
               aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
               inherit.aes=FALSE,
               linewidth=0.5, color="black", alpha=0.8) +
  scale_color_viridis_c(option="inferno", name="Pseudotime\norder",
                        guide=guide_colorbar(barwidth=unit(3,"mm"),
                                             barheight=unit(18,"mm"),
                                             frame.colour="black",
                                             frame.linewidth=0.3)) +
  coord_equal() +
  labs(x=NULL, y=NULL,
       caption="Pseudotime order = inferred ranking, not a temporal sequence of state transitions.") +
  theme_pub(base_size=9) +
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.caption=element_text(size=5.5, color="grey30",
                                  hjust=0, lineheight=1.0)) +
  umap_arrow_axes(cells_pt, "UMAP1", "UMAP2")

ggsave("fig2g_pseudotime_continuum.pdf", p2g, width=90, height=80, units="mm")
ggsave("fig2g_pseudotime_continuum.png", p2g, width=90, height=80, units="mm", dpi=600)
cat("  2G pseudotime continuum (monocle3) saved\n")


cells_mp <- cells[cells$dominant_MP %in% c("MP1","MP2","MP3","MP4"), ]
cells_mp$dominant_MP <- factor(cells_mp$dominant_MP, levels=c("MP1","MP2","MP3","MP4"))
set.seed(42)
cells_mp <- cells_mp[sample(nrow(cells_mp)), ]

p2h <- ggplot(cells_mp, aes(x=UMAP1, y=UMAP2))

if(has_ggrastr) {
  p2h <- p2h + ggrastr::rasterise(
    geom_point(aes(color=dominant_MP), size=0.18, stroke=0, alpha=0.45, shape=16), dpi=600)
} else {
  p2h <- p2h + geom_point(aes(color=dominant_MP), size=0.18, stroke=0, alpha=0.45, shape=16)
}

p2h <- p2h +
  geom_segment(data=seg,
               aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
               inherit.aes=FALSE,
               linewidth=0.5, color="black", alpha=0.8) +
  scale_color_manual(values=mp_colors, labels=mp_labels, name="Dominant MP") +
  guides(color=guide_legend(override.aes=list(size=2.5, alpha=1))) +
  coord_equal() +
  labs(x=NULL, y=NULL,
       caption="MP-coloured pseudotime continuum (Monocle3). Cells are not assumed to transition between MPs.") +
  theme_pub(base_size=9) +
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        plot.caption=element_text(size=5.5, color="grey30",
                                  hjust=0, lineheight=1.0)) +
  umap_arrow_axes(cells_mp, "UMAP1", "UMAP2")

ggsave("fig2h_continuum_by_mp.pdf", p2h, width=100, height=80, units="mm")
ggsave("fig2h_continuum_by_mp.png", p2h, width=100, height=80, units="mm", dpi=600)
cat("  2H MP continuum (monocle3) saved\n")

# Cleanup obsolete trajectory-named outputs from this and earlier blocks
for (b in c("fig2g_pseudotime_trajectory","fig2g_trajectory_by_mp",
            "fig2h_trajectory_by_mp","fig2g_trajectory_v3",
            "fig2h_pseudotime_curves_v3"))
  for (ext in c("pdf","png"))
    if (file.exists(paste0(b, ".", ext))) file.remove(paste0(b, ".", ext))

message("=== Fig 2G/H continuum done ===")
