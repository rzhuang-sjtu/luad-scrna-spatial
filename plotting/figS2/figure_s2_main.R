## Supplementary Figure 2  |  CNA Profiling (3 panels)

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(R.utils)
suppressPackageStartupMessages({
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig_s2_cna" else
  "${WORK_ROOT}/luad_figures/fig_s2_cna")

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


# S2A Average CNA by sample (, =bin, =sample

cat("── S2A ──\n")

sample_cna <- fread("cna_by_sample.csv.gz")
meta <- fread("cna_sample_metadata.csv")

coord_cols <- c("chr","start","end","abspos")
coord_cols <- intersect(coord_cols, names(sample_cna))
cat("  coord columns:", coord_cols, "\n")

chr_vec <- sample_cna$chr
mat <- as.matrix(sample_cna[, !coord_cols, with=FALSE])
rownames(mat) <- NULL

cat("  matrix:", nrow(mat), "bins x", ncol(mat), "samples\n")
cat("  CNA range:", round(range(mat, na.rm=TRUE), 3), "\n")

chr_vec <- factor(chr_vec, levels=c(paste0("chr", 1:22), "chrX","chrY"))
if(any(grepl("^[0-9]+$", as.character(sample_cna$chr)))) {
  chr_vec <- factor(sample_cna$chr, levels=c(1:22, "X","Y"))
}

# MP majority
meta <- meta[match(colnames(mat), meta$sample_id)]
if(any(is.na(meta$sample_id))) {
  # patient_id
  meta <- fread("cna_sample_metadata.csv")
  meta <- meta[match(colnames(mat), meta$patient_id)]
}

# annotation: dominant MP
if("dominant_MP_majority" %in% names(meta)) {
  sample_mp <- meta$dominant_MP_majority
} else {
  sample_mp <- rep("Unknown", ncol(mat))
}

col_order <- order(factor(sample_mp, levels=c("MP1","MP2","MP3","MP4")))
mat <- mat[, col_order]
sample_mp <- sample_mp[col_order]

# annotation
ha_top <- HeatmapAnnotation(
  MP = sample_mp,
  col = list(MP = mp_colors),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize=7, fontfamily=my_font),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=8, fontfamily=my_font),
    labels_gp = gpar(fontsize=7, fontfamily=my_font)
  ),
  height = unit(3, "mm")
)

col_cna <- colorRamp2(c(-0.3, -0.1, 0, 0.1, 0.3),
                      c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B"))

draw_s2a <- function() {
  ht <- Heatmap(mat,
                name = "Avg. CNA",
                col = col_cna,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize=3, fontfamily=my_font),
                column_names_rot = 45,
                top_annotation = ha_top,
                row_split = chr_vec,
                row_gap = unit(0.3, "mm"),
                row_title_gp = gpar(fontsize=5, fontfamily=my_font),
                row_title_rot = 0,
                use_raster = TRUE,
                raster_quality = 3,
                # Force a square plotting body physical heatmap area is
                # 160 × 160 mm regardless of how many samples / bins.
                width  = unit(160, "mm"),
                height = unit(160, "mm"),
                column_title = "Average CNA by Sample",
                column_title_gp = gpar(fontsize=10, fontfamily=my_font, fontface="bold"),
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=8, fontfamily=my_font),
                  labels_gp = gpar(fontsize=7, fontfamily=my_font),
                  legend_height = unit(20, "mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig_s2a_cna_by_sample.pdf", width=8, height=8)
draw_s2a()
dev.off()

png("fig_s2a_cna_by_sample.png", width=8, height=8, units="in", res=300)
draw_s2a()
dev.off()
cat("  S2A saved\n")


## S2B  Patient CNA hierarchical clustering dendrogram

cat("── S2B ──\n")

dist_df <- fread("cna_patient_dist.csv")
dist_ids <- dist_df[[1]]
dist_mat <- as.matrix(dist_df[, -1, with=FALSE])
rownames(dist_mat) <- dist_ids
colnames(dist_mat) <- dist_ids

d <- as.dist(dist_mat)
hc <- hclust(d, method="ward.D2")

if(requireNamespace("dendextend", quietly=TRUE)) {
  library(dendextend)
  dend <- as.dendrogram(hc)
  
  leaf_ids <- labels(dend)
  meta_full <- fread("cna_sample_metadata.csv")
  
  # sample_id patient_id
  if(all(leaf_ids %in% meta_full$sample_id)) {
    leaf_mp <- meta_full$dominant_MP_majority[match(leaf_ids, meta_full$sample_id)]
  } else if(all(leaf_ids %in% meta_full$patient_id)) {
    leaf_mp <- meta_full$dominant_MP_majority[match(leaf_ids, meta_full$patient_id)]
  } else {
    leaf_mp <- rep("Unknown", length(leaf_ids))
  }
  
  leaf_col <- mp_colors[leaf_mp]
  leaf_col[is.na(leaf_col)] <- "grey50"
  labels_colors(dend) <- leaf_col
  labels_cex(dend) <- 0.5
  
  pdf("fig_s2b_cna_dendrogram.pdf", width=10, height=4)
  par(mar=c(5, 4, 2, 1), family=my_font)
  plot(dend,
       main="Hierarchical clustering of patients based on CNA profiles",
       ylab="Height", cex.main=0.9)
  legend("topright", legend=names(mp_colors), fill=mp_colors,
         border=NA, bty="n", cex=0.7)
  dev.off()
  
  png("fig_s2b_cna_dendrogram.png", width=10, height=4, units="in", res=300)
  par(mar=c(5, 4, 2, 1), family=my_font)
  plot(dend,
       main="Hierarchical clustering of patients based on CNA profiles",
       ylab="Height", cex.main=0.9)
  legend("topright", legend=names(mp_colors), fill=mp_colors,
         border=NA, bty="n", cex=0.7)
  dev.off()
} else {
  pdf("fig_s2b_cna_dendrogram.pdf", width=10, height=4)
  par(mar=c(5, 4, 2, 1), family=my_font)
  plot(hc, main="Hierarchical clustering of patients based on CNA profiles",
       xlab="", ylab="Height", cex=0.5)
  dev.off()
  png("fig_s2b_cna_dendrogram.png", width=10, height=4, units="in", res=300)
  par(mar=c(5, 4, 2, 1), family=my_font)
  plot(hc, main="Hierarchical clustering of patients based on CNA profiles",
       xlab="", ylab="Height", cex=0.5)
  dev.off()
}
cat("  S2B saved\n")


# S2C Average CNA by MP (, =bin, =MP1-MP4

cat("── S2C ──\n")

mp_cna <- fread("cna_by_mp.csv.gz")

coord_cols_mp <- intersect(c("chr","start","end","abspos"), names(mp_cna))
chr_vec_mp <- mp_cna$chr

# factor
if(any(grepl("^[0-9]+$", as.character(chr_vec_mp)))) {
  chr_vec_mp <- factor(chr_vec_mp, levels=c(1:22, "X","Y"))
} else {
  chr_vec_mp <- factor(chr_vec_mp, levels=c(paste0("chr", 1:22), "chrX","chrY"))
}

mp_mat <- as.matrix(mp_cna[, .(MP1, MP2, MP3, MP4)])

cat("  MP CNA range:", round(range(mp_mat, na.rm=TRUE), 3), "\n")

# MP annotation (=MP
ha_top_mp <- HeatmapAnnotation(
  MP = c("MP1","MP2","MP3","MP4"),
  col = list(MP = mp_colors),
  show_legend = FALSE,
  height = unit(3, "mm")
)

draw_s2c <- function() {
  ht <- Heatmap(mp_mat,
                name = "Avg. CNA",
                col = col_cna,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize=10, fontfamily=my_font, fontface="bold"),
                column_names_rot = 45,
                top_annotation = ha_top_mp,
                row_split = chr_vec_mp,
                row_gap = unit(0.3, "mm"),
                row_title_gp = gpar(fontsize=5, fontfamily=my_font),
                row_title_rot = 0,
                use_raster = TRUE,
                raster_quality = 3,
                # Force a square plotting body 130 × 130 mm regardless of
                # the small column count
                width  = unit(130, "mm"),
                height = unit(130, "mm"),
                column_title = "Average CNA by MP",
                column_title_gp = gpar(fontsize=10, fontfamily=my_font, fontface="bold"),
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=8, fontfamily=my_font),
                  labels_gp = gpar(fontsize=7, fontfamily=my_font),
                  legend_height = unit(20, "mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig_s2c_cna_by_mp.pdf", width=6, height=6)
draw_s2c()
dev.off()

png("fig_s2c_cna_by_mp.png", width=6, height=6, units="in", res=300)
draw_s2c()
dev.off()
cat("  S2C saved\n")


message("=== Supp Fig 2 (CNA) done ===")
