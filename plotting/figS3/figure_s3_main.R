# Supplementary Figure 3 | cNMF QC ( HCC Supp Fig 3
## 7 panels: A-G

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)
library(R.utils)
suppressPackageStartupMessages({
  if (requireNamespace("showtext", quietly = TRUE)) library(showtext)
  if (requireNamespace("sysfonts", quietly = TRUE)) library(sysfonts)
})

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig_s2" else
  "${WORK_ROOT}/luad_figures/fig_s2")

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

mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488",
               "MP5"="#F39B7F")

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

shannon_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log2(p))
}

# ── 770 GEP → 15 consensus cluster ──
gep_anno <- fread("../fig2/gep_mp_annotation.csv")
gep_anno <- gep_anno[MP != "MP5"]

corr_770 <- fread("gep_spearman_corr.csv")
gep_ids_770 <- corr_770[[1]]
corr_mat_770 <- as.matrix(corr_770[, -1, with=FALSE])
rownames(corr_mat_770) <- gep_ids_770

keep <- intersect(gep_anno$gep_id, gep_ids_770)
corr_sub <- corr_mat_770[keep, keep]
hc_770 <- hclust(as.dist(1 - corr_sub), method="average")
clust_15 <- cutree(hc_770, k=15)

gep_anno <- gep_anno[gep_id %in% keep]
gep_anno[, cluster := paste0("cNMF_", clust_15[match(gep_id, names(clust_15))])]


# S3A K selection curves (per-patient

cat("── S3A ──\n")

kcurves <- fread("s2_k_error_curves.csv")
k_sel   <- fread("s2_k_selection_summary.csv")

p_s3a <- ggplot(kcurves, aes(x=K)) +
  geom_line(aes(y=stability, color="Stability"), linewidth=0.4) +
  geom_point(aes(y=stability, color="Stability"), size=0.6) +
  geom_line(aes(y=prediction_error / max(prediction_error, na.rm=TRUE),
                color="Prediction Error"), linewidth=0.4) +
  geom_point(aes(y=prediction_error / max(prediction_error, na.rm=TRUE),
                 color="Prediction Error"), size=0.6) +
  geom_vline(data=k_sel, aes(xintercept=selected_k),
             linetype="dashed", color="grey40", linewidth=0.2) +
  facet_wrap(~patient_key, scales="free_y", ncol=8) +
  scale_color_manual(values=c("Stability"="#3C5488","Prediction Error"="#E64B35")) +
  labs(x="Number of Components (K)", y="Stability / Scaled Error") +
  theme_pub(base_size=6) +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(3,"mm"),
        strip.text=element_text(size=4.5),
        axis.text=element_text(size=4),
        panel.spacing=unit(0.5,"mm"))

ggsave("fig_s3a_k_curves.pdf", p_s3a, width=350, height=220, units="mm")
ggsave("fig_s3a_k_curves.png", p_s3a, width=350, height=220, units="mm", dpi=300)
cat("  S3A saved\n")


# S3B Distance matrix of 15 consensus programs

cat("── S3B ──\n")

corr_15 <- fread("../fig2/cnmf_consensus_corr.csv")
prog_ids <- corr_15[[1]]
corr_mat_15 <- as.matrix(corr_15[, -1, with=FALSE])
rownames(corr_mat_15) <- prog_ids
colnames(corr_mat_15) <- prog_ids

dist_15 <- as.matrix(dist(t(corr_mat_15)))  # program 
# 1-corr dissimilarity
dissim_15 <- 1 - corr_mat_15

hc_15 <- hclust(as.dist(dissim_15), method="average")

col_dist <- colorRamp2(c(0, 0.5, 1, 1.5),
                       c("#FDE725","#5DC863","#21908C","#440154"))

draw_s3b <- function() {
  ht <- Heatmap(dissim_15,
                name = "Distance\n(1 - Spearman)",
                col = col_dist,
                cluster_rows = hc_15,
                cluster_columns = hc_15,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize=8, fontfamily=my_font),
                column_names_gp = gpar(fontsize=8, fontfamily=my_font),
                column_names_rot = 45,
                row_dend_width = unit(12, "mm"),
                column_dend_height = unit(12, "mm"),
                row_dend_gp = gpar(lwd=0.5),
                column_dend_gp = gpar(lwd=0.5),
                rect_gp = gpar(col="white", lwd=0.3),
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize=8, fontfamily=my_font),
                  labels_gp = gpar(fontsize=7, fontfamily=my_font),
                  legend_height = unit(20, "mm")
                ))
  ComplexHeatmap::draw(ht)
}

pdf("fig_s3b_distance_matrix.pdf", width=6, height=5.5)
draw_s3b()
dev.off()

png("fig_s3b_distance_matrix.png", width=6, height=5.5, units="in", res=300)
draw_s3b()
dev.off()
cat("  S3B saved\n")


## S3C  Proportion of 15 consensus clusters across samples

cat("── S3C ──\n")

# sample (patient) cluster GEP
# patient_key sample ID
clust_by_pt <- gep_anno[, .N, by=.(patient_key, cluster)]
clust_by_pt[, total := sum(N), by=patient_key]
clust_by_pt[, pct := N / total]

n_clust <- length(unique(clust_by_pt$cluster))
clust_levels <- paste0("cNMF_", 1:n_clust)
clust_levels <- intersect(clust_levels, unique(clust_by_pt$cluster))
clust_levels <- sort(unique(clust_by_pt$cluster))

# Paired + Set3
clust_pal <- c(RColorBrewer::brewer.pal(12, "Paired"),
               RColorBrewer::brewer.pal(8, "Set3"))[1:length(clust_levels)]
names(clust_pal) <- clust_levels

clust_by_pt[, cluster := factor(cluster, levels=rev(clust_levels))]

# sample : dominant cluster
dom_clust <- clust_by_pt[, .SD[which.max(pct)], by=patient_key]
pt_order <- dom_clust[order(cluster, -pct)]$patient_key
clust_by_pt[, patient_key := factor(patient_key, levels=pt_order)]

p_s3c <- ggplot(clust_by_pt, aes(x=patient_key, y=pct, fill=cluster)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values=clust_pal, name="cNMF_cluster") +
  scale_y_continuous(labels=percent, expand=c(0,0)) +
  labs(x=NULL, y="Proportion") +
  theme_pub(base_size=6) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=3.5),
        legend.position="right",
        legend.key.size=unit(3,"mm"),
        legend.text=element_text(size=5))

ggsave("fig_s3c_cluster_proportion.pdf", p_s3c, width=250, height=65, units="mm")
ggsave("fig_s3c_cluster_proportion.png", p_s3c, width=250, height=65, units="mm", dpi=300)
cat("  S3C saved\n")


## S3D  Patient mixing entropy per consensus cluster (15 bars)

cat("── S3D ──\n")

cluster_ent <- gep_anno[, {
  pt_counts <- table(patient_key)
  .(entropy=shannon_entropy(pt_counts),
    n_patients=length(unique(patient_key)),
    n_geps=.N)
}, by=cluster]

cluster_ent <- cluster_ent[order(entropy, decreasing=TRUE)]
cluster_ent[, cluster := factor(cluster, levels=cluster)]

mean_ent <- mean(cluster_ent$entropy)

p_s3d <- ggplot(cluster_ent, aes(x=entropy, y=cluster)) +
  geom_col(fill="#4DBBD5", width=0.7) +
  geom_vline(xintercept=mean_ent, linetype="dashed", color="#E64B35", linewidth=0.5) +
  labs(x="Mixing Entropy", y="Cluster",
       title="Patient Mixing Entropy per cNMF Cluster") +
  theme_pub(base_size=8) +
  theme(plot.title=element_text(size=9))

ggsave("fig_s3d_cluster_entropy.pdf", p_s3d, width=90, height=80, units="mm")
ggsave("fig_s3d_cluster_entropy.png", p_s3d, width=90, height=80, units="mm", dpi=300)
cat("  S3D saved\n")


## S3E  Hierarchical clustering dendrogram of 15 programs

cat("── S3E ──\n")

anno_15 <- fread("../fig2/cnmf_consensus_mp_annotation.csv")
anno_15 <- anno_15[match(prog_ids, anno_15$program_id)]

if(requireNamespace("dendextend", quietly=TRUE)) {
  library(dendextend)
  dend <- as.dendrogram(hc_15)
  leaf_mp <- anno_15$MP[order.dendrogram(dend)]
  leaf_col <- mp_colors[leaf_mp]
  leaf_col[is.na(leaf_col)] <- "grey60"
  labels_colors(dend) <- leaf_col
  labels_cex(dend) <- 0.8
  
  pdf("fig_s3e_dendrogram.pdf", width=6, height=3.5)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(dend, main="Hierarchical Clustering of cNMF Programs",
       ylab="Height", cex.main=0.9)
  legend("topright", legend=names(mp_colors)[1:4], fill=mp_colors[1:4],
         border=NA, bty="n", cex=0.7)
  dev.off()
  
  png("fig_s3e_dendrogram.png", width=6, height=3.5, units="in", res=300)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(dend, main="Hierarchical Clustering of cNMF Programs",
       ylab="Height", cex.main=0.9)
  legend("topright", legend=names(mp_colors)[1:4], fill=mp_colors[1:4],
         border=NA, bty="n", cex=0.7)
  dev.off()
} else {
  pdf("fig_s3e_dendrogram.pdf", width=6, height=3.5)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(hc_15, main="Hierarchical Clustering of cNMF Programs",
       xlab="", ylab="Height", cex=0.8)
  dev.off()
  png("fig_s3e_dendrogram.png", width=6, height=3.5, units="in", res=300)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(hc_15, main="Hierarchical Clustering of cNMF Programs",
       xlab="", ylab="Height", cex=0.8)
  dev.off()
}
cat("  S3E saved\n")


# S3F MP proportion across samples

cat("── S3F ──\n")

umap <- fread("../fig2/malignant_umap_metadata.csv.gz")
umap_mp <- umap[dominant_MP %in% c("MP1","MP2","MP3","MP4")]

sample_col <- intersect(c("sample_id","patient_id"), names(umap_mp))[1]
cat("  grouping by:", sample_col, "\n")

mp_by_sample <- umap_mp[, .N, by=c(sample_col, "dominant_MP")]
mp_by_sample[, total := sum(N), by=sample_col]
mp_by_sample[, pct := N / total]
mp_by_sample[, dominant_MP := factor(dominant_MP, levels=c("MP4","MP3","MP2","MP1"))]

mp1_pct <- mp_by_sample[dominant_MP=="MP1", .(mp1_pct=pct), by=sample_col]
sample_order <- mp1_pct[order(mp1_pct)][[sample_col]]
all_samples <- unique(mp_by_sample[[sample_col]])
sample_order <- c(sample_order, setdiff(all_samples, sample_order))
mp_by_sample[, (sample_col) := factor(get(sample_col), levels=sample_order)]

mp_labels_full <- c("MP1"="MP1","MP2"="MP2","MP3"="MP3","MP4"="MP4")

p_s3f <- ggplot(mp_by_sample, aes(x=get(sample_col), y=pct, fill=dominant_MP)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values=mp_colors, labels=mp_labels_full, name="Meta-Program") +
  scale_y_continuous(labels=percent, expand=c(0,0)) +
  labs(x=NULL, y="Proportion") +
  theme_pub(base_size=6) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=3.5),
        legend.position="right",
        legend.key.size=unit(3,"mm"))

ggsave("fig_s3f_mp_proportion.pdf", p_s3f, width=250, height=60, units="mm")
ggsave("fig_s3f_mp_proportion.png", p_s3f, width=250, height=60, units="mm", dpi=300)
cat("  S3F saved\n")


## S3G  cNMF cluster vs MP mixing entropy (boxplot)

cat("── S3G ──\n")

# MP entropy
mp_ent <- gep_anno[, {
  pt_counts <- table(patient_key)
  .(entropy=shannon_entropy(pt_counts))
}, by=MP]

ent_compare <- rbind(
  data.frame(level="cNMF", entropy=cluster_ent$entropy),
  data.frame(level="MP",   entropy=mp_ent$entropy)
)
ent_compare$level <- factor(ent_compare$level, levels=c("cNMF","MP"))

wt <- wilcox.test(cluster_ent$entropy, mp_ent$entropy)
cat("  Wilcoxon p:", wt$p.value, "\n")

p_s3g <- ggplot(ent_compare, aes(x=level, y=entropy, fill=level)) +
  geom_boxplot(width=0.5, outlier.size=0.8) +
  geom_jitter(width=0.15, size=1.5, alpha=0.6) +
  scale_fill_manual(values=c("cNMF"="#4DBBD5","MP"="#E64B35")) +
  annotate("text", x=1.5, y=max(ent_compare$entropy)*1.05,
           label=sprintf("p = %.2g", wt$p.value),
           size=3, family=my_font) +
  labs(x=NULL, y="Patient Mixing Entropy",
       title="Comparison of Mixing Entropy") +
  theme_pub(base_size=9) +
  theme(legend.position="none")

ggsave("fig_s3g_entropy_compare.pdf", p_s3g, width=70, height=65, units="mm")
ggsave("fig_s3g_entropy_compare.png", p_s3g, width=70, height=65, units="mm", dpi=300)
cat("  S3G saved\n")


message("=== Supp Fig 3 (cNMF QC) all 7 panels done ===")


# S3A / S3C / S3D

library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig_s2" else "${WORK_ROOT}/luad_figures/fig_s2")

# (override removed: keep Arial from top)

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}


# S3A y + per-patient

cat("── S3A ──\n")

kcurves <- fread("s2_k_error_curves.csv")
k_sel   <- fread("s2_k_selection_summary.csv")

# patient stability error [0,1
kcurves[, stab_norm := (stability - min(stability)) / (max(stability) - min(stability) + 1e-10),
        by=patient_key]
kcurves[, err_norm := (prediction_error - min(prediction_error)) / 
          (max(prediction_error) - min(prediction_error) + 1e-10),
        by=patient_key]

p_s3a <- ggplot(kcurves, aes(x=K)) +
  # Stability
  geom_line(aes(y=stab_norm, color="Stability"), linewidth=0.4) +
  geom_point(aes(y=stab_norm, color="Stability"), size=0.8) +
  # Error error = y
  geom_line(aes(y=err_norm, color="Error"), linewidth=0.4) +
  geom_point(aes(y=err_norm, color="Error"), size=0.8) +
  # selected K
  geom_vline(data=k_sel, aes(xintercept=selected_k),
             linetype="dashed", color="grey30", linewidth=0.3) +
  facet_wrap(~patient_key, scales="free_x", ncol=8) +
  scale_color_manual(values=c("Stability"="#3C5488","Error"="#E64B35"),
                     labels=c("Stability"="Stability","Error"="Error")) +
  scale_y_continuous(name="Stability", 
                     sec.axis=sec_axis(~., name="Error")) +
  labs(x="Number of Components") +
  theme_pub(base_size=6) +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(3,"mm"),
        strip.text=element_text(size=4.5),
        axis.text=element_text(size=3.5),
        axis.title.y.left=element_text(color="#3C5488", size=6),
        axis.title.y.right=element_text(color="#E64B35", size=6),
        panel.spacing=unit(0.5,"mm"))

ggsave("fig_s3a_k_curves.pdf", p_s3a, width=350, height=220, units="mm")
ggsave("fig_s3a_k_curves.png", p_s3a, width=350, height=220, units="mm", dpi=300)
cat("  S3A saved\n")


# S3C 15 cluster proportion 1-15 +

cat("── S3C ──\n")

# gep_anno + cutree
gep_anno <- fread("../fig2/gep_mp_annotation.csv")
gep_anno <- gep_anno[MP != "MP5"]

corr_770 <- fread("gep_spearman_corr.csv")
gep_ids_770 <- corr_770[[1]]
corr_mat_770 <- as.matrix(corr_770[, -1, with=FALSE])
rownames(corr_mat_770) <- gep_ids_770

keep <- intersect(gep_anno$gep_id, gep_ids_770)
corr_sub <- corr_mat_770[keep, keep]
hc_770 <- hclust(as.dist(1 - corr_sub), method="average")
clust_15 <- cutree(hc_770, k=15)

gep_anno <- gep_anno[gep_id %in% keep]
gep_anno[, cluster := paste0("cNMF_", clust_15[match(gep_id, names(clust_15))])]

clust_by_pt <- gep_anno[, .N, by=.(patient_key, cluster)]
clust_by_pt[, total := sum(N), by=patient_key]
clust_by_pt[, pct := N / total]

# cNMF_1, cNMF_2, ..., cNMF_15
clust_levels <- paste0("cNMF_", 1:15)
clust_levels <- intersect(clust_levels, unique(clust_by_pt$cluster))

clust_pal <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set3"))[1:length(clust_levels)]
names(clust_pal) <- clust_levels

clust_by_pt[, cluster := factor(cluster, levels=rev(clust_levels))]

# sample
dom_clust <- clust_by_pt[, .SD[which.max(pct)], by=patient_key]
pt_order <- dom_clust[order(cluster, -pct)]$patient_key
clust_by_pt[, patient_key := factor(patient_key, levels=pt_order)]

p_s3c <- ggplot(clust_by_pt, aes(x=patient_key, y=pct, fill=cluster)) +
  geom_bar(stat="identity", width=1, color="black", linewidth=0.1) +  # 
  scale_fill_manual(values=clust_pal, name="cNMF_cluster",
                    breaks=clust_levels) +   #  1-15 
  scale_y_continuous(labels=percent, expand=c(0,0)) +
  labs(x=NULL, y="Proportion") +
  theme_pub(base_size=6) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=3.5),
        legend.position="right",
        legend.key.size=unit(3,"mm"),
        legend.text=element_text(size=5))

ggsave("fig_s3c_cluster_proportion.pdf", p_s3c, width=250, height=65, units="mm")
ggsave("fig_s3c_cluster_proportion.png", p_s3c, width=250, height=65, units="mm", dpi=300)
cat("  S3C saved\n")


# S3D Entropy 0 cluster

cat("── S3D ──\n")

shannon_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log2(p))
}

cluster_ent <- gep_anno[, {
  pt_counts <- table(patient_key)
  .(entropy=shannon_entropy(pt_counts),
    n_patients=length(unique(patient_key)),
    n_geps=.N)
}, by=cluster]

# 15 cluster
all_clusters <- paste0("cNMF_", 1:15)
missing <- setdiff(all_clusters, cluster_ent$cluster)
if(length(missing) > 0) {
  cluster_ent <- rbind(cluster_ent,
                       data.table(cluster=missing, entropy=0, n_patients=0, n_geps=0))
}

cluster_ent <- cluster_ent[order(entropy, decreasing=TRUE)]
cluster_ent[, cluster := factor(cluster, levels=cluster)]

mean_ent <- mean(cluster_ent$entropy)

# n_patients
p_s3d <- ggplot(cluster_ent, aes(x=entropy, y=cluster)) +
  geom_col(fill="#4DBBD5", width=0.7) +
  geom_vline(xintercept=mean_ent, linetype="dashed", color="#E64B35", linewidth=0.5) +
  geom_text(aes(label=paste0("n=", n_patients)),
            hjust=-0.1, size=2, family=my_font) +
  scale_x_continuous(expand=expansion(mult=c(0, 0.15))) +
  labs(x="Mixing Entropy", y="Cluster",
       title="Patient Mixing Entropy per cNMF Cluster") +
  theme_pub(base_size=8) +
  theme(plot.title=element_text(size=9))

ggsave("fig_s3d_cluster_entropy.pdf", p_s3d, width=100, height=80, units="mm")
ggsave("fig_s3d_cluster_entropy.png", p_s3d, width=100, height=80, units="mm", dpi=300)
cat("  S3D saved\n")


message("=== S3A/C/D fix done ===")


# S3A v2

library(data.table)
library(ggplot2)

setwd(if (.Platform$OS.type == "windows") "${WORK_ROOT}/luad_figures/fig_s2" else "${WORK_ROOT}/luad_figures/fig_s2")
# (override removed: keep Arial from top)

kcurves <- fread("s2_k_error_curves.csv")

# patient [0,1
kcurves[, stab_01 := (stability - min(stability)) / 
          (max(stability) - min(stability) + 1e-10), by=patient_key]
kcurves[, err_01 := (prediction_error - min(prediction_error)) / 
          (max(prediction_error) - min(prediction_error) + 1e-10), by=patient_key]

summary_k <- kcurves[, .(
  stab_med  = median(stab_01),
  stab_q25  = quantile(stab_01, 0.25),
  stab_q75  = quantile(stab_01, 0.75),
  err_med   = median(err_01),
  err_q25   = quantile(err_01, 0.25),
  err_q75   = quantile(err_01, 0.75)
), by=K]

p <- ggplot(summary_k, aes(x=K)) +
  # Stability ribbon + line
  geom_ribbon(aes(ymin=stab_q25, ymax=stab_q75), fill="#3C5488", alpha=0.15) +
  geom_line(aes(y=stab_med, color="Stability"), linewidth=0.9) +
  geom_point(aes(y=stab_med, color="Stability"), size=2.5) +
  # Error ribbon + line
  geom_ribbon(aes(ymin=err_q25, ymax=err_q75), fill="#E64B35", alpha=0.15) +
  geom_line(aes(y=err_med, color="Error"), linewidth=0.9) +
  geom_point(aes(y=err_med, color="Error"), size=2.5) +
  scale_color_manual(values=c("Stability"="#3C5488","Error"="#E64B35")) +
  scale_x_continuous(breaks=unique(summary_k$K)) +
  scale_y_continuous(name="Stability (normalized)",
                     sec.axis=sec_axis(~., name="Error (normalized)")) +
  labs(x="Number of Components") +
  theme_classic(base_family=my_font, base_size=10) +
  theme(
    axis.text=element_text(color="black"),
    axis.title.y.left=element_text(color="#3C5488", face="bold"),
    axis.title.y.right=element_text(color="#E64B35", face="bold"),
    axis.text.y.left=element_text(color="#3C5488"),
    axis.text.y.right=element_text(color="#E64B35"),
    legend.position="top",
    legend.title=element_blank()
  )

ggsave("fig_s3a_summary.pdf", p, width=100, height=75, units="mm")
ggsave("fig_s3a_summary.png", p, width=100, height=75, units="mm", dpi=300)
cat("S3A summary v2 saved\n")
