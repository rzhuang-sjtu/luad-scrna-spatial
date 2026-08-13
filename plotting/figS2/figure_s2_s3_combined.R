# Figure S2 | cNMF QC (7 panels

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)

setwd("${WORK_ROOT}/luad_figures/fig_s2")

my_font <- "sans"
mp_colors <- c("MP1"="#E64B35","MP2"="#4DBBD5","MP3"="#00A087","MP4"="#3C5488")

theme_pub <- function(base_size=10) {
  theme_classic(base_family=my_font, base_size=base_size) +
    theme(axis.text=element_text(color="black"),
          plot.title=element_text(face="bold", size=base_size+1))
}

shannon_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log2(p))
}


# S2A K selection curves (per-patient

cat("── S2A ──\n")

kcurves <- fread("s2_k_error_curves.csv")
k_sel   <- fread("s2_k_selection_summary.csv")

p_s2a <- ggplot(kcurves, aes(x=K)) +
  geom_line(aes(y=stability, color="Stability"), linewidth=0.4) +
  geom_point(aes(y=stability, color="Stability"), size=0.6) +
  geom_line(aes(y=prediction_error / max(prediction_error, na.rm=TRUE),
                color="Prediction Error"), linewidth=0.4) +
  geom_point(aes(y=prediction_error / max(prediction_error, na.rm=TRUE),
                 color="Prediction Error"), size=0.6) +
  geom_vline(data=k_sel, aes(xintercept=selected_k),
             linetype="dashed", color="grey40", linewidth=0.2) +
  facet_wrap(~patient_key, scales="free_y", ncol=8) +
  scale_color_manual(values=c("Stability"="#E64B35","Prediction Error"="#3C5488")) +
  labs(x="Number of Components (K)", y="Stability / Scaled Error") +
  theme_pub(base_size=6) +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.key.size=unit(3,"mm"),
        strip.text=element_text(size=4.5),
        axis.text=element_text(size=4),
        panel.spacing=unit(0.5,"mm"))

ggsave("fig_s2a_k_curves.pdf", p_s2a, width=350, height=220, units="mm")
ggsave("fig_s2a_k_curves.png", p_s2a, width=350, height=220, units="mm", dpi=300)
cat("  S2A saved\n")


## S2B  K histogram

cat("── S2B ──\n")

khist <- fread("s2_k_histogram.csv")

p_s2b <- ggplot(khist, aes(x=factor(K), y=count)) +
  geom_bar(stat="identity", fill="#3C5488", width=0.7) +
  geom_text(aes(label=count), vjust=-0.3, size=2.5, family=my_font) +
  labs(x="Selected K", y="Number of Patients") +
  theme_pub(base_size=9)

ggsave("fig_s2b_k_histogram.pdf", p_s2b, width=70, height=65, units="mm")
ggsave("fig_s2b_k_histogram.png", p_s2b, width=70, height=65, units="mm", dpi=300)
cat("  S2B saved\n")


## S2C  Stability vs n_cells scatter

cat("── S2C ──\n")

stab <- fread("s2_stability_vs_ncells.csv")
rho <- cor(stab$n_cells, stab$stability, method="spearman")

p_s2c <- ggplot(stab, aes(x=n_cells, y=stability)) +
  geom_point(aes(color=factor(selected_k)), size=1.5, alpha=0.8) +
  geom_smooth(method="lm", se=TRUE, color="grey30", linewidth=0.5,
              linetype="dashed") +
  scale_color_brewer(palette="Set2", name="Selected K") +
  scale_x_log10(labels=comma) +
  labs(x="Number of Malignant Cells (log)", y="Stability") +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
           label=sprintf("Spearman \u03C1 = %.2f", rho),
           size=3, family=my_font, fontface="italic") +
  theme_pub(base_size=8) +
  theme(legend.position="right",
        legend.key.size=unit(3,"mm"))

ggsave("fig_s2c_stability_vs_ncells.pdf", p_s2c, width=90, height=65, units="mm")
ggsave("fig_s2c_stability_vs_ncells.png", p_s2c, width=90, height=65, units="mm", dpi=300)
cat("  S2C saved\n")


# S2D Patient mixing entropy per MP (, MP5

cat("── S2D ──\n")

mp_mix <- fread("mp_patient_mixing.csv")
mp_mix <- mp_mix[MP != "MP5"]
mp_mix <- mp_mix[order(entropy)]
mp_mix[, MP := factor(MP, levels=MP)]

p_s2d <- ggplot(mp_mix, aes(x=entropy, y=MP, fill=MP)) +
  geom_col(width=0.6) +
  scale_fill_manual(values=mp_colors) +
  labs(x="Mixing Entropy", y=NULL,
       title="Patient Mixing Entropy per Meta-Program") +
  theme_pub(base_size=9) +
  theme(legend.position="none",
        plot.title=element_text(size=9))

ggsave("fig_s2d_entropy.pdf", p_s2d, width=90, height=50, units="mm")
ggsave("fig_s2d_entropy.png", p_s2d, width=90, height=50, units="mm", dpi=300)
cat("  S2D saved\n")


## S2E  Hierarchical clustering dendrogram of 15 consensus cNMF

cat("── S2E ──\n")

corr <- fread("../fig2/cnmf_consensus_corr.csv")
prog_ids <- corr[[1]]
corr_mat <- as.matrix(corr[, -1, with=FALSE])
rownames(corr_mat) <- prog_ids

anno <- fread("../fig2/cnmf_consensus_mp_annotation.csv")
anno <- anno[match(prog_ids, anno$program_id)]

dist_mat <- as.dist(1 - corr_mat)
hc <- hclust(dist_mat, method="average")

if(requireNamespace("dendextend", quietly=TRUE)) {
  library(dendextend)
  dend <- as.dendrogram(hc)
  leaf_mp <- anno$MP[order.dendrogram(dend)]
  leaf_col <- mp_colors[leaf_mp]
  leaf_col[is.na(leaf_col)] <- "grey60"
  labels_colors(dend) <- leaf_col
  labels_cex(dend) <- 0.8
  
  pdf("fig_s2e_dendrogram.pdf", width=6, height=3.5)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(dend, main="Hierarchical Clustering of cNMF Programs",
       ylab="Height (1 - Spearman \u03C1)", cex.main=0.9)
  legend("topright", legend=names(mp_colors), fill=mp_colors,
         border=NA, bty="n", cex=0.7)
  dev.off()
  
  png("fig_s2e_dendrogram.png", width=6, height=3.5, units="in", res=300)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(dend, main="Hierarchical Clustering of cNMF Programs",
       ylab="Height (1 - Spearman \u03C1)", cex.main=0.9)
  legend("topright", legend=names(mp_colors), fill=mp_colors,
         border=NA, bty="n", cex=0.7)
  dev.off()
} else {
  pdf("fig_s2e_dendrogram.pdf", width=6, height=3.5)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(hc, main="Hierarchical Clustering of cNMF Programs",
       xlab="", ylab="Height", cex=0.8)
  dev.off()
  png("fig_s2e_dendrogram.png", width=6, height=3.5, units="in", res=300)
  par(mar=c(4, 3, 2, 1), family=my_font)
  plot(hc, main="Hierarchical Clustering of cNMF Programs",
       xlab="", ylab="Height", cex=0.8)
  dev.off()
}
cat("  S2E saved\n")


# S2F MP proportion across samples

cat("── S2F ──\n")

umap <- fread("../fig2/malignant_umap_metadata.csv.gz")
umap_mp <- umap[dominant_MP %in% c("MP1","MP2","MP3","MP4")]

sample_col <- intersect(c("sample_id","patient_id"), names(umap_mp))[1]
cat("  grouping by:", sample_col, "\n")

mp_by_sample <- umap_mp[, .N, by=c(sample_col, "dominant_MP")]
mp_by_sample[, total := sum(N), by=sample_col]
mp_by_sample[, pct := N / total]
mp_by_sample[, dominant_MP := factor(dominant_MP,
                                     levels=c("MP4","MP3","MP2","MP1"))]

mp1_pct <- mp_by_sample[dominant_MP=="MP1", .(mp1_pct=pct), by=sample_col]
sample_order <- mp1_pct[order(mp1_pct)][[sample_col]]
all_samples <- unique(mp_by_sample[[sample_col]])
sample_order <- c(sample_order, setdiff(all_samples, sample_order))
mp_by_sample[, (sample_col) := factor(get(sample_col), levels=sample_order)]

mp_labels_full <- c("MP1"="MP1: Stress/AP-1","MP2"="MP2: Proliferative",
                    "MP3"="MP3: EMT/IFN","MP4"="MP4: AT2-like")

p_s2f <- ggplot(mp_by_sample, aes(x=get(sample_col), y=pct, fill=dominant_MP)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values=mp_colors, labels=mp_labels_full, name="Meta-Program") +
  scale_y_continuous(labels=percent, expand=c(0,0)) +
  labs(x=NULL, y="Proportion",
       title="Proportion of Meta-Program Clusters across Samples") +
  theme_pub(base_size=6) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=4),
        legend.position="right",
        legend.key.size=unit(3,"mm"),
        plot.title=element_text(size=8))

ggsave("fig_s2f_mp_proportion_samples.pdf", p_s2f, width=250, height=60, units="mm")
ggsave("fig_s2f_mp_proportion_samples.png", p_s2f, width=250, height=60, units="mm", dpi=300)
cat("  S2F saved\n")


# S2G Consensus cluster vs MP entropy (boxplot, MP5

cat("── S2G ──\n")

gep_anno <- fread("../fig2/gep_mp_annotation.csv")
gep_anno <- gep_anno[MP != "MP5"]

# 770 GEP 15 cluster
corr_full <- fread("gep_spearman_corr.csv")
gep_ids_full <- corr_full[[1]]
corr_mat_full <- as.matrix(corr_full[, -1, with=FALSE])
rownames(corr_mat_full) <- gep_ids_full

# MP5 GEP
keep_gep <- gep_anno$gep_id
corr_sub <- corr_mat_full[keep_gep, keep_gep]
hc_full <- hclust(as.dist(1 - corr_sub), method="average")
clusters_15 <- cutree(hc_full, k=15)

gep_anno[, consensus_cluster := paste0("C", clusters_15[match(gep_id, names(clusters_15))])]

# consensus cluster entropy
cluster_ent <- gep_anno[, {
  pt_counts <- table(patient_key)
  .(entropy=shannon_entropy(pt_counts),
    n_patients=length(unique(patient_key)),
    n_geps=.N)
}, by=consensus_cluster]

cat("  clusters:", nrow(cluster_ent), "\n")
cat("  cluster entropy range:", round(range(cluster_ent$entropy), 2), "\n")

# MP entropy
mp_ent <- gep_anno[, {
  pt_counts <- table(patient_key)
  .(entropy=shannon_entropy(pt_counts),
    n_patients=length(unique(patient_key)),
    n_geps=.N)
}, by=MP]

cat("  MP entropy range:", round(range(mp_ent$entropy), 2), "\n")

ent_compare <- rbind(
  data.frame(level="cNMF", entropy=cluster_ent$entropy),
  data.frame(level="MP",   entropy=mp_ent$entropy)
)
ent_compare$level <- factor(ent_compare$level, levels=c("cNMF","MP"))

wt <- wilcox.test(cluster_ent$entropy, mp_ent$entropy)
cat("  Wilcoxon p:", wt$p.value, "\n")

p_s2g <- ggplot(ent_compare, aes(x=level, y=entropy, fill=level)) +
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

ggsave("fig_s2g_entropy_compare.pdf", p_s2g, width=70, height=65, units="mm")
ggsave("fig_s2g_entropy_compare.png", p_s2g, width=70, height=65, units="mm", dpi=300)
cat("  S2G saved\n")


message("=== Fig S2 all 7 panels done ===")