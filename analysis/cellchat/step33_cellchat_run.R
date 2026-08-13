#!/usr/bin/env Rscript
# step33: run CellChat on 4 inputs (all + Normal/Tumor/Metastasis)
# Saves rds objects to ~/luad/data/processed/cellchat_<group>.rds
# Then plot Fig 6 panels A-F to ~/luad/results/fig6_panels/
# Sync the panels to ${WORK_ROOT}/luad_figures/fig6/panels/

suppressPackageStartupMessages({
  library(CellChat)
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Matrix)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(showtext)
  library(sysfonts)
})

DAT  <- "${PROJECT_ROOT}/data/processed"
RES  <- "${PROJECT_ROOT}/results"
OUT  <- file.path(RES, "fig6_panels")
RDS  <- file.path(DAT, "cellchat_rds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS, showWarnings = FALSE, recursive = TRUE)

arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial",
    regular   = arial_p,
    bold      = "~/.local/share/fonts/arialbd.ttf",
    italic    = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300)
  my_font <- "Arial"
} else { my_font <- "sans" }
cat("[font]", my_font, "\n")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)
# Use multicore (fork, shared memory) on Linux to avoid 1GB globals serialization
# fallback to sequential if not supported.
if (.Platform$OS.type == "unix" && !interactive()) {
  tryCatch(future::plan("multicore", workers = 4),
           error = function(e) future::plan("sequential"))
} else {
  future::plan("sequential")
}
cat("[future]", class(future::plan())[1], "\n")

neu_colors <- c(
  "Neu_Inflammatory" = "#E64B35", "Neu_Angiogenic" = "#F39B7F",
  "Neu_Metastatic" = "#3C5488", "Neu_ECM_remodeling" = "#4DBBD5",
  "Neu_OSM_priming" = "#00A087", "Neu_OSM_low" = "#8491B4",
  "Neu_IFN_response" = "#91D1C2"
)
mp_colors <- c("Mal_MP1" = "#E64B35", "Mal_MP2" = "#4DBBD5",
               "Mal_MP3" = "#00A087", "Mal_MP4" = "#3C5488")
macro_colors <- c(
  "Macro_C1QC" = "#4DBBD5", "Macro_FCN1" = "#E64B35",
  "Macro_FOLR2" = "#00A087", "Macro_MARCO" = "#3C5488",
  "Macro_SPP1" = "#F39B7F", "Macro_general" = "#8491B4",
  "Macro_prolif" = "#91D1C2"
)
other_colors <- c(
  "Fibroblast" = "#B09C85", "Endothelial" = "#7E57C2",
  "T_NK" = "#FFCB5C", "B" = "#9C27B0", "Plasma" = "#D9D9D9",
  "Mast" = "#F0A23B", "Epithelial_Normal" = "#D2B48C",
  "Mono_nonclassical" = "#999999",
  "cDC1" = "#80CBC4", "cDC2" = "#7570B3",
  "cDC_LAMP3" = "#984EA3", "pDC" = "#B2B2B2"
)
all_colors <- c(neu_colors, mp_colors, macro_colors, other_colors)

# ── helper: build CellChat object from h5ad ──
build_cellchat <- function(h5ad_path, label_col = "cellchat_label",
                            min_cells_per_group = 10) {
  cat("  reading", h5ad_path, "\n")
  sce <- readH5AD(h5ad_path, raw = FALSE, use_hdf5 = FALSE)
  expr <- assay(sce, "X")
  if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "CsparseMatrix")
  meta <- as.data.frame(colData(sce))
  meta$labels <- as.character(meta[[label_col]])
  # drop groups with too few cells
  tab <- table(meta$labels)
  keep_lbl <- names(tab)[tab >= min_cells_per_group]
  keep_idx <- meta$labels %in% keep_lbl
  expr <- expr[, keep_idx]
  meta <- meta[keep_idx, , drop = FALSE]
  cat("  after group-filter:", ncol(expr), "cells across",
      length(unique(meta$labels)), "groups\n")
  cat("  groups:\n"); print(sort(table(meta$labels)))
  cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE,
                                  population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}

groups <- list(
  all = file.path(DAT, "cellchat_input_all.h5ad"),
  Normal = file.path(DAT, "cellchat_input_Normal.h5ad"),
  Tumor = file.path(DAT, "cellchat_input_Tumor.h5ad"),
  Metastasis = file.path(DAT, "cellchat_input_Metastasis.h5ad")
)

cellchat_list <- list()
for (g in names(groups)) {
  rds_path <- file.path(RDS, sprintf("cellchat_%s.rds", g))
  if (file.exists(rds_path)) {
    cat(sprintf("[%s] loading cached %s\n", g, rds_path))
    cellchat_list[[g]] <- readRDS(rds_path)
  } else if (file.exists(groups[[g]])) {
    cat(sprintf("[%s] building from %s\n", g, groups[[g]]))
    t1 <- Sys.time()
    cc <- build_cellchat(groups[[g]])
    saveRDS(cc, rds_path)
    cellchat_list[[g]] <- cc
    cat(sprintf("[%s] done in %.1f min\n", g,
                as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  } else {
    cat(sprintf("[%s] missing input — skip\n", g))
  }
}

# Panel 6A — incoming vs outgoing scatter (Tumor + Metastasis)
cat("\n[6A] netAnalysis_signalingRole_scatter\n")
if (!is.null(cellchat_list$Tumor)) {
  pdf(file.path(OUT, "fig6a_scatter_tumor.pdf"), width = 5.5, height = 4.5)
  print(netAnalysis_signalingRole_scatter(cellchat_list$Tumor,
        title = NULL, font.size = 8, label.size = 2.5))
  dev.off()
  png(file.path(OUT, "fig6a_scatter_tumor.png"),
      width = 5.5, height = 4.5, units = "in", res = 300)
  print(netAnalysis_signalingRole_scatter(cellchat_list$Tumor,
        title = NULL, font.size = 8, label.size = 2.5))
  dev.off()
}
if (!is.null(cellchat_list$Metastasis)) {
  pdf(file.path(OUT, "fig6a_scatter_met.pdf"), width = 5.5, height = 4.5)
  print(netAnalysis_signalingRole_scatter(cellchat_list$Metastasis,
        title = NULL, font.size = 8, label.size = 2.5))
  dev.off()
  png(file.path(OUT, "fig6a_scatter_met.png"),
      width = 5.5, height = 4.5, units = "in", res = 300)
  print(netAnalysis_signalingRole_scatter(cellchat_list$Metastasis,
        title = NULL, font.size = 8, label.size = 2.5))
  dev.off()
}

# Panel 6B — differential heatmap (Tumor vs Metastasis)
cat("\n[6B] mergeCellChat + netVisual_heatmap (count + weight)\n")
if (!is.null(cellchat_list$Tumor) && !is.null(cellchat_list$Metastasis)) {
  cc_tumor <- cellchat_list$Tumor
  cc_met   <- cellchat_list$Metastasis
  # CellChat::mergeCellChat requires identical group sets across objects — lift to union
  union_groups <- sort(union(levels(cc_tumor@idents), levels(cc_met@idents)))
  cat("  union groups:", length(union_groups), "\n")
  cc_tumor <- tryCatch(liftCellChat(cc_tumor, group.new = union_groups),
                       error = function(e) { cat("  liftCellChat tumor failed:", conditionMessage(e), "\n"); cc_tumor })
  cc_met   <- tryCatch(liftCellChat(cc_met,   group.new = union_groups),
                       error = function(e) { cat("  liftCellChat met failed:", conditionMessage(e), "\n"); cc_met })
  obj_list <- list(Tumor = cc_tumor, Met = cc_met)
  cellchat_merged <- mergeCellChat(obj_list, add.names = names(obj_list),
                                    cell.prefix = FALSE)
  # 6B-1: differential interaction count
  pdf(file.path(OUT, "fig6b_diff_count.pdf"), width = 6.5, height = 6.5)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2),
                           measure = "count", title.name = NULL,
                           font.size = 6, font.size.title = 8,
                           color.heatmap = c("#3C5488", "#E64B35")))
  dev.off()
  png(file.path(OUT, "fig6b_diff_count.png"),
      width = 6.5, height = 6.5, units = "in", res = 300)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2),
                           measure = "count", title.name = NULL,
                           font.size = 6, font.size.title = 8,
                           color.heatmap = c("#3C5488", "#E64B35")))
  dev.off()
  # 6B-2: differential interaction strength
  pdf(file.path(OUT, "fig6b_diff_weight.pdf"), width = 6.5, height = 6.5)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2),
                           measure = "weight", title.name = NULL,
                           font.size = 6, font.size.title = 8,
                           color.heatmap = c("#3C5488", "#E64B35")))
  dev.off()
  png(file.path(OUT, "fig6b_diff_weight.png"),
      width = 6.5, height = 6.5, units = "in", res = 300)
  print(netVisual_heatmap(cellchat_merged, comparison = c(1, 2),
                           measure = "weight", title.name = NULL,
                           font.size = 6, font.size.title = 8,
                           color.heatmap = c("#3C5488", "#E64B35")))
  dev.off()
  saveRDS(cellchat_merged, file.path(RDS, "cellchat_merged_TM.rds"))
}

# Panel 6C — same as 6B but show raw weight side-by-side
cat("\n[6C] side-by-side weight heatmaps Normal vs Tumor vs Met\n")
if (length(cellchat_list) >= 2) {
  vis_groups <- intersect(c("Normal","Tumor","Metastasis"), names(cellchat_list))
  pdf(file.path(OUT, "fig6c_weight_per_tissue.pdf"),
      width = 6 * length(vis_groups), height = 5.5)
  par(mfrow = c(1, length(vis_groups)))
  for (g in vis_groups) {
    print(netVisual_heatmap(cellchat_list[[g]], measure = "weight",
                             title.name = g, font.size = 5,
                             font.size.title = 8,
                             color.heatmap = "Reds"))
  }
  dev.off()
  png(file.path(OUT, "fig6c_weight_per_tissue.png"),
      width = 6 * length(vis_groups), height = 5.5, units = "in", res = 300)
  par(mfrow = c(1, length(vis_groups)))
  for (g in vis_groups) {
    print(netVisual_heatmap(cellchat_list[[g]], measure = "weight",
                             title.name = g, font.size = 5,
                             font.size.title = 8,
                             color.heatmap = "Reds"))
  }
  dev.off()
}

# Panel 6D — Neu subtypes → Mal_MP* bubble (EMT pathways focus)
cat("\n[6D] netVisual_bubble Neu_* → Mal_MP*\n")
cc_focus <- cellchat_list$all
if (!is.null(cc_focus)) {
  ident_levels <- levels(cc_focus@idents)
  neu_present <- intersect(
    c("Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low",
      "Neu_ECM_remodeling","Neu_Angiogenic","Neu_Metastatic","Neu_IFN_response"),
    ident_levels)
  mal_present <- intersect(c("Mal_MP1","Mal_MP2","Mal_MP3","Mal_MP4"), ident_levels)
  emt_pw <- c("TGFb","IL1","TNF","OSM","CXCL","CCL","VEGF","FN1","SPP1","MMP","IL6","EGF")
  if (length(neu_present) > 0 && length(mal_present) > 0) {
    pdf(file.path(OUT, "fig6d_bubble_neu_to_mal.pdf"),
        width = 13, height = 10)
    print(netVisual_bubble(cc_focus,
        sources.use = neu_present, targets.use = mal_present,
        signaling = intersect(emt_pw, cc_focus@netP$pathways),
        sort.by.source = TRUE, font.size = 7, font.size.title = 9,
        angle.x = 45))
    dev.off()
    png(file.path(OUT, "fig6d_bubble_neu_to_mal.png"),
        width = 13, height = 10, units = "in", res = 300)
    print(netVisual_bubble(cc_focus,
        sources.use = neu_present, targets.use = mal_present,
        signaling = intersect(emt_pw, cc_focus@netP$pathways),
        sort.by.source = TRUE, font.size = 7, font.size.title = 9,
        angle.x = 45))
    dev.off()
  }

  # Panel 6E — Neu_Inflammatory only → Mal_MP1 vs MP3 contrast
  cat("\n[6E] netVisual_bubble Neu_Inflammatory → Mal_MP1 vs Mal_MP3\n")
  if ("Neu_Inflammatory" %in% ident_levels &&
      length(intersect(c("Mal_MP1","Mal_MP3"), ident_levels)) == 2) {
    pdf(file.path(OUT, "fig6e_bubble_inflam_mp1_vs_mp3.pdf"),
        width = 9, height = 9)
    print(netVisual_bubble(cc_focus,
        sources.use = "Neu_Inflammatory",
        targets.use = c("Mal_MP1", "Mal_MP3"),
        signaling = intersect(emt_pw, cc_focus@netP$pathways),
        sort.by.source = TRUE, font.size = 7, font.size.title = 9,
        angle.x = 45))
    dev.off()
    png(file.path(OUT, "fig6e_bubble_inflam_mp1_vs_mp3.png"),
        width = 9, height = 9, units = "in", res = 300)
    print(netVisual_bubble(cc_focus,
        sources.use = "Neu_Inflammatory",
        targets.use = c("Mal_MP1", "Mal_MP3"),
        signaling = intersect(emt_pw, cc_focus@netP$pathways),
        sort.by.source = TRUE, font.size = 7, font.size.title = 9,
        angle.x = 45))
    dev.off()
  }

  # Panel 6F — bar: per-sender total EMT pathway prob (Neu vs Macro vs Fib)
  cat("\n[6F] EMT pathway sender-strength bar\n")
  emt_pw_present <- intersect(emt_pw, cc_focus@netP$pathways)
  net <- subsetCommunication(cc_focus, signaling = emt_pw_present)
  if (nrow(net) > 0) {
    senders_focus <- intersect(c(
      neu_present,
      "Macro_SPP1","Macro_FCN1","Macro_C1QC","Macro_general",
      "Macro_FOLR2","Macro_MARCO","Macro_prolif",
      "Fibroblast", "Endothelial"), levels(cc_focus@idents))
    mal_target <- intersect(c("Mal_MP1","Mal_MP3"), levels(cc_focus@idents))
    sub_net <- net %>%
      filter(source %in% senders_focus, target %in% mal_target) %>%
      group_by(source) %>%
      summarise(total_prob = sum(prob, na.rm = TRUE),
                n_pairs = n()) %>%
      arrange(desc(total_prob))
    write.csv(sub_net, file.path(OUT, "fig6f_emt_pathway_strength.csv"), row.names = FALSE)
    sub_net$class <- ifelse(sub_net$source %in% neu_present, "Neutrophil",
                     ifelse(grepl("^Macro_", sub_net$source), "Macrophage",
                     ifelse(sub_net$source == "Fibroblast", "Fibroblast",
                     ifelse(sub_net$source == "Endothelial", "Endothelial", "Other"))))
    sub_net$source <- factor(sub_net$source, levels = sub_net$source)
    pal_fill <- c("Neutrophil" = "#E64B35", "Macrophage" = "#3C5488",
                  "Fibroblast" = "#B09C85", "Endothelial" = "#7E57C2",
                  "Other" = "grey60")
    p6f <- ggplot(sub_net, aes(source, total_prob, fill = class)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = pal_fill, name = NULL) +
      coord_flip() +
      labs(x = NULL, y = "EMT pathway communication prob.\n(sum, sender → Mal_MP1+MP3)") +
      theme_classic(base_size = 8, base_family = my_font) +
      theme(axis.text = element_text(color = "black"),
            axis.line = element_line(linewidth = 0.4, color = "black"),
            plot.title = element_blank(),
            legend.position = "right",
            legend.key.size = unit(3, "mm"))
    ggsave(file.path(OUT, "fig6f_emt_pathway_bar.pdf"), p6f,
           width = 130, height = 110, units = "mm")
    ggsave(file.path(OUT, "fig6f_emt_pathway_bar.png"), p6f,
           width = 130, height = 110, units = "mm", dpi = 300)
  }
}

cat("\n=== DONE ===\n")
print(list.files(OUT))
