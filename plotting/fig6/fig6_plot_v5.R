#!/usr/bin/env Rscript
# fig6_plot_v5.R — final pass:
#   6A: sqrt-axis (kept from v4)
#   6B/6C: hard-subset @net matrices to top-N active groups, then netVisual_heatmap
#   6D: square cells (kept from v4)
#   6F: free_y + fresh palette (kept from v4)

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
  library(showtext)
  library(sysfonts)
  library(grid)
  library(scales)
})

DAT <- "${PROJECT_ROOT}/data/processed"
RES <- "${PROJECT_ROOT}/results"
OUT <- file.path(RES, "fig6_panels")
RDS <- file.path(DAT, "cellchat_rds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

arial_p <- "~/.local/share/fonts/arial.ttf"
if (file.exists(path.expand(arial_p))) {
  sysfonts::font_add("Arial", regular = arial_p,
    bold = "~/.local/share/fonts/arialbd.ttf",
    italic = "~/.local/share/fonts/ariali.ttf")
  showtext_auto(); showtext_opts(dpi = 300)
  my_font <- "Arial"
} else { my_font <- "sans" }

neu_colors <- c(
  "Neu_Inflammatory"="#E64B35","Neu_Angiogenic"="#F39B7F",
  "Neu_Metastatic"="#3C5488","Neu_ECM_remodeling"="#4DBBD5",
  "Neu_OSM_priming"="#00A087","Neu_OSM_low"="#8491B4",
  "Neu_IFN_response"="#91D1C2")
mp_colors <- c("Mal_MP1"="#E64B35","Mal_MP2"="#4DBBD5",
                "Mal_MP3"="#00A087","Mal_MP4"="#3C5488")
macro_colors <- c("Macro_C1QC"="#4DBBD5","Macro_FCN1"="#E64B35",
  "Macro_FOLR2"="#00A087","Macro_MARCO"="#3C5488","Macro_SPP1"="#F39B7F",
  "Macro_general"="#8491B4","Macro_prolif"="#91D1C2")
other_colors <- c("Fibroblast"="#B09C85","Endothelial"="#7E57C2",
  "T_NK"="#FFCB5C","B"="#9C27B0","Plasma"="#D9D9D9","Mast"="#F0A23B",
  "Epithelial_Normal"="#D2B48C","Mono_nonclassical"="#999999",
  "cDC1"="#80CBC4","cDC2"="#7570B3","cDC_LAMP3"="#984EA3","pDC"="#B2B2B2")
all_colors <- c(neu_colors, mp_colors, macro_colors, other_colors)

cellchat_list <- list()
for (g in c("all","Normal","Tumor","Metastasis")) {
  p <- file.path(RDS, sprintf("cellchat_%s.rds", g))
  if (file.exists(p)) cellchat_list[[g]] <- readRDS(p)
}

# ── helper: hard-subset a CellChat object's @net to a subset of groups ──
subset_cc_groups <- function(cc, groups) {
  groups <- intersect(groups, levels(cc@idents))
  for (slot_name in c("net","netP")) {
    if (.hasSlot(cc, slot_name)) {
      sl <- slot(cc, slot_name)
      if (!is.null(sl$weight)) sl$weight <- sl$weight[groups, groups]
      if (!is.null(sl$count))  sl$count  <- sl$count[groups, groups]
      if (!is.null(sl$sum))    sl$sum    <- sl$sum[groups, groups]
      if (!is.null(sl$prob)) {
        if (length(dim(sl$prob)) == 3) sl$prob <- sl$prob[groups, groups, , drop = FALSE]
        else                            sl$prob <- sl$prob[groups, groups]
      }
      if (!is.null(sl$pval)) {
        if (length(dim(sl$pval)) == 3) sl$pval <- sl$pval[groups, groups, , drop = FALSE]
        else                            sl$pval <- sl$pval[groups, groups]
      }
      slot(cc, slot_name) <- sl
    }
  }
  # update idents — works for both single and merged objects
  if (is.list(cc@idents)) {
    if (!is.null(cc@idents$joint)) {
      cc@idents$joint <- factor(as.character(cc@idents$joint),
                                 levels = groups)
    }
    for (k in names(cc@idents)) {
      if (k != "joint" && is.factor(cc@idents[[k]])) {
        cc@idents[[k]] <- factor(as.character(cc@idents[[k]]),
                                  levels = groups)
      }
    }
  } else if (is.factor(cc@idents)) {
    cc@idents <- factor(as.character(cc@idents), levels = groups)
  }
  cc
}

# ── helper: hard-subset a merged CellChat object (per-condition net lists) ──
subset_merged_groups <- function(cc_merged, groups) {
  for (cond in names(cc_merged@net)) {
    sl <- cc_merged@net[[cond]]
    if (!is.null(sl$weight)) sl$weight <- sl$weight[groups, groups]
    if (!is.null(sl$count))  sl$count  <- sl$count[groups, groups]
    if (!is.null(sl$sum))    sl$sum    <- sl$sum[groups, groups]
    if (!is.null(sl$prob)) {
      if (length(dim(sl$prob)) == 3) sl$prob <- sl$prob[groups, groups, , drop = FALSE]
      else                            sl$prob <- sl$prob[groups, groups]
    }
    if (!is.null(sl$pval)) {
      if (length(dim(sl$pval)) == 3) sl$pval <- sl$pval[groups, groups, , drop = FALSE]
      else                            sl$pval <- sl$pval[groups, groups]
    }
    cc_merged@net[[cond]] <- sl
  }
  for (cond in names(cc_merged@netP)) {
    sl <- cc_merged@netP[[cond]]
    if (!is.null(sl$prob)) {
      if (length(dim(sl$prob)) == 3) sl$prob <- sl$prob[groups, groups, , drop = FALSE]
      else                            sl$prob <- sl$prob[groups, groups]
    }
    cc_merged@netP[[cond]] <- sl
  }
  if (!is.null(cc_merged@idents$joint)) {
    cc_merged@idents$joint <- factor(as.character(cc_merged@idents$joint),
                                      levels = groups)
  }
  for (cond in setdiff(names(cc_merged@idents), "joint")) {
    if (is.factor(cc_merged@idents[[cond]])) {
      cc_merged@idents[[cond]] <- factor(as.character(cc_merged@idents[[cond]]),
                                          levels = groups)
    }
  }
  cc_merged
}

# 6A — sqrt-axis scatter
#       v6 fix: dramatic label de-clutter + categorical colour by lineage
#       Only top-N by total signal gets a label; everything else stays as
#       a small grey dot. ggrepel with curved segments + larger spacing.
cat("\n[6A] sqrt-axis scatter (de-cluttered v6)\n")

# Lineage palette: warm = Neu, blue = Mal, blue-cyan = Macro, purple = Stromal,
# yellow = T_NK / DC, grey = everything else
lineage_pal <- c(
  "Neu"     = "#E64B35",
  "Mal"     = "#3C5488",
  "Macro"   = "#4DBBD5",
  "Stromal" = "#7E57C2",
  "T_NK"    = "#FFCB5C",
  "DC"      = "#984EA3",
  "Other"   = "grey60"
)
classify_lineage <- function(x) {
  ifelse(grepl("^Neu_", x), "Neu",
  ifelse(grepl("^Mal_", x), "Mal",
  ifelse(grepl("^Macro_", x), "Macro",
  ifelse(x %in% c("Fibroblast", "Endothelial"), "Stromal",
  ifelse(x %in% c("T_NK"), "T_NK",
  ifelse(grepl("^cDC|^pDC", x), "DC", "Other"))))))
}

draw_scatter_sqrt <- function(cc, fname, n_label = 14) {
  W <- cc@net$weight; C <- cc@net$count
  df <- data.frame(group = rownames(W),
                   out = rowSums(W), inc = colSums(W),
                   count = rowSums(C) + colSums(C))
  df <- df[df$count > 0, ]
  df$lineage <- factor(classify_lineage(df$group),
                       levels = names(lineage_pal))
  # Total signal = out + inc (used both for label selection and bubble size)
  df$total_signal <- df$out + df$inc
  # Label only top-N by total signal; force-keep all Neu_* and Mal_* always
  always_label <- df$group[grepl("^(Neu_|Mal_)", df$group)]
  top_label    <- df$group[order(-df$total_signal)][1:min(n_label, nrow(df))]
  df$show_label <- df$group %in% union(always_label, top_label)
  df$lab <- ifelse(df$show_label, df$group, NA_character_)

  p <- ggplot(df, aes(out, inc)) +
    geom_point(aes(size = count, fill = lineage),
               shape = 21, color = "black", stroke = 0.3, alpha = 0.92) +
    ggrepel::geom_text_repel(
      data = subset(df, show_label),
      aes(label = lab, color = lineage),
      size = 2.5, family = my_font, fontface = "bold",
      max.overlaps = Inf,
      box.padding = unit(0.55, "lines"),
      point.padding = unit(0.35, "lines"),
      min.segment.length = 0.15,
      segment.size = 0.2, segment.alpha = 0.6,
      segment.curvature = -0.15, segment.angle = 25,
      segment.ncp = 4,
      force = 4, force_pull = 0.4,
      seed = 42,
      bg.color = "white", bg.r = 0.12,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = lineage_pal, name = "Lineage", drop = FALSE,
                      guide = guide_legend(override.aes = list(size = 4),
                                           ncol = 1)) +
    scale_color_manual(values = lineage_pal, guide = "none") +
    scale_size_continuous(range = c(2, 8), name = "LR count",
                          breaks = pretty(df$count, 4)) +
    scale_x_continuous(trans = "sqrt",
                       breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08),
                       expand = expansion(mult = c(0.04, 0.08))) +
    scale_y_continuous(trans = "sqrt",
                       breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08),
                       expand = expansion(mult = c(0.04, 0.08))) +
    labs(x = "Outgoing interaction strength (√-scale)",
         y = "Incoming interaction strength (√-scale)") +
    theme_classic(base_size = 9, base_family = my_font) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(size = 9),
          axis.line = element_line(linewidth = 0.4),
          legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7),
          legend.key.size = unit(3, "mm"),
          plot.margin = margin(4, 6, 4, 4, "pt"))
  ggsave(file.path(OUT, paste0(fname, ".pdf")), p,
         width = 160, height = 130, units = "mm")
  ggsave(file.path(OUT, paste0(fname, ".png")), p,
         width = 160, height = 130, units = "mm", dpi = 300)
}
if (!is.null(cellchat_list$Tumor))      draw_scatter_sqrt(cellchat_list$Tumor,      "fig6a_scatter_tumor")
if (!is.null(cellchat_list$Metastasis)) draw_scatter_sqrt(cellchat_list$Metastasis, "fig6a_scatter_met")

# 6B — V1 style netVisual_heatmap, BUT with hard-subsetted @net to top-15 groups
cat("\n[6B] netVisual_heatmap on subsetted merged object\n")
if (!is.null(cellchat_list$Tumor) && !is.null(cellchat_list$Metastasis)) {
  cc_t <- cellchat_list$Tumor
  cc_m <- cellchat_list$Metastasis
  union_g <- sort(union(levels(cc_t@idents), levels(cc_m@idents)))
  cc_t_l <- liftCellChat(cc_t, group.new = union_g)
  cc_m_l <- liftCellChat(cc_m, group.new = union_g)
  diffW <- cc_m_l@net$weight - cc_t_l@net$weight
  diffC <- cc_m_l@net$count  - cc_t_l@net$count
  totW <- rowSums(abs(diffW)) + colSums(abs(diffW))
  top_g <- names(sort(totW, decreasing = TRUE))[1:min(15, length(totW))]
  top_g <- sort(top_g)
  cat(sprintf("  top-15: %s\n", paste(top_g, collapse=", ")))

  cellchat_merged <- mergeCellChat(list(Tumor = cc_t_l, Met = cc_m_l),
                                    add.names = c("Tumor","Met"),
                                    cell.prefix = FALSE)
  cellchat_merged_sub <- subset_merged_groups(cellchat_merged, top_g)

  for (m in c("count","weight")) {
    pdf(file.path(OUT, sprintf("fig6b_diff_%s.pdf", m)), width = 6, height = 6)
    print(netVisual_heatmap(cellchat_merged_sub, comparison = c(1, 2),
            measure = m, title.name = NULL,
            font.size = 8, font.size.title = 10,
            color.heatmap = c("#3C5488", "#E64B35")))
    dev.off()
    png(file.path(OUT, sprintf("fig6b_diff_%s.png", m)),
        width = 6, height = 6, units = "in", res = 300)
    print(netVisual_heatmap(cellchat_merged_sub, comparison = c(1, 2),
            measure = m, title.name = NULL,
            font.size = 8, font.size.title = 10,
            color.heatmap = c("#3C5488", "#E64B35")))
    dev.off()
  }
}

# 6C — per-tissue netVisual_heatmap on hard-subsetted top-15 groups
cat("\n[6C] per-tissue netVisual_heatmap (hard-subset top-15)\n")
for (g_old in c("Normal","Tumor","Metastasis")) {
  for (ext in c("pdf","png")) {
    f <- file.path(OUT, sprintf("fig6c_weight_%s.%s", g_old, ext))
    if (file.exists(f)) file.remove(f)
  }
}
for (g in intersect(c("Normal","Tumor","Metastasis"), names(cellchat_list))) {
  cc <- cellchat_list[[g]]
  W <- cc@net$weight
  diag(W) <- 0  # exclude self-loops from ranking AND display (they dominate the color scale)
  totals <- rowSums(W) + colSums(W)
  top <- names(sort(totals, decreasing = TRUE))[1:min(15, length(totals))]
  top <- sort(top)
  cat(sprintf("  [%s] top-15: %s\n", g, paste(top, collapse=", ")))
  cc_sub <- subset_cc_groups(cc, top)
  # zero out diagonal in the subsetted object so color scale is set by inter-group signal
  diag(cc_sub@net$weight) <- 0
  diag(cc_sub@net$count)  <- 0
  pdf(file.path(OUT, sprintf("fig6c_weight_%s.pdf", g)),
      width = 5.5, height = 5.5)
  print(netVisual_heatmap(cc_sub, measure = "weight",
            title.name = g, font.size = 8, font.size.title = 11,
            color.heatmap = "Reds"))
  dev.off()
  png(file.path(OUT, sprintf("fig6c_weight_%s.png", g)),
      width = 5.5, height = 5.5, units = "in", res = 300)
  print(netVisual_heatmap(cc_sub, measure = "weight",
            title.name = g, font.size = 8, font.size.title = 11,
            color.heatmap = "Reds"))
  dev.off()
}

# 6D — square bubble (unchanged from v4)
cat("\n[6D] square bubble\n")
cc_focus <- cellchat_list$all
if (!is.null(cc_focus)) {
  ident_levels <- levels(cc_focus@idents)
  neu_present <- intersect(c("Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low",
    "Neu_ECM_remodeling","Neu_Angiogenic","Neu_Metastatic","Neu_IFN_response"),
    ident_levels)
  mal_present <- intersect(c("Mal_MP1","Mal_MP2","Mal_MP3","Mal_MP4"), ident_levels)
  emt_pw <- c("TGFb","IL1","TNF","OSM","CCL","VEGF","FN1","SPP1","EGF")
  emt_pw_present <- intersect(emt_pw, cc_focus@netP$pathways)
  if (length(neu_present) > 0 && length(mal_present) > 0) {
    net_lr <- subsetCommunication(cc_focus,
                  sources.use = neu_present, targets.use = mal_present,
                  signaling = emt_pw_present)
    if (nrow(net_lr) > 0) {
      net_lr$pair <- paste(net_lr$source, net_lr$target, sep = " → ")

      # Drop pairs / LR with no significant signal — these are the rows /
      # cols that previously rendered as huge swathes of empty grey tiles.
      net_sig <- net_lr %>% filter(pval < 0.05 & prob > 0)
      keep_pairs <- net_sig %>% group_by(pair) %>%
        summarise(n = n()) %>% filter(n >= 1) %>% pull(pair)
      keep_lr    <- net_sig %>% group_by(interaction_name_2) %>%
        summarise(n = n()) %>% filter(n >= 1) %>% pull(interaction_name_2)
      net_show <- net_lr %>%
        filter(pair %in% keep_pairs, interaction_name_2 %in% keep_lr,
               !is.na(prob))

      pair_order <- intersect(
        as.vector(outer(neu_present, mal_present, paste, sep = " → ")),
        keep_pairs)
      inter_order <- net_show %>% group_by(interaction_name_2) %>%
        summarise(s = sum(prob, na.rm = TRUE)) %>%
        arrange(desc(s)) %>% pull(interaction_name_2)

      net_show$pair <- factor(net_show$pair, levels = pair_order)
      net_show$interaction_name_2 <- factor(net_show$interaction_name_2,
                                             levels = rev(inter_order))
      net_show$mlog_p <- -log10(pmax(net_show$pval, 1e-6))

      n_x <- length(pair_order); n_y <- length(inter_order)
      cat(sprintf("  6D: %d pairs × %d LR after filter (was %d × %d full grid)\n",
                  n_x, n_y, length(keep_pairs), length(keep_lr)))

      p6d <- ggplot(net_show, aes(pair, interaction_name_2)) +
        geom_point(aes(size = mlog_p, color = prob)) +
        scale_color_gradientn(colors = c("#4575B4","#FDDA82","#D73027"),
                              name = "Commun.\nProb.",
                              guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                     barheight = unit(20, "mm"),
                                                     frame.colour = "black",
                                                     frame.linewidth = 0.3,
                                                     ticks.colour = "black")) +
        scale_size_continuous(range = c(0.6, 6),
                              name = expression(-log[10]~italic(P)),
                              breaks = c(2, 3, 4, 5)) +
        coord_fixed(ratio = 1) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 9, base_family = my_font) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                         color = "black", size = 7),
              axis.text.y = element_text(color = "black", size = 7),
              axis.line = element_blank(),
              axis.ticks = element_line(linewidth = 0.3),
              panel.border = element_rect(color = "black", fill = NA,
                                          linewidth = 0.4),
              panel.grid.major = element_line(color = "grey92",
                                              linewidth = 0.2),
              legend.key.size = unit(3, "mm"),
              legend.title = element_text(size = 7, face = "bold"),
              legend.text = element_text(size = 6),
              legend.box = "vertical",
              legend.spacing.y = unit(2, "mm"),
              plot.margin = margin(4, 6, 4, 4, "pt"))

      # 7 mm per cell + margin for legend & rotated x-axis labels
      w_mm <- max(150, 7 * n_x + 70)
      h_mm <- max(120, 7 * n_y + 50)
      ggsave(file.path(OUT, "fig6d_bubble_neu_to_mal.pdf"), p6d,
             width = w_mm, height = h_mm, units = "mm", limitsize = FALSE)
      ggsave(file.path(OUT, "fig6d_bubble_neu_to_mal.png"), p6d,
             width = w_mm, height = h_mm, units = "mm", dpi = 300,
             limitsize = FALSE)

      # ── Heatmap version (square tiles, Fig-3/4/5-consistent style) ──
      # Build the full pair × LR grid with NA fill so the panel stays square.
      grid_full <- expand.grid(
        pair = pair_order,
        interaction_name_2 = rev(inter_order),
        stringsAsFactors = FALSE)
      grid_full <- dplyr::left_join(
        grid_full,
        net_show %>% select(pair, interaction_name_2, prob, pval, mlog_p) %>%
          mutate(pair = as.character(pair),
                 interaction_name_2 = as.character(interaction_name_2)),
        by = c("pair", "interaction_name_2"))
      grid_full$pair <- factor(grid_full$pair, levels = pair_order)
      grid_full$interaction_name_2 <- factor(grid_full$interaction_name_2,
                                              levels = rev(inter_order))
      grid_full$sig_lab <- ifelse(is.na(grid_full$pval), "",
                            ifelse(grid_full$pval < 1e-4, "***",
                            ifelse(grid_full$pval < 1e-3, "**",
                            ifelse(grid_full$pval < 5e-2, "*", ""))))

      p6d_hm <- ggplot(grid_full, aes(pair, interaction_name_2, fill = prob)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = sig_lab), size = 2.2, vjust = 0.75,
                  family = my_font, color = "black") +
        scale_fill_gradientn(colors = c("#4575B4","#FDDA82","#D73027"),
                             name = "Commun.\nProb.",
                             na.value = "grey94",
                             guide = guide_colorbar(barwidth = unit(2.5, "mm"),
                                                    barheight = unit(20, "mm"),
                                                    frame.colour = "black",
                                                    frame.linewidth = 0.3,
                                                    ticks.colour = "black")) +
        coord_fixed(ratio = 1) +
        labs(x = NULL, y = NULL,
             caption = "* p<0.05, ** p<0.001, *** p<1e-4") +
        theme_classic(base_size = 9, base_family = my_font) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                         color = "black", size = 7),
              axis.text.y = element_text(color = "black", size = 7),
              axis.line = element_blank(),
              axis.ticks = element_line(linewidth = 0.3),
              panel.border = element_rect(color = "black", fill = NA,
                                          linewidth = 0.4),
              legend.key.size = unit(3, "mm"),
              legend.title = element_text(size = 7, face = "bold"),
              legend.text = element_text(size = 6),
              plot.caption = element_text(size = 6, hjust = 0,
                                          color = "grey30"),
              plot.margin = margin(4, 6, 4, 4, "pt"))

      ggsave(file.path(OUT, "fig6d_heatmap_neu_to_mal.pdf"), p6d_hm,
             width = w_mm, height = h_mm, units = "mm", limitsize = FALSE)
      ggsave(file.path(OUT, "fig6d_heatmap_neu_to_mal.png"), p6d_hm,
             width = w_mm, height = h_mm, units = "mm", dpi = 300,
             limitsize = FALSE)
    }
  }

  # 6E — Neu_Inflammatory → MP1 vs MP3, paired dodged bar (original v2)
  cat("\n[6E] paired bar Neu_Inflammatory → MP1 vs MP3\n")
  for (ext in c("pdf", "png")) {
    f <- file.path(OUT, sprintf("fig6e_inflam_mp1_vs_mp3_dumbbell.%s", ext))
    if (file.exists(f)) file.remove(f)
  }
  if ("Neu_Inflammatory" %in% ident_levels &&
      length(intersect(c("Mal_MP1", "Mal_MP3"), ident_levels)) == 2) {
    net_pE <- subsetCommunication(cc_focus, slot.name = "netP",
                  sources.use = "Neu_Inflammatory",
                  targets.use = c("Mal_MP1", "Mal_MP3"))
    if (nrow(net_pE) > 0) {
      df_e <- net_pE %>%
        group_by(pathway_name, target) %>%
        summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop")
      pw_keep <- df_e %>% group_by(pathway_name) %>%
        summarise(s = sum(prob)) %>% filter(s > 0) %>% pull(pathway_name)
      df_e <- df_e %>% filter(pathway_name %in% pw_keep)
      pw_order <- df_e %>% group_by(pathway_name) %>%
        summarise(s = sum(prob)) %>% arrange(desc(s)) %>% pull(pathway_name)
      df_e$pathway_name <- factor(df_e$pathway_name, levels = pw_order)
      n_pw <- nlevels(df_e$pathway_name)
      # Per-pathway slot ~22mm (2 dodged bars look slim, not blocky); +35mm
      # for y-axis label & legend + plot margin.
      bar_dodge <- 0.55
      p6e <- ggplot(df_e, aes(pathway_name, prob, fill = target)) +
        geom_col(position = position_dodge(width = bar_dodge + 0.05),
                 width = bar_dodge, color = "black", linewidth = 0.25) +
        scale_fill_manual(values = c("Mal_MP1" = "#E64B35",
                                     "Mal_MP3" = "#00A087"), name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
        labs(x = NULL,
             y = "Pathway communication probability\n(Neu_Inflammatory → target)") +
        theme_classic(base_size = 9, base_family = my_font) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(linewidth = 0.4),
              legend.position = "top",
              legend.key.size = unit(3, "mm"),
              plot.margin = margin(4, 6, 4, 4, "pt"))
      w_e <- max(75, 22 * n_pw + 35)
      ggsave(file.path(OUT, "fig6e_inflam_mp1_vs_mp3_bar.pdf"), p6e,
             width = w_e, height = 90, units = "mm")
      ggsave(file.path(OUT, "fig6e_inflam_mp1_vs_mp3_bar.png"), p6e,
             width = w_e, height = 90, units = "mm", dpi = 300)
      cat(sprintf("  6E saved (%d pathways, %d mm wide)\n", n_pw, w_e))
    }
  }

  # 6F — free_y absolute panel + fresh palette (unchanged from v4)
  cat("\n[6F] pathway specificity\n")
  pathway_groups <- list(
    "Neutrophil-leaning"  = c("OSM","IL1","PLAU"),
    "Macrophage-dominant" = c("FN1","SPP1","COLLAGEN","TGFb"),
    "Shared / Other"      = c("VEGF","TNF","CCL")
  )
  all_pw_present <- intersect(unlist(pathway_groups, use.names = FALSE),
                               cc_focus@netP$pathways)
  net_p <- subsetCommunication(cc_focus, slot.name = "netP",
                                signaling = all_pw_present)
  if (nrow(net_p) > 0) {
    sender_class <- function(x) {
      ifelse(grepl("^Neu_", x), "Neutrophil",
      ifelse(grepl("^Macro_", x), "Macrophage",
      ifelse(x %in% c("Fibroblast","Endothelial"), "Stromal", NA_character_)))
    }
    df <- net_p %>% mutate(sender_class = sender_class(source)) %>%
      filter(!is.na(sender_class)) %>%
      group_by(pathway_name, sender_class) %>%
      summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop")
    pw_to_group <- stack(pathway_groups); names(pw_to_group) <- c("pathway_name","group")
    df <- df %>% left_join(pw_to_group, by = "pathway_name")
    df$group <- factor(df$group, levels = names(pathway_groups))
    pw_order <- df %>% group_by(group, pathway_name) %>%
      summarise(s = sum(prob), .groups = "drop") %>%
      arrange(group, desc(s)) %>% pull(pathway_name)
    df$pathway_name <- factor(df$pathway_name, levels = pw_order)
    sender_pal <- c("Neutrophil"="#E64B35","Macrophage"="#4DBBD5","Stromal"="#7E57C2")
    df$sender_class <- factor(df$sender_class,
                               levels = c("Neutrophil","Macrophage","Stromal"))
    df_prop <- df %>% group_by(pathway_name) %>%
      mutate(prop = prob / sum(prob, na.rm = TRUE))

    p_abs <- ggplot(df, aes(pathway_name, prob, fill = sender_class)) +
      geom_col(position = "stack", width = 0.55,
               color = "black", linewidth = 0.25) +
      scale_fill_manual(values = sender_pal, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
      facet_wrap(~ group, scales = "free", nrow = 1) +
      labs(x = NULL, y = "Comm. prob. (absolute)") +
      theme_classic(base_size = 9, base_family = my_font) +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(size = 8),
            axis.line = element_line(linewidth = 0.4),
            strip.text = element_text(face = "bold", size = 9, family = my_font),
            strip.background = element_rect(fill = "grey92", color = NA),
            panel.spacing.x = unit(4, "mm"),
            legend.position = "top",
            legend.key.size = unit(3, "mm"),
            plot.margin = margin(4, 6, 2, 4, "pt"))
    p_prop <- ggplot(df_prop, aes(pathway_name, prop, fill = sender_class)) +
      geom_col(position = "stack", width = 0.55,
               color = "black", linewidth = 0.25) +
      scale_fill_manual(values = sender_pal, name = NULL) +
      facet_wrap(~ group, scales = "free_x", nrow = 1) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(x = NULL, y = "% of pathway total") +
      theme_classic(base_size = 9, base_family = my_font) +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(size = 8),
            axis.line = element_line(linewidth = 0.4),
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.spacing.x = unit(4, "mm"),
            legend.position = "none",
            plot.margin = margin(2, 6, 4, 4, "pt"))
    p6f <- p_abs / p_prop + plot_layout(heights = c(1, 1))
    ggsave(file.path(OUT, "fig6f_pathway_specificity.pdf"), p6f,
           width = 170, height = 130, units = "mm")
    ggsave(file.path(OUT, "fig6f_pathway_specificity.png"), p6f,
           width = 170, height = 130, units = "mm", dpi = 300)
    write.csv(df_prop, file.path(OUT, "fig6f_pathway_specificity.csv"),
              row.names = FALSE)
  }
}

cat("\n=== DONE ===\n")
print(list.files(OUT))
