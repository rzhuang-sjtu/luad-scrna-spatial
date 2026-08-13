#!/usr/bin/env Rscript
# Shared CellChat build helper — sourced by cellchat_{all,normal,tumor,met}.R
# Each per-tissue script sets GROUP + H5AD + RDS_OUT, then sources this file.

suppressPackageStartupMessages({
  library(CellChat)
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)
# Each per-tissue process gets its own R session — keep its internal future
# plan sequential so 3 parallel processes don't fight for cores.
future::plan("sequential")
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

build_cellchat <- function(h5ad_path, label_col = "cellchat_label",
                            min_cells_per_group = 10) {
  cat("  reading", h5ad_path, "\n")
  sce <- readH5AD(h5ad_path, raw = FALSE, use_hdf5 = FALSE)
  expr <- assay(sce, "X")
  if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "CsparseMatrix")
  meta <- as.data.frame(colData(sce))
  meta$labels <- as.character(meta[[label_col]])
  tab <- table(meta$labels)
  keep_lbl <- names(tab)[tab >= min_cells_per_group]
  keep_idx <- meta$labels %in% keep_lbl
  expr <- expr[, keep_idx]
  meta <- meta[keep_idx, , drop = FALSE]
  cat("  after group-filter:", ncol(expr), "cells across",
      length(unique(meta$labels)), "groups\n")
  print(sort(table(meta$labels)))
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

run_one <- function(group_name, h5ad_path, rds_out, force = FALSE) {
  if (file.exists(rds_out) && !force) {
    cat(sprintf("[%s] cached %s — skipping (set FORCE=1 to overwrite)\n",
                group_name, rds_out))
    return(invisible(NULL))
  }
  cat(sprintf("[%s] building from %s\n", group_name, h5ad_path))
  t1 <- Sys.time()
  cc <- build_cellchat(h5ad_path)
  dir.create(dirname(rds_out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cc, rds_out)
  cat(sprintf("[%s] done in %.1f min → %s (%.0f MB)\n", group_name,
              as.numeric(difftime(Sys.time(), t1, units = "mins")),
              rds_out, file.info(rds_out)$size / 1e6))
}
