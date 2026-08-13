#!/usr/bin/env Rscript
# D2 — run decontX per sample (celda package, Yang et al. 2020 Genome Biology).
#
# Model: each cell's counts are a mixture of the true cell-type expression profile and a sample-shared background profile;
# variational Bayes estimates per-cell contamination and returns a decontaminated count matrix.
# Requires all cell types from the same sample as cluster prior (z); a single cell type cannot identify the background.
#
# Output (per sample):
#   decontaminated.mtx  decontaminated counts (genes x cells, neutrophil columns only, to save space)
#   contamination.tsv   per-cell contamination fraction (all cells)
suppressMessages({library(celda); library(Matrix); library(SingleCellExperiment)})

IN  <- "${PROJECT_ROOT}/results/decontx/input"
OUT <- "${PROJECT_ROOT}/results/decontx/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
samples <- list.dirs(IN, full.names = FALSE, recursive = FALSE)
cat("n samples:", length(samples), "\n")

done <- 0; failed <- character()
t0 <- Sys.time()
for (s in samples) {
  outd <- file.path(OUT, s)
  if (file.exists(file.path(outd, "contamination.tsv"))) { done <- done + 1; next }
  res <- tryCatch({
    m  <- readMM(file.path(IN, s, "matrix.mtx"))
    g  <- readLines(file.path(IN, s, "genes.tsv"))
    bc <- readLines(file.path(IN, s, "barcodes.tsv"))
    lb <- readLines(file.path(IN, s, "labels.tsv"))
    m  <- as(m, "CsparseMatrix"); rownames(m) <- g; colnames(m) <- bc
    # decontX cannot estimate background with fewer than 2 cell types
    if (length(unique(lb)) < 2) stop("cell types < 2")
    sce <- SingleCellExperiment(list(counts = m))
    sce <- decontX(sce, z = lb, verbose = FALSE)
    dir.create(outd, showWarnings = FALSE, recursive = TRUE)
    cont <- data.frame(barcode = bc, label = lb,
                       contamination = colData(sce)$decontX_contamination)
    write.table(cont, file.path(outd, "contamination.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # Save decontaminated counts for neutrophil columns only (other types unused; save disk)
    keep <- which(lb == "Neutrophil")
    if (length(keep) > 0) {
      dec <- round(decontXcounts(sce)[, keep, drop = FALSE])
      writeMM(as(dec, "CsparseMatrix"), file.path(outd, "decontaminated.mtx"))
      writeLines(bc[keep], file.path(outd, "neu_barcodes.tsv"))
    }
    TRUE
  }, error = function(e) { failed <<- c(failed, paste0(s, ": ", conditionMessage(e))); FALSE })
  done <- done + 1
  if (done %% 10 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%3d/%d] elapsed %.1f min, ~%.1f min remaining\n",
                done, length(samples), el, el / done * (length(samples) - done)))
    flush.console()
  }
}
cat("\nDone", done, "个样本；失败", length(failed), "个\n")
if (length(failed)) writeLines(failed, file.path(OUT, "failed.txt"))
cat("Output directory:", OUT, "\n")
