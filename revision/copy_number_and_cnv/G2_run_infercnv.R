#!/usr/bin/env Rscript
# G2 — run inferCNV per patient to independently recheck CopyKAT malignant calls (independent CNV cross-check).
#
# Usage: Rscript G2_run_infercnv.R <start_index> <end_index>
# Invoked by G3_launch_infercnv.sh as 12 parallel slices.
#
# Disk policy (data volume only 50G; full outputs for 89 patients would fill it):
#   After each patient finishes → immediately extract per-cell summary (a few MB) → delete all raw outputs for that patient.
#   Kept summaries suffice for agreement with CopyKAT; raw heatmaps are not for the paper and can be regenerated for a single patient if needed.
#
# Per-cell continuous summaries (threshold-free; compare at multiple thresholds locally to avoid cherry-picking):
#   cnv_var       variance of the cell CNV profile
#   cnv_mean_abs  mean of |x-1|
#   cnv_rms       sqrt(mean((x-1)^2))
#   cnv_p90_abs   90th percentile of |x-1| (most deviant decile of genes only)
#   cnv_frac_dev  fraction of genes off the reference plateau; see below
#   frac_altered  fraction of the genome called non-diploid by HMM (state != 3); subclone-level
#
# Do not use median absolute deviation (cnv_mad): on P38 (101 cells) it took only one value.
# clear_noise_via_ref_quantiles collapses genes inside the reference quantile band to a cell-independent
# plateau value; once more than half the genes are flattened the median equals the plateau and carries no cell-level information.
# cnv_frac_dev is designed for this: define a per-gene plateau from reference cells, then count genes that escape it per cell.
#
# Also keep metrics on reference cells themselves (T/NK + B) as an internal negative control —
# shows what the same scale reports for non-epithelial cells.
suppressMessages({library(infercnv); library(Matrix); library(matrixStats)})
options(scipen = 100)   # without this under subclusters mode, hclust may error (infercnv warning)

# Bayesian post-hoc filtering (BayesNet MCMC) is the main cost and only affects frac_altered;
# cnv_mad / cnv_var come from the denoised residual matrix and are finished before this step.
# Setting to 0 reports raw HMM calls without post-filtering — more conservative for independent validation.
BAYES <- as.numeric(Sys.getenv("BAYES_MAX_PNORMAL", "0"))

# Two call modes:
#   Rscript G2_run_infercnv.R --single <patient>     dispatched one-by-one by G3 (with memory guardrails)
#   Rscript G2_run_infercnv.R <start_index> <end_index>  manual batch (not recommended; no guardrails)
args <- commandArgs(trailingOnly = TRUE)
SINGLE <- NULL
if (length(args) >= 2 && args[1] == "--single") {
  SINGLE <- args[2]
} else {
  i0 <- as.integer(args[1]); i1 <- as.integer(args[2])
}
# inferCNV was run on a rented GPU node; ROOT there was a scratch path.
# Set INFERCNV_ROOT to point at the scratch directory on your own machine.
ROOT <- Sys.getenv("INFERCNV_ROOT",
                   unset = file.path(Sys.getenv("WORK_ROOT", unset = "."),
                                     "infercnv"))
IN   <- file.path(ROOT, "input")
OUT  <- file.path(ROOT, "output")
SUM  <- file.path(ROOT, "summary")     # keep only this directory
GO   <- file.path(ROOT, "gene_order.tsv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(SUM, showWarnings = FALSE, recursive = TRUE)

if (!is.null(SINGLE)) {
  pats <- SINGLE
  cat(sprintf("[%s] single-patient mode: %s\n", format(Sys.time(), "%H:%M:%S"), SINGLE))
} else {
  pats <- sort(list.dirs(IN, full.names = FALSE, recursive = FALSE))
  pats <- pats[i0:min(i1, length(pats))]
  cat(sprintf("[%s] handling %d patients: %s ... %s\n", format(Sys.time(), "%H:%M:%S"),
              length(pats), pats[1], tail(pats, 1)))
}

summarise_and_clean <- function(p, od) {
  # 1) Denoised observation matrix from run.final.infercnv_obj.
  #    Do not use infercnv.observations.txt — from infercnv::run source, write_expr_matrix
  #    is only passed to plot_cnv, and plot_cnv is entirely inside if(!no_plot),
  #    so when no_plot=TRUE that observations txt is never written
  f_rds <- file.path(od, "run.final.infercnv_obj")
  if (!file.exists(f_rds)) return(FALSE)
  o <- readRDS(f_rds)
  obs <- unlist(o@observation_grouped_cell_indices, use.names = FALSE)
  ref <- unlist(o@reference_grouped_cell_indices, use.names = FALSE)
  E <- o@expr.data

  # Reference plateau: per-gene median over reference cells. Denoising collapses in-band genes to a cell-independent value;
  # reference cells are near-diploid by construction, so they give the most stable plateau.
  plateau <- matrixStats::rowMedians(E[, ref, drop = FALSE], na.rm = TRUE)

  cell_metrics <- function(idx, grp) {
    M  <- E[, idx, drop = FALSE]
    D  <- abs(M - 1)
    n  <- nrow(M)
    v  <- matrixStats::colVars(M, na.rm = TRUE)
    mu <- matrixStats::colMeans2(M, na.rm = TRUE)
    data.frame(
      cell         = colnames(M),
      group        = grp,
      cnv_var      = v,
      cnv_mean_abs = matrixStats::colMeans2(D, na.rm = TRUE),
      # mean((x-1)^2) = var*(n-1)/n + (mean-1)^2; avoids allocating another full matrix
      cnv_rms      = sqrt(v * (n - 1) / n + (mu - 1)^2),
      cnv_p90_abs  = as.numeric(matrixStats::colQuantiles(D, probs = 0.9, na.rm = TRUE)),
      cnv_frac_dev = colMeans(abs(M - plateau) > 1e-9, na.rm = TRUE),
      row.names    = NULL)
  }
  res <- rbind(cell_metrics(obs, "epithelial"), cell_metrics(ref, "reference"))
  rm(E, o); invisible(gc())

  # 2) HMM-called non-diploid fraction.
  #    pred_cnv_genes.dat lists only genes in non-neutral segments (states 1,2,4,5,6; no state 3),
  #    so per-subclone frac_altered = unique listed genes / total genes_used.
  #    Subclone-level quantity; all cells in a subclone share the value — as expected under i6-subclusters.
  f_grp <- list.files(od, pattern = "cell_groupings$", full.names = TRUE)
  f_gen <- list.files(od, pattern = "pred_cnv_genes\\.dat$", full.names = TRUE)
  f_use <- list.files(od, pattern = "genes_used\\.dat$", full.names = TRUE)
  res$frac_altered <- NA_real_
  hmm_ok <- "none"
  if (length(f_grp) && length(f_gen) && length(f_use)) {
    grp <- read.table(f_grp[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    gen <- read.table(f_gen[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    ntot <- nrow(read.table(f_use[1], header = TRUE, sep = "\t", row.names = 1))
    gen <- gen[gen$state != 3, , drop = FALSE]
    # Subclones with no altered segments get no tapply entry; direct lookup becomes NA.
    # True value is 0 (genome-wide diploid) and must be filled — otherwise the cleanest reference subclones
    # are dropped and control readouts are inflated. On P38, 60 of 70 reference cells hit this.
    all_grp <- unique(grp$cell_group_name)
    fa <- setNames(rep(0, length(all_grp)), all_grp)
    tab <- tapply(gen$gene, gen$cell_group_name,
                  function(g) length(unique(g)) / ntot)
    fa[names(tab)] <- as.numeric(tab)
    res$frac_altered <- as.numeric(fa[grp$cell_group_name[match(res$cell, grp$cell)]])
    hmm_ok <- sprintf("%s (%d genes, %d subclones)", basename(f_gen[1]), ntot, length(all_grp))
    # Archive the three small HMM files compressed so frac_altered can be recomputed locally later
    # without re-running multi-hour inferCNV for a definition change. ~100 KB per patient after compression.
    hd <- file.path(SUM, paste0(p, "_hmm"))
    dir.create(hd, showWarnings = FALSE)
    for (f in c(f_grp[1], f_gen[1], f_use[1])) {
      file.copy(f, file.path(hd, basename(f)), overwrite = TRUE)
    }
    system2("gzip", c("-f", shQuote(list.files(hd, full.names = TRUE))))
  }
  res$patient <- p
  write.csv(res, file.path(SUM, paste0(p, "_cells.csv")), row.names = FALSE)
  # Leave a file listing before deletion: if extraction fails, check whether names changed or files were never written,
  # without re-running the whole patient for diagnosis.
  writeLines(c(paste0("# HMM source:", hmm_ok),
               paste0("# frac_altered non-NA:", sum(!is.na(res$frac_altered)),
                      " / ", nrow(res)),
               list.files(od, recursive = TRUE)),
             file.path(SUM, paste0(p, "_files.txt")))
  unlink(od, recursive = TRUE)      # clear raw outputs immediately after the summary is written
  TRUE
}

for (p in pats) {
  if (file.exists(file.path(SUM, paste0(p, "_cells.csv")))) {
    cat(sprintf("skip %s (summary already present)\n", p)); next
  }
  od <- file.path(OUT, p); t0 <- Sys.time()
  ok <- tryCatch({
    m  <- readMM(file.path(IN, p, "counts.mtx"))
    g  <- readLines(file.path(IN, p, "genes.tsv"))
    bc <- readLines(file.path(IN, p, "cells.tsv"))
    m  <- as(m, "CsparseMatrix"); rownames(m) <- g; colnames(m) <- bc
    dir.create(od, showWarnings = FALSE, recursive = TRUE)
    obj <- CreateInfercnvObject(raw_counts_matrix = m,
                                annotations_file = file.path(IN, p, "annot.tsv"),
                                delim = "\t", gene_order_file = GO,
                                ref_group_names = c("reference"))
    obj <- infercnv::run(obj, cutoff = 0.1, out_dir = od,
                         cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE,
                         num_threads = 2, no_plot = TRUE,
                         BayesMaxPNormal = BAYES,
                         # drop intermediate rds (disk), but keep the final object —
                         # the denoised matrix can only be taken from it; see comments in summarise_and_clean.
                         save_rds = FALSE, save_final_rds = TRUE)
    summarise_and_clean(p, od)
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(SUM, paste0("FAILED_", p, ".txt")))
    unlink(od, recursive = TRUE)
    FALSE
  })
  cat(sprintf("[%s] %s %s (%.1f min)\n", format(Sys.time(), "%H:%M:%S"), p,
              ifelse(isTRUE(ok), "done", "失败"),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] process finished\n", format(Sys.time(), "%H:%M:%S")))
