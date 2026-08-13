#!/usr/bin/env Rscript
# Step 4c: Run CopyKAT for a single patient
# Usage: Rscript data_prep/copykat/04c_run_copykat.R <patient_dir> <n_cores>

# user library (avoid site-library permission issues)
user_lib <- file.path(Sys.getenv("HOME"), "R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
    library(copykat)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 04c_run_copykat.R <patient_dir> <n_cores>")
patient_dir <- normalizePath(args[1], mustWork = TRUE)
n_cores <- as.integer(args[2])
sam_name <- basename(patient_dir)

t0 <- Sys.time()
cat(sprintf("[%s] START %s (n_cores=%d)\n", format(t0, "%H:%M:%S"), sam_name, n_cores))

# ---- Read inputs ----
counts <- readMM(file.path(patient_dir, "counts.mtx"))     # genes × cells
barcodes <- readLines(file.path(patient_dir, "barcodes.tsv"))
genes    <- readLines(file.path(patient_dir, "genes.tsv"))
ref_cells <- readLines(file.path(patient_dir, "ref_cells.txt"))
ref_cells <- ref_cells[nzchar(ref_cells)]

stopifnot(nrow(counts) == length(genes), ncol(counts) == length(barcodes))

rawmat <- as.matrix(counts)
rownames(rawmat) <- genes
colnames(rawmat) <- barcodes
rm(counts); gc(verbose = FALSE)

cat(sprintf("  matrix: %d genes × %d cells,  ref_cells: %d\n",
            nrow(rawmat), ncol(rawmat), length(ref_cells)))

# ---- cd to patient dir (CopyKAT writes temp files to cwd) ----
old_wd <- getwd()
setwd(patient_dir)
on.exit(setwd(old_wd), add = TRUE)

# ---- Run CopyKAT ----
common_args <- list(
    rawmat       = rawmat,
    id.type      = "S",
    ngene.chr    = 5,
    win.size     = 25,
    KS.cut       = 0.1,
    sam.name     = sam_name,
    distance     = "euclidean",
    output.seg   = FALSE,
    plot.genes   = FALSE,
    genome       = "hg20",
    n.cores      = n_cores
)
if (length(ref_cells) > 0) {
    common_args$norm.cell.names <- ref_cells
    cat(sprintf("  mode: WITH reference (%d cells)\n", length(ref_cells)))
} else {
    cat("  mode: REF-FREE\n")
}

result <- tryCatch(
    do.call(copykat, common_args),
    error = function(e) {
        msg <- conditionMessage(e)
        cat(sprintf("copykat() failed: %s\n", msg))
        writeLines(msg, file.path(patient_dir, "copykat_error.txt"))
        NULL
    }
)

if (is.null(result)) {
    cat(sprintf("[%s] FAILED %s\n", format(Sys.time(), "%H:%M:%S"), sam_name))
    quit(status = 2)
}

# ---- Save slim results ----
pred <- result$prediction
stopifnot("cell.names" %in% colnames(pred), "copykat.pred" %in% colnames(pred))
write.csv(pred, file.path(patient_dir, "copykat_prediction.csv"), row.names = FALSE)

if (!is.null(result$CNAmat)) {
    saveRDS(result$CNAmat, file.path(patient_dir, "copykat_CNA.rds"))
}

# Integrity check
n_pred <- nrow(pred)
n_in <- ncol(rawmat)
cat(sprintf("  prediction rows: %d  (input cells: %d)\n", n_pred, n_in))
if (n_pred < n_in * 0.9) {
    cat(sprintf("Predicted cells far fewer than input (%.1f%%)\n", 100 * n_pred / n_in))
}

# Label distribution
pred_counts <- table(pred$copykat.pred)
cat("  label distribution:\n")
for (lab in names(pred_counts)) {
    cat(sprintf("    %s: %d\n", lab, pred_counts[[lab]]))
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[%s] DONE  %s  elapsed=%.1f min\n",
            format(Sys.time(), "%H:%M:%S"), sam_name, elapsed))
