#!/usr/bin/env Rscript
# Recompute UMAP using the 4-D MP score matrix as input.
# Output: malignant_mp_umap_metadata.csv.gz with new MP_UMAP1/2 columns.
# Use case: Fig 2C MP-space variant where MPs are guaranteed to separate.
suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
})

setwd(if (.Platform$OS.type == "windows")
  "${WORK_ROOT}/luad_figures/fig2" else
  "${WORK_ROOT}/luad_figures/fig2")

cat("[load] malignant_umap_metadata.csv.gz\n")
meta <- fread("malignant_umap_metadata.csv.gz")
cat(sprintf("  %d cells\n", nrow(meta)))

# 4 MP score columns (drop MP5 user uses 1-4 only
score_cols <- c("MP1_score", "MP2_score", "MP3_score", "MP4_score")
stopifnot(all(score_cols %in% names(meta)))
M <- as.matrix(meta[, ..score_cols])

cat("[run] uwot::umap on 4-D MP score matrix\n")
set.seed(42)
um <- uwot::umap(M,
                 n_neighbors = 30,
                 min_dist    = 0.3,
                 n_components = 2,
                 metric = "cosine",
                 verbose = TRUE)

meta$MP_UMAP1 <- um[, 1]
meta$MP_UMAP2 <- um[, 2]

out_path <- "malignant_mp_umap_metadata.csv.gz"
fwrite(meta, out_path)
cat(sprintf("[save] %s  (%d rows)\n", out_path, nrow(meta)))
