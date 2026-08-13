#!/usr/bin/env Rscript
# Build CellChat object for the pooled "all" sample (independent process).
# This is the heaviest of the four (50k cells) — needed for panels 6D/6E/6F.
source("${PROJECT_ROOT}/analysis/cellchat/_cellchat_build_helper.R")
DAT <- "${PROJECT_ROOT}/data/processed"
RDS <- file.path(DAT, "cellchat_rds")
FORCE <- nzchar(Sys.getenv("FORCE"))
run_one("all",
        file.path(DAT, "cellchat_input_all.h5ad"),
        file.path(RDS, "cellchat_all.rds"),
        force = FORCE)
