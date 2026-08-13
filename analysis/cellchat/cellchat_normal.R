#!/usr/bin/env Rscript
# Build CellChat object for Normal tissue group (independent process).
source("${PROJECT_ROOT}/analysis/cellchat/_cellchat_build_helper.R")
DAT <- "${PROJECT_ROOT}/data/processed"
RDS <- file.path(DAT, "cellchat_rds")
FORCE <- nzchar(Sys.getenv("FORCE"))
run_one("Normal",
        file.path(DAT, "cellchat_input_Normal.h5ad"),
        file.path(RDS, "cellchat_Normal.rds"),
        force = FORCE)
