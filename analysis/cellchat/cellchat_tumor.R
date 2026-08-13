#!/usr/bin/env Rscript
# Build CellChat object for Tumor tissue group (independent process).
source("${PROJECT_ROOT}/analysis/cellchat/_cellchat_build_helper.R")
DAT <- "${PROJECT_ROOT}/data/processed"
RDS <- file.path(DAT, "cellchat_rds")
FORCE <- nzchar(Sys.getenv("FORCE"))
run_one("Tumor",
        file.path(DAT, "cellchat_input_Tumor.h5ad"),
        file.path(RDS, "cellchat_Tumor.rds"),
        force = FORCE)
