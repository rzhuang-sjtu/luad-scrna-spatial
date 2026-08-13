#!/usr/bin/env Rscript
# Build CellChat object for Metastasis tissue group (independent process).
source("${PROJECT_ROOT}/analysis/cellchat/_cellchat_build_helper.R")
DAT <- "${PROJECT_ROOT}/data/processed"
RDS <- file.path(DAT, "cellchat_rds")
FORCE <- nzchar(Sys.getenv("FORCE"))
run_one("Metastasis",
        file.path(DAT, "cellchat_input_Metastasis.h5ad"),
        file.path(RDS, "cellchat_Metastasis.rds"),
        force = FORCE)
