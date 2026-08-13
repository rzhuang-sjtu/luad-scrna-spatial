suppressPackageStartupMessages(library(CellChat))
cc <- readRDS("${PROJECT_ROOT}/data/processed/cellchat_rds/cellchat_all.rds")
pw <- cc@netP$pathways
cat("class:", class(pw), " length:", length(pw), "\n")
cat("sorted:\n", paste(sort(pw), collapse = ", "), "\n\n")
target <- c("IL1","OSM","PLAU","TGFb","FN1","COLLAGEN","SPP1","VEGF","TNF","CXCL","CCL")
cat("requested vs present:\n")
for (p in target) cat(sprintf("  %-10s : %s\n", p, p %in% pw))
