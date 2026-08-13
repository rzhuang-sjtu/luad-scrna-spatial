## Step 7b: MISTy multi-view RF on PROGENy targets, parallelized across sections.
suppressPackageStartupMessages({
  library(mistyR); library(dplyr); library(tidyr); library(readr); library(future)
})

DATA_DIR <- "${DATA_ROOT}/ST/results/step07_misty/data"
OUT_DIR  <- "${DATA_ROOT}/ST/results/step07_misty"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- gsub("_intra.csv$", "", list.files(DATA_DIR, pattern = "_intra.csv$"))
cat("[info] sections:", length(samples), "->", paste(samples, collapse=" "), "\n")

# Use future for parallel section runs
n_workers <- min(length(samples), parallel::detectCores() - 2L)
cat("[info] workers:", n_workers, "\n")
plan(multisession, workers = n_workers)

# We pick PROGENy pathway targets relevant to Fig 7
TARGETS <- c("progeny_NFkB", "progeny_JAK.STAT", "progeny_TNFa", "progeny_TGFb",
             "progeny_Hypoxia", "progeny_MAPK", "progeny_EGFR", "progeny_p53", "progeny_VEGF")

run_one <- function(s) {
  intra   <- read.csv(file.path(DATA_DIR, paste0(s, "_intra.csv")), row.names = 1, check.names = FALSE)
  progeny <- read.csv(file.path(DATA_DIR, paste0(s, "_progeny.csv")), row.names = 1, check.names = FALSE)
  coords  <- read.csv(file.path(DATA_DIR, paste0(s, "_coords.csv")), row.names = 1)
  # Some R column names get . substitution; column "JAK-STAT" becomes "JAK.STAT" — fine
  # Ensure progeny has the named columns
  cols_have <- intersect(TARGETS, colnames(progeny))
  if (length(cols_have) == 0) {
    return(list(sample = s, status = "no_targets",
                cols = paste(colnames(progeny), collapse = ",")))
  }
  # Combine cell-type abundances + pathway scores into a single matrix.
  # MISTy treats every column as a potential predictor and any subset as targets.
  full <- cbind(intra, progeny[, cols_have, drop = FALSE])
  full[is.na(full)] <- 0
  pathway_cols <- cols_have

  # Build MISTy views: intra + juxta + para (same predictor matrix, different spatial scopes)
  views <- create_initial_view(full) |>
    add_juxtaview(positions = coords, neighbor.thr = 130) |>           # immediate neighbors
    add_paraview(positions = coords, l = 220, family = "gaussian")    # broader gaussian

  result_dir <- file.path(OUT_DIR, "misty_runs", s)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  # Run MISTy: targets = pathway columns (cell types serve only as predictors)
  run_misty(views, results.folder = result_dir, model.function = random_forest_model,
            num.trees = 100, importance.cutoff = 0,
            target.subset = pathway_cols)

  # Collect results
  res <- collect_results(result_dir)

  list(
    sample = s,
    status = "ok",
    importance = res$importances.aggregated,
    r2 = res$improvements.stats,
    contrib = res$contributions.stats
  )
}

t0 <- Sys.time()
fut_results <- future.apply::future_lapply(samples, run_one, future.seed = TRUE)
t1 <- Sys.time()
cat("[info] all sections done in", as.numeric(difftime(t1, t0, units="mins")), "min\n")

# Aggregate
all_imp <- do.call(rbind, lapply(fut_results, function(x) {
  if (x$status != "ok") return(NULL)
  if (is.null(x$importance)) return(NULL)
  cbind(sample = x$sample, x$importance)
}))
write.csv(all_imp, file.path(OUT_DIR, "per_section_importance.csv"), row.names = FALSE)

all_r2 <- do.call(rbind, lapply(fut_results, function(x) {
  if (x$status != "ok") return(NULL)
  cbind(sample = x$sample, x$r2)
}))
write.csv(all_r2, file.path(OUT_DIR, "per_section_r2.csv"), row.names = FALSE)

# Mean importance per (view, target, predictor)
agg <- all_imp |>
  group_by(view, Target, Predictor) |>
  summarize(mean_importance = mean(Importance, na.rm = TRUE), .groups = "drop")
write.csv(agg, file.path(OUT_DIR, "aggregated_importance.csv"), row.names = FALSE)
cat("[info] saved per_section_importance.csv, per_section_r2.csv, aggregated_importance.csv\n")

# Final misty-style heatmap from a virtual aggregate folder
agg_dir <- file.path(OUT_DIR, "misty_runs")
all_subdirs <- list.dirs(agg_dir, recursive = FALSE)
cat("[info] misty subdirs:", length(all_subdirs), "\n")
agg_res <- collect_results(all_subdirs)

png(file.path(OUT_DIR, "heatmap_misty_importances.png"), width = 1400, height = 1000, res = 130)
print(plot_interaction_heatmap(agg_res, view = "intra", clean = TRUE) %>%
        suppressWarnings())
dev.off()
png(file.path(OUT_DIR, "heatmap_misty_juxta.png"), width = 1400, height = 1000, res = 130)
print(plot_interaction_heatmap(agg_res, view = "juxta.130", clean = TRUE) %>%
        suppressWarnings())
dev.off()
png(file.path(OUT_DIR, "heatmap_misty_para.png"), width = 1400, height = 1000, res = 130)
print(plot_interaction_heatmap(agg_res, view = "para.220", clean = TRUE) %>%
        suppressWarnings())
dev.off()
cat("[done]\n")
