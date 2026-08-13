## Step 9c: MISTy on Okamura LUAD No.1-5 cohort.
## Same logic as step07b_misty.R but pointing at step09 paths.
suppressPackageStartupMessages({
  library(mistyR); library(dplyr); library(tidyr); library(readr); library(future); library(future.apply)
})

ROOT <- "${DATA_ROOT}/ST/results/step09_okamura_validation"
DATA_DIR <- file.path(ROOT, "misty_data")
OUT_DIR  <- file.path(ROOT, "misty")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- gsub("_intra.csv$", "", list.files(DATA_DIR, pattern = "_intra.csv$"))
cat("[info] sections:", length(samples), ":", paste(samples, collapse=" "), "\n")

n_workers <- min(length(samples), parallel::detectCores() - 2L)
cat("[info] workers:", n_workers, "\n")
plan(multisession, workers = n_workers)

TARGETS <- c("progeny_NFkB","progeny_JAK.STAT","progeny_TNFa","progeny_TGFb",
             "progeny_Hypoxia","progeny_MAPK","progeny_EGFR","progeny_p53","progeny_VEGF")

run_one <- function(s) {
  intra   <- read.csv(file.path(DATA_DIR, paste0(s, "_intra.csv")), row.names = 1, check.names = FALSE)
  progeny <- read.csv(file.path(DATA_DIR, paste0(s, "_progeny.csv")), row.names = 1, check.names = FALSE)
  coords  <- read.csv(file.path(DATA_DIR, paste0(s, "_coords.csv")), row.names = 1)
  cols_have <- intersect(TARGETS, colnames(progeny))
  if (length(cols_have) == 0) return(list(sample=s, status="no_targets"))
  full <- cbind(intra, progeny[, cols_have, drop = FALSE])
  full[is.na(full)] <- 0
  views <- create_initial_view(full) |>
    add_juxtaview(positions = coords, neighbor.thr = 130) |>
    add_paraview(positions = coords, l = 220, family = "gaussian")
  result_dir <- file.path(OUT_DIR, "misty_runs", s)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  run_misty(views, results.folder = result_dir, model.function = random_forest_model,
            num.trees = 100, importance.cutoff = 0,
            target.subset = cols_have)
  res <- collect_results(result_dir)
  list(sample = s, status = "ok",
       importance = res$importances.aggregated,
       r2 = res$improvements.stats)
}

t0 <- Sys.time()
fut_results <- future.apply::future_lapply(samples, run_one, future.seed = TRUE)
t1 <- Sys.time()
cat("[info] all sections done in", as.numeric(difftime(t1, t0, units="mins")), "min\n")

all_imp <- do.call(rbind, lapply(fut_results, function(x) {
  if (x$status != "ok" || is.null(x$importance)) return(NULL)
  cbind(sample = x$sample, x$importance)
}))
write.csv(all_imp, file.path(OUT_DIR, "per_section_importance.csv"), row.names = FALSE)
all_r2 <- do.call(rbind, lapply(fut_results, function(x) {
  if (x$status != "ok") return(NULL)
  cbind(sample = x$sample, x$r2)
}))
write.csv(all_r2, file.path(OUT_DIR, "per_section_r2.csv"), row.names = FALSE)
agg <- all_imp |>
  group_by(view, Target, Predictor) |>
  summarize(mean_importance = mean(Importance, na.rm = TRUE), .groups = "drop")
write.csv(agg, file.path(OUT_DIR, "aggregated_importance.csv"), row.names = FALSE)

# Aggregate plots from misty results folders
all_subdirs <- list.dirs(file.path(OUT_DIR, "misty_runs"), recursive = FALSE)
agg_res <- collect_results(all_subdirs)
for (vw in c("intra","juxta.130","para.220")) {
  png(file.path(OUT_DIR, paste0("heatmap_", gsub("\\.","_", vw), ".png")), width=1400, height=1000, res=130)
  print(suppressWarnings(plot_interaction_heatmap(agg_res, view = vw, clean = TRUE)))
  dev.off()
}
cat("[done]\n")
