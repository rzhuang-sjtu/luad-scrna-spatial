#!/bin/bash
# Step 4a: Install R + CopyKAT environment
# Run: bash scripts/04a_setup_R.sh
set -euo pipefail

# 1. Check R
# This script does not install system packages for you. Installing R and its
# system libraries needs root, and a published script should not change
# someone else's machine without being asked.
if ! command -v Rscript &> /dev/null; then
    echo "Rscript not found on PATH." >&2
    echo "Install R and the system libraries it needs, then re-run. On Debian/Ubuntu:" >&2
    echo "  apt-get install r-base r-base-dev \\" >&2
    echo "    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \\" >&2
    echo "    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \\" >&2
    echo "    libtiff5-dev libjpeg-dev" >&2
    exit 1
fi
echo "R found:"
Rscript --version

# 2. Prepare user library (avoid site-library permission issues)
USER_LIB="$HOME/R/library"
mkdir -p "$USER_LIB"
export R_LIBS_USER="$USER_LIB"
echo "R_LIBS_USER=$USER_LIB"

# 3. Install CopyKAT and dependencies in R
echo ""
echo "Check / install R packages..."
Rscript - <<'RSCRIPT'
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 8)   # parallel compile
user_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(user_lib, .libPaths()))
cat(".libPaths:\n"); print(.libPaths())

# Core deps (lightweight remotes instead of devtools)
deps <- c("remotes", "Matrix", "Rtsne", "Rcpp")
for (pkg in deps) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, lib = user_lib)
    } else {
        cat("Already installed", pkg, "\n")
    }
}

# CopyKAT
if (!requireNamespace("copykat", quietly = TRUE)) {
    cat("Installing copykat from GitHub (navinlabcode/copykat)...\n")
    remotes::install_github("navinlabcode/copykat",
                            upgrade = "never", lib = user_lib, dependencies = TRUE)
} else {
    cat("copykat already installed\n")
}

# Verify
library(copykat)
library(Matrix)
cat("\n=== Versions ===\n")
cat("copykat:", as.character(packageVersion("copykat")), "\n")
cat("Matrix :", as.character(packageVersion("Matrix")),  "\n")
cat("R      :", R.version.string, "\n")
cat("\n R + CopyKAT ready\n")
RSCRIPT
