# Installation

This project uses **two conda environments** for Python analysis plus **system R 4.3.3**
for figure rendering and a handful of analytic R scripts (CellChat, MisTy, NMF).

## 1. Conda envs

### `scst` — main scanpy/scVI environment (Figs 1–7 analysis, Fig 8 downstream)

```bash
mamba env create -f scst.yml -n scst
mamba activate scst

# The *-pip-freeze.txt files are a provenance snapshot of the environment
# actually used, not an installable requirements file: some entries point at
# local build trees. Create the environment from the .yml above, then consult
# the freeze for exact versions if a result does not reproduce.
```

Key Python packages: `scanpy`, `scvi-tools`, `anndata`, `squidpy`, `decoupler`,
`liana`, `cell2location`, `commot`, `rapids-singlecell`, `cuML`, `numba`, `pytables`.

### `geneformer` — Fig 8 in-silico perturbation only

```bash
mamba env create -f geneformer.yml -n geneformer
mamba activate geneformer
pip install -r geneformer-pip-freeze.txt
```

This is intentionally separated from `scst` because Geneformer V2-104M pins
`transformers`/`torch` versions that conflict with `scvi-tools`.

## 2. R environment

R 4.3.3 (Linux x86_64). A curated list of the packages required by scripts
in `plotting/` and `analysis/cellchat/` is in `R_packages.tsv` with versions
captured from the production environment.

Install with:

```r
# 1. CRAN packages
install.packages(c(
  "Seurat", "SeuratObject", "survival", "survminer",
  "ggplot2", "ggpubr", "ggrepel", "ggalluvial", "patchwork",
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "VennDiagram", "pheatmap", "viridis", "RColorBrewer",
  "future", "future.apply",
  "ggridges", "ggsci", "ggsignif", "R.utils", "jsonlite", "yaml",
  "data.table", "eulerr", "ggdist", "circlize",
  "scales", "showtext", "sysfonts", "ggbeeswarm", "ggnewscale",
  "ggrastr", "ggtext", "ragg", "png", "uwot", "enrichR",
  "matrixStats", "dendextend", "BiocManager"
))

# 2. Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "ComplexHeatmap", "SingleCellExperiment", "SummarizedExperiment",
  "edgeR", "limma", "DESeq2", "GSVA", "GSEABase",
  "clusterProfiler", "enrichplot", "msigdbr", "DOSE",
  "org.Hs.eg.db", "AnnotationDbi", "decoupleR", "OmnipathR",
  "progeny", "dorothea", "slingshot", "presto", "zellkonverter"))

# 3. CellChat (GitHub master ≥ 2.2)
remotes::install_github("jinworks/CellChat")

# 4. MisTy (Bioconductor or GitHub)
BiocManager::install("mistyR")

# 5. NMF (CRAN, for cNMF post-processing)
install.packages("NMF")

# 6. monocle3 (cole-trapnell-lab)
remotes::install_github("cole-trapnell-lab/monocle3")

# 7. GeneSwitches (SGDDNB)
remotes::install_github("SGDDNB/GeneSwitches")
```

## 3. Hardware

| Stage | GPU | Notes |
|---|---|---|
| scVI / scANVI integration (10b, step25c) | required (16 GB+) | rapids-singlecell + cuML supported |
| cNMF (05b/05c per-patient, 18 global) | optional | CPU acceptable; ~6 h on 16 cores |
| Geneformer in-silico delete (fig8/03_perturb) | required (24 GB+) | V2-104M; ~24–48 h for 3 transitions × 500 cells |
| cell2location deconv (fig7/step03) | required (24 GB+) | NB model fitting |
| COMMOT cell-cell signaling (fig7/step06) | recommended | OT solver |
| All R plotting | none | seconds to minutes per figure |

Total disk for intermediate h5ads: ~120 GB. Final figure outputs: ~200 MB.

## 4. Pip-freeze caveat

`scst-pip-freeze.txt` (~460 packages) captures the **exact** transitive closure
that produced the published results, including version-pinned scientific stack
(`scvi-tools==1.x`, `scanpy==1.10.x`, `rapids-singlecell==0.10.x`, etc.).
For a less brittle install, the `*.yml` files give the user-requested top-level
deps and let mamba's resolver pick compatible versions.

## 5. R version note

R 4.3.3 (2024-02-29) was used throughout. Newer R 4.4.x has not been validated
against CellChat 2.2.0.9001 + Seurat 5.5.0 for the cellchat_*.R scripts.


## Environments added during revision

Three additional environments were used for the analyses added in revision.

`revision-python-pip-freeze.txt` — the Python environment in which every
`revision/**/*.py` script was run (a superset of `scst`; adds `lifelines`,
`matplotlib-venn`, `statsmodels` and `decoupler`).

`revision-R-packages.tsv` — R 4.5.3 with the package versions used for the
manuscript panels regenerated during revision (`ggplot2` 4.0.3, `patchwork`
1.3.2, `showtext` 0.9.8, `survival` 3.8.6, `survminer` 0.5.2) and for decontX
(`celda` 1.18.2).

`infercnv.yml` — inferCNV 1.26.0 on R 4.5.3, installed from conda binaries. The
source build of inferCNV did not compile on the development machine; the conda
binary package including `rjags` is what was used. This analysis ran on a
rented compute node that no longer exists, so the file is a specification
rather than an export.
