"""Per-cell iLISI and cLISI before and after Harmony integration (Fig. S1F).

Recomputed during revision: LISI is reported split into iLISI (batch mixing)
and cLISI (cell-type purity) rather than as one aggregate score.
Writes s1_lisi_scores.csv and s1_lisi_per_cell.csv.gz, both read by
plotting/figS1/figure_s1_main.R.
"""
import numpy as np
import pandas as pd
import scanpy as sc
from harmonypy.lisi import compute_lisi

adata = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
print(f"Full adata: {adata.n_obs} cells")

N_SUBSAMPLE = 100_000
rng = np.random.default_rng(0)
if adata.n_obs > N_SUBSAMPLE:
    idx = rng.choice(adata.n_obs, size=N_SUBSAMPLE, replace=False)
    idx.sort()
else:
    idx = np.arange(adata.n_obs)
print(f"Subsampled to {len(idx)} cells")

pca_pre = adata.obsm["X_pca"][idx, :30]
pca_post = adata.obsm["X_pca_harmony"][idx, :30]
meta = pd.DataFrame({"dataset": adata.obs["dataset"].values[idx]},
                    index=adata.obs_names[idx])

print("Computing pre-Harmony LISI...")
lisi_pre = compute_lisi(pca_pre, meta, ["dataset"])
print("Computing post-Harmony LISI...")
lisi_post = compute_lisi(pca_post, meta, ["dataset"])

out = pd.DataFrame({
    "stage": ["Before Harmony", "After Harmony"],
    "median_LISI": [np.median(lisi_pre), np.median(lisi_post)],
    "mean_LISI": [np.mean(lisi_pre), np.mean(lisi_post)],
    "n_cells": [len(idx), len(idx)],
})
print(out)

out.to_csv("${WORK_ROOT}/luad_figures/fig_s1/s1_lisi_scores.csv", index=False)

lisi_df = pd.DataFrame({
    "LISI_pre": lisi_pre.flatten(),
    "LISI_post": lisi_post.flatten(),
})
lisi_df.to_csv("${WORK_ROOT}/luad_figures/fig_s1/s1_lisi_per_cell.csv.gz",
               compression="gzip", index=False)
print("Saved outputs.")
