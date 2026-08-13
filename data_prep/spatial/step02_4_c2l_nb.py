"""
Step 2.4: Train cell2location NB regression on the unified reference to obtain
per-gene per-cell-type expression signatures (inf_aver) for spot deconvolution.

Inputs:
  ${DATA_ROOT}/ST/results/step02_reference/unified_reference.h5ad  (raw counts)
Outputs:
  ${DATA_ROOT}/ST/results/step02_reference/c2l_reg_model/          trained model dir
  ${DATA_ROOT}/ST/results/step02_reference/sc_with_signatures.h5ad reference + posterior
  ${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv            cell-type × gene signature matrix
  ${DATA_ROOT}/ST/results/step02_reference/elbo.png                training curve

Following cell2location reference tutorial:
  - filter_genes via cell2location.utils.filtering.filter_genes (count>=5, pct>=0.03, nonz_mean>=1.12)
  - batch_key = 'dataset', labels_key = 'cell_type_fine'
  - categorical_covariate_keys = ['patient_id'] for patient-level batch effect
  - max_epochs=250, batch_size=2500, lr=0.002, GPU
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF_IN  = "${DATA_ROOT}/ST/results/step02_reference/unified_reference.h5ad"
OUT     = Path("${DATA_ROOT}/ST/results/step02_reference")
MODEL_D = OUT / "c2l_reg_model"
MODEL_D.mkdir(parents=True, exist_ok=True)

t0 = time.time()
print(f"[load] {REF_IN}")
adata = sc.read_h5ad(REF_IN)
print(f"   loaded: {adata.shape}")
print(f"   X type: {type(adata.X).__name__}, dtype={adata.X.dtype}")

# Coerce X to sparse int (c2l expects integer counts)
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)
adata.X = adata.X.astype("float32")  # c2l accepts float32 ints

# c2l-style gene filter
import cell2location
from cell2location.utils.filtering import filter_genes
print("[filter_genes] cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12")
selected = filter_genes(adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata = adata[:, selected].copy()
print(f"   after filter: {adata.shape}")

# Setup AnnData for RegressionModel
from cell2location.models import RegressionModel

# fix dtype of categorical covariates
adata.obs["dataset"] = adata.obs["dataset"].astype(str)
adata.obs["patient_id"] = adata.obs["patient_id"].astype(str)
adata.obs["cell_type_fine"] = adata.obs["cell_type_fine"].astype(str)

print("[setup_anndata]")
RegressionModel.setup_anndata(
    adata=adata,
    batch_key="dataset",
    labels_key="cell_type_fine",
    categorical_covariate_keys=["patient_id"],
)

print("[train] RegressionModel max_epochs=250, batch_size=2500, lr=0.002, GPU")
mod = RegressionModel(adata)
import torch
use_gpu = torch.cuda.is_available()
print(f"   CUDA available: {use_gpu}")

train_kwargs = dict(max_epochs=250, batch_size=2500, train_size=1.0, lr=0.002)
# scvi-tools 1.x: accelerator only; do NOT pass devices (c2l passes its own); do NOT pass use_gpu
try:
    mod.train(accelerator="gpu" if use_gpu else "cpu", **train_kwargs)
except TypeError:
    # very old c2l fallback
    mod.train(**train_kwargs)

# Save ELBO curve
fig, ax = plt.subplots(figsize=(5,3))
mod.plot_history(20)
fig.tight_layout()
fig.savefig(OUT / "elbo.png", dpi=150)
plt.close(fig)
print(f"[plot] elbo saved to {OUT/'elbo.png'}")

# Posterior
print("[export_posterior]")
adata = mod.export_posterior(
    adata,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500},
)

# Save model + reference + signature
print("[save model]")
mod.save(str(MODEL_D), overwrite=True)

# Extract inf_aver: per-gene per-cell-type expression
key = "means_per_cluster_mu_fg"
if key in adata.varm:
    inf_aver = adata.varm[key].copy()
    inf_aver.columns = [c.replace("means_per_cluster_mu_fg_","")
                        for c in inf_aver.columns]
elif f"q05_per_cluster_mu_fg" in adata.varm:
    inf_aver = adata.varm["q05_per_cluster_mu_fg"].copy()
else:
    # fallback path used by some c2l versions
    inf_aver = pd.DataFrame(
        adata.varm[list(adata.varm.keys())[0]],
        index=adata.var_names,
    )

inf_aver.to_csv(OUT / "inf_aver.csv")
print(f"[save] inf_aver.csv ({inf_aver.shape[0]} genes × {inf_aver.shape[1]} cell types)")
print(f"   first cell types: {list(inf_aver.columns)[:8]}")

# Save reference with posterior
ref_out = OUT / "sc_with_signatures.h5ad"
adata.write_h5ad(str(ref_out), compression="gzip")
print(f"[save] {ref_out}  ({ref_out.stat().st_size/1e9:.2f} GB)")

print(f"\n[done] elapsed: {(time.time()-t0)/60:.1f} min")
