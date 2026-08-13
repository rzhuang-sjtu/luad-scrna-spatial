#!/usr/bin/env python
"""
Fig S1 Panel I extra data: UMAP after regressing out n_counts
Compare before vs after regress to show technical variation did not drive celltype clustering

Output: ${WORK_ROOT}/luad_figures/fig_s1/s1_regress_umap.csv.gz
"""
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
from datetime import datetime

IN_PATH = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
OUT_DIR = Path("${WORK_ROOT}/luad_figures/fig_s1")


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log("Reading h5ad")
    adata = sc.read_h5ad(IN_PATH)
    log(f"  shape: {adata.shape}")

    # -- regress out n_counts → re-run PCA → Harmony → UMAP --
    log("Prepare regress_out analysis")

    # Start from lognorm layer
    adata.X = adata.layers["lognorm"].copy()

    # Need total counts per cell
    import scipy.sparse as sp
    raw = adata.layers["counts"]
    if sp.issparse(raw):
        n_counts = np.array(raw.sum(axis=1)).flatten()
    else:
        n_counts = raw.sum(axis=1)
    adata.obs["n_counts"] = n_counts

    # HVG subset
    hvg_mask = adata.var["highly_variable"].values
    adata_hvg = adata[:, hvg_mask].copy()
    log(f"  HVG subset: {adata_hvg.shape}")

    # regress out
    log("regress_out(['n_counts'])...")
    sc.pp.regress_out(adata_hvg, ["n_counts"])
    log("regress_out done")

    # scale
    sc.pp.scale(adata_hvg, max_value=10)

    # PCA
    log("PCA (regressed)")
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="randomized")

    # Harmony — call harmonypy directly (scanpy.external's .T breaks with PyTorch backend)
    log("Harmony (regressed)")
    import harmonypy as hm
    ho = hm.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, "dataset")
    z = np.asarray(ho.Z_corr)
    if z.shape[0] != adata_hvg.n_obs:
        z = z.T
    adata_hvg.obsm["X_pca_harmony"] = z

    # UMAP
    log("UMAP (regressed)")
    sc.pp.neighbors(adata_hvg, use_rep="X_pca_harmony", n_neighbors=30)
    sc.tl.umap(adata_hvg)

    umap_regress = adata_hvg.obsm["X_umap"]
    del adata_hvg; gc.collect()

    log("Export s1_regress_umap.csv.gz")
    df = pd.DataFrame({
        "UMAP_regress_1": umap_regress[:, 0],
        "UMAP_regress_2": umap_regress[:, 1],
        "celltype_coarse": adata.obs["celltype_coarse"].values,
    }, index=adata.obs_names)
    df.to_csv(OUT_DIR / "s1_regress_umap.csv.gz", compression="gzip")
    log(f"  {len(df):,} cells saved")
    log("Done")


if __name__ == "__main__":
    main()
