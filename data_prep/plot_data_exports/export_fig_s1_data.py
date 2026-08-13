#!/usr/bin/env python
"""
Supplementary Figure 1 data export
Extra QC metrics + pre-Harmony UMAP + kBET + cluster×sample proportions

Input: ~/luad/data/processed/luad_integrated.h5ad (read-only)
Output: ${WORK_ROOT}/luad_figures/fig_s1/*.csv
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path
from datetime import datetime


IN_PATH = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
OUT_DIR = Path("${WORK_ROOT}/luad_figures/fig_s1")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log("Reading h5ad")
    adata = sc.read_h5ad(IN_PATH)
    log(f"  shape: {adata.shape}")

    # 1. QC metrics (recompute from raw counts)
    log("Compute QC metrics")

    X = adata.layers["counts"]  # raw counts

    # total counts per cell
    if sp.issparse(X):
        n_counts = np.array(X.sum(axis=1)).flatten()
        n_genes = np.array((X > 0).sum(axis=1)).flatten()
    else:
        n_counts = X.sum(axis=1)
        n_genes = (X > 0).sum(axis=1)

    # mito genes
    mito_genes = adata.var_names.str.startswith("MT-")
    n_mito = mito_genes.sum()
    log(f"MT- genes: {n_mito}")

    if n_mito > 0:
        if sp.issparse(X):
            mito_counts = np.array(X[:, mito_genes].sum(axis=1)).flatten()
        else:
            mito_counts = X[:, mito_genes].sum(axis=1)
        pct_mito = (mito_counts / (n_counts + 1e-10)) * 100
    else:
        pct_mito = np.zeros(adata.n_obs)
        log("No MT- genes found; pct_mito all 0")

    # ribo genes
    ribo_genes = adata.var_names.str.match("^RP[SL]\\d")
    n_ribo = ribo_genes.sum()
    if n_ribo > 0 and sp.issparse(X):
        ribo_counts = np.array(X[:, ribo_genes].sum(axis=1)).flatten()
        pct_ribo = (ribo_counts / (n_counts + 1e-10)) * 100
    else:
        pct_ribo = np.zeros(adata.n_obs)

    log(f"  n_counts range: {n_counts.min():.0f} - {n_counts.max():.0f}")
    log(f"  n_genes range: {n_genes.min()} - {n_genes.max()}")
    log(f"  pct_mito range: {pct_mito.min():.2f} - {pct_mito.max():.2f}")

    # 2. Pre-Harmony UMAP (Panel F, left)
    log("Compute pre-Harmony UMAP")

    # Neighbors → UMAP from uncorrected PCA (X_pca)
    import anndata as ad
    adata_tmp = ad.AnnData(
        obs=pd.DataFrame(index=adata.obs_names),
        obsm={"X_pca": adata.obsm["X_pca"]}
    )
    sc.pp.neighbors(adata_tmp, use_rep="X_pca", n_neighbors=30)
    sc.tl.umap(adata_tmp)
    umap_pre = adata_tmp.obsm["X_umap"]
    del adata_tmp; gc.collect()
    log(f"  pre-Harmony UMAP done")

    # 3. Assemble full metadata
    log("Assemble S1 metadata")

    meta = adata.obs[["dataset", "sample_id", "patient_id", "tissue_type",
                       "celltype_coarse", "leiden_1.0", "leiden_2.0"]].copy()

    meta["n_counts"] = n_counts
    meta["n_genes"] = n_genes
    meta["pct_mito"] = np.round(pct_mito, 2)
    meta["pct_ribo"] = np.round(pct_ribo, 2)

    # doublet_score (if present)
    if "doublet_score" in adata.obs.columns:
        meta["doublet_score"] = adata.obs["doublet_score"].values

    # post-Harmony UMAP
    meta["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    meta["UMAP_2"] = adata.obsm["X_umap"][:, 1]

    # pre-Harmony UMAP
    meta["UMAP_pre_1"] = umap_pre[:, 0]
    meta["UMAP_pre_2"] = umap_pre[:, 1]

    # 4. Export

    # 4a. Full metadata
    log("Export s1_cell_metadata.csv.gz")
    meta.to_csv(OUT_DIR / "s1_cell_metadata.csv.gz", compression="gzip")
    log(f"  {len(meta):,} cells × {len(meta.columns)} cols")

    # 4b. QC summary per sample (Panels A–C)
    log("Export s1_qc_per_sample.csv")
    qc_sample = meta.groupby("sample_id").agg(
        dataset=("dataset", "first"),
        tissue_type=("tissue_type", "first"),
        n_cells=("n_counts", "size"),
        median_counts=("n_counts", "median"),
        median_genes=("n_genes", "median"),
        median_pct_mito=("pct_mito", "median"),
    ).reset_index()
    qc_sample.to_csv(OUT_DIR / "s1_qc_per_sample.csv", index=False)
    log(f"  {len(qc_sample)} samples")

    # 4c. Cluster × sample proportions (Panel H)
    log("Export s1_cluster_by_sample.csv")
    ct = pd.crosstab(meta["sample_id"], meta["leiden_1.0"])
    ct_pct = ct.div(ct.sum(axis=1), axis=0)
    ct_pct.to_csv(OUT_DIR / "s1_cluster_by_sample.csv")
    log(f"  {ct_pct.shape}")

    # 4d. Cluster × dataset proportions (Panel H backup)
    log("Export s1_cluster_by_dataset.csv")
    ct_ds = pd.crosstab(meta["dataset"], meta["leiden_1.0"])
    ct_ds_pct = ct_ds.div(ct_ds.sum(axis=1), axis=0)
    ct_ds_pct.to_csv(OUT_DIR / "s1_cluster_by_dataset.csv")

    # 4e. Simplified kBET: per-cluster batch entropy
    # (true kBET needs an R package; Shannon entropy used as approximation)
    log("Compute batch mixing entropy (kBET approximation)")
    from scipy.stats import entropy

    def batch_entropy(group):
        counts = group["dataset"].value_counts()
        p = counts / counts.sum()
        return entropy(p, base=2)

    # per-cluster entropy
    cluster_entropy = meta.groupby("leiden_1.0").apply(batch_entropy)
    max_entropy = np.log2(meta["dataset"].nunique())  # max entropy under perfect mixing

    ent_df = pd.DataFrame({
        "cluster": cluster_entropy.index,
        "batch_entropy": cluster_entropy.values,
        "max_entropy": max_entropy,
        "normalized_entropy": cluster_entropy.values / max_entropy,
    })
    ent_df.to_csv(OUT_DIR / "s1_batch_entropy.csv", index=False)
    log(f"  mean normalized entropy: {ent_df['normalized_entropy'].mean():.3f}")

    # 4f. Global kBET approximation (pre vs post Harmony)
    log("Compute global batch mixing (pre vs post Harmony)")

    from sklearn.neighbors import NearestNeighbors
    from scipy.stats import chi2 as chi2_dist

    def compute_kbet_like(coords, batch_labels, n_neighbors=50, n_sample=5000):
        """Simplified kBET: sample → k nearest neighbours → chi-square on neighbour batch fractions vs global"""
        np.random.seed(42)
        idx = np.random.choice(len(coords), min(n_sample, len(coords)), replace=False)
        coords_sub = coords[idx]
        batch_sub = batch_labels[idx]

        nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=-1)
        nn.fit(coords)
        neighbors = nn.kneighbors(coords_sub, return_distance=False)

        global_freq = pd.Series(batch_labels).value_counts(normalize=True)
        df = len(global_freq) - 1
        rejections = 0

        for i in range(len(idx)):
            neighbor_batches = batch_labels[neighbors[i]]
            local_freq = pd.Series(neighbor_batches).value_counts(normalize=True)
            chi2_stat = 0
            for b in global_freq.index:
                obs = local_freq.get(b, 0)
                exp = global_freq[b]
                chi2_stat += ((obs - exp) ** 2) / exp
            p_value = 1 - chi2_dist.cdf(chi2_stat * n_neighbors, df)
            if p_value < 0.05:
                rejections += 1

        return rejections / len(idx)

    batch_arr = meta["dataset"].values
    pre_coords = meta[["UMAP_pre_1", "UMAP_pre_2"]].values
    post_coords = meta[["UMAP_1", "UMAP_2"]].values

    # PCA space is more accurate
    pre_pca = adata.obsm["X_pca"]
    post_pca = adata.obsm["X_pca_harmony"]

    log("  computing pre-Harmony kBET...")
    kbet_pre = compute_kbet_like(pre_pca, batch_arr)
    log(f"  pre-Harmony rejection rate: {kbet_pre:.3f}")

    log("  computing post-Harmony kBET...")
    kbet_post = compute_kbet_like(post_pca, batch_arr)
    log(f"  post-Harmony rejection rate: {kbet_post:.3f}")

    kbet_df = pd.DataFrame({
        "stage": ["Before Harmony", "After Harmony"],
        "rejection_rate": [kbet_pre, kbet_post],
    })
    kbet_df.to_csv(OUT_DIR / "s1_kbet_rejection.csv", index=False)

    log("=" * 50)
    log("S1 export done:")
    for f in sorted(OUT_DIR.iterdir()):
        sz = f.stat().st_size / 1e6
        log(f"  {f.name}: {sz:.1f} MB")
    log("=" * 50)


if __name__ == "__main__":
    main()
