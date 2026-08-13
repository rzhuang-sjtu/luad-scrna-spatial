#!/usr/bin/env python3
"""
Step 3b: Preprocess + Harmony integration + Leiden clustering

Input:  ${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad  (853,469 × 9,881, X=raw counts)
Output:  ${PROJECT_ROOT}/data/processed/luad_integrated.h5ad

Constraints:
  - X stays raw counts (float32 CSR) for downstream CopyKAT / cNMF
  - layers['counts']  = raw counts (redundant copy for downstream)
  - layers['lognorm'] = normalize_total(1e4)+log1p
  - obsm['X_pca'], obsm['X_pca_harmony'], obsm['X_umap']
  - obs['leiden_1.0'], obs['leiden_2.0']

Expected runtime:
  - GPU (rapids-singlecell):  3-5 min
  - CPU (fallback):           15-25 min
Peak memory ~30GB (CPU) / ~15GB RAM + 15GB VRAM (GPU).
"""
import gc
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad

# --- Try GPU (rapids-singlecell) ---
try:
    import rapids_singlecell as rsc
    import cupy as cp
    GPU = True
    _gpu_msg = f"rapids-singlecell {rsc.__version__} + cupy {cp.__version__}"
except ImportError as _e:
    GPU = False
    _gpu_msg = f"rapids-singlecell not installed ({_e.name}); fallback to CPU"

IN_PATH = Path("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad")
OUT_PATH = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")

N_TOP_GENES = 3000
N_PCS = 50
N_NEIGHBORS = 30
LEIDEN_RES = [1.0, 2.0]

sc.settings.verbosity = 3
sc.settings.n_jobs = 16


def _to_cpu_array(x):
    """cupy ndarray / cupy sparse → numpy / scipy sparse; leave others unchanged"""
    if not GPU:
        return x
    if isinstance(x, cp.ndarray):
        return x.get()
    try:
        import cupyx.scipy.sparse as cps
        if isinstance(x, cps.spmatrix):
            return x.get()
    except ImportError:
        pass
    return x


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = datetime.now()
    log("=" * 60)
    log(f"Step 3b Integrate started: {t0}")
    log(f"GPU: {GPU}  ({_gpu_msg})")
    log("=" * 60)

    log(f"Reading {IN_PATH}")
    adata = sc.read_h5ad(IN_PATH)
    log(f"  shape: {adata.shape}")
    log(f"  X dtype: {adata.X.dtype}, sparse: {hasattr(adata.X, 'toarray')}")
    log(f"  datasets: {dict(adata.obs['dataset'].value_counts())}")

    log("Save raw counts to layers['counts']")
    adata.layers["counts"] = adata.X.copy()

    log("normalize_total(target_sum=1e4) + log1p")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    log("Save log-normalized matrix to layers['lognorm']")
    adata.layers["lognorm"] = adata.X.copy()

    log(f"HVG n_top_genes={N_TOP_GENES}, batch_key='dataset', subset=False")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=N_TOP_GENES, batch_key="dataset", subset=False
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    log(f"  HVG count: {n_hvg}")

    log("HVG subset → scale(max_value=10) → PCA")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log(f"HVG subset shape: {adata_hvg.shape}")

    if GPU:
        log("  [GPU] anndata_to_GPU → rsc.pp.scale → rsc.pp.pca")
        rsc.get.anndata_to_GPU(adata_hvg)
        rsc.pp.scale(adata_hvg, max_value=10)
        rsc.pp.pca(adata_hvg, n_comps=N_PCS)
        adata.obsm["X_pca"] = _to_cpu_array(adata_hvg.obsm["X_pca"]).astype(np.float32)
        if "pca" in adata_hvg.uns:
            pca_uns = adata_hvg.uns["pca"]
            adata.uns["pca"] = {
                "variance": _to_cpu_array(pca_uns.get("variance", np.array([]))),
                "variance_ratio": _to_cpu_array(pca_uns.get("variance_ratio", np.array([]))),
            }
    else:
        log("  [CPU] sc.pp.scale → sc.tl.pca(svd_solver=randomized)")
        sc.pp.scale(adata_hvg, max_value=10)
        log(f"scale done, dtype={adata_hvg.X.dtype}")
        sc.tl.pca(adata_hvg, n_comps=N_PCS, zero_center=True, svd_solver="randomized")
        adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
        adata.uns["pca"] = {
            "variance": adata_hvg.uns["pca"]["variance"].copy(),
            "variance_ratio": adata_hvg.uns["pca"]["variance_ratio"].copy(),
        }
    log(f"PCA done, X_pca shape: {adata.obsm['X_pca'].shape}")

    # Free scaled matrix immediately (including GPU memory)
    del adata_hvg
    gc.collect()
    if GPU:
        cp.get_default_memory_pool().free_all_blocks()
    log("Freed scaled matrix")

    if GPU:
        log("[GPU] rsc.pp.harmony_integrate (key='dataset')")
        rsc.pp.harmony_integrate(adata, key="dataset", basis="X_pca",
                                 adjusted_basis="X_pca_harmony")
    else:
        log("[CPU] sce.pp.harmony_integrate (key='dataset')")
        import scanpy.external as sce
        sce.pp.harmony_integrate(adata, key="dataset", basis="X_pca",
                                 adjusted_basis="X_pca_harmony")
    # If cupy array, convert to numpy
    adata.obsm["X_pca_harmony"] = _to_cpu_array(adata.obsm["X_pca_harmony"])
    log(f"  X_pca_harmony shape: {adata.obsm['X_pca_harmony'].shape}")

    if GPU:
        log(f"[GPU] rsc.pp.neighbors (use_rep='X_pca_harmony', n_neighbors={N_NEIGHBORS})")
        rsc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=N_NEIGHBORS)
        log("[GPU] rsc.tl.umap")
        rsc.tl.umap(adata)
        log(f"  X_umap shape: {adata.obsm['X_umap'].shape}")
        for res in LEIDEN_RES:
            key = f"leiden_{res}"
            log(f"[GPU] rsc.tl.leiden resolution={res} → obs['{key}']")
            rsc.tl.leiden(adata, resolution=res, key_added=key)
            n_clusters = adata.obs[key].nunique()
            log(f"  {key}: {n_clusters} clusters")
        # rsc may store obsm / obsp as cupy; convert back to CPU for h5ad write
        for k in list(adata.obsm.keys()):
            adata.obsm[k] = _to_cpu_array(adata.obsm[k])
        for k in list(adata.obsp.keys()):
            adata.obsp[k] = _to_cpu_array(adata.obsp[k])
        cp.get_default_memory_pool().free_all_blocks()
    else:
        log(f"[CPU] sc.pp.neighbors (use_rep='X_pca_harmony', n_neighbors={N_NEIGHBORS})")
        sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=N_NEIGHBORS)
        log("[CPU] sc.tl.umap")
        sc.tl.umap(adata)
        log(f"  X_umap shape: {adata.obsm['X_umap'].shape}")
        for res in LEIDEN_RES:
            key = f"leiden_{res}"
            log(f"[CPU] sc.tl.leiden resolution={res} → obs['{key}']")
            sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph",
                         n_iterations=2, directed=False)
            n_clusters = adata.obs[key].nunique()
            log(f"  {key}: {n_clusters} clusters")

    log("Restore X ← layers['counts'] (raw counts)")
    adata.X = adata.layers["counts"].copy()
    log(f"  X dtype: {adata.X.dtype}, sparse: {hasattr(adata.X, 'toarray')}")

    log(f"Writing {OUT_PATH} (gzip)")
    adata.write_h5ad(OUT_PATH, compression="gzip")
    size_gb = OUT_PATH.stat().st_size / 1e9
    log(f"Done: {OUT_PATH} ({size_gb:.2f} GB)")

    log("=" * 60)
    log("Final checks")
    log("=" * 60)
    log(f"shape: {adata.shape}")
    log(f"X = raw counts, dtype={adata.X.dtype}, sparse={hasattr(adata.X, 'toarray')}")
    log(f"layers: {list(adata.layers.keys())}")
    log(f"obsm: {list(adata.obsm.keys())}")
    for res in LEIDEN_RES:
        key = f"leiden_{res}"
        nc = adata.obs[key].nunique()
        log(f"{key}: {nc} clusters")
    log("Cells per dataset:")
    for ds, n in adata.obs["dataset"].value_counts().items():
        log(f"  {ds}: {n:,}")

    elapsed = (datetime.now() - t0).total_seconds() / 60
    log(f"\nElapsed: {elapsed:.1f} min")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
