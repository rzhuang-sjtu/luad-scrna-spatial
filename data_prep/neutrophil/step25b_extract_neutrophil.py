"""Extract Neutrophil subset from luad_myeloid.h5ad → small h5ad for fast iteration.

Output: data/processed/luad_neutrophil_own_raw.h5ad
Keeps raw counts in .X (from layers['counts']), drops other layers/obsm to save space.
"""
import anndata as ad
import numpy as np
import scanpy as sc
import os, sys, time

t0 = time.time()
SRC = "${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad"
DST = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad"

print(f"[load] {SRC}")
a = sc.read_h5ad(SRC)
print(f"  full myeloid shape: {a.shape}")
print(f"  myeloid_subtype counts:\n{a.obs['myeloid_subtype'].value_counts()}")

mask = a.obs["myeloid_subtype"].astype(str) == "Neutrophil"
print(f"\n[subset] Neutrophil: n_cells = {int(mask.sum())}")
n = a[mask].copy()

# move raw counts → .X, drop log layer
if "counts" in n.layers:
    n.X = n.layers["counts"].copy()
    print("  set .X = layers['counts'] (raw)")
else:
    print("  WARN: no 'counts' layer — assuming .X is already raw")

# strip layers / obsm / uns to make file small + portable
n.layers.clear()
for k in list(n.obsm.keys()):
    del n.obsm[k]
for k in list(n.uns.keys()):
    del n.uns[k]

# verify integer-valued
xs = n.X[:200, :200].toarray() if hasattr(n.X, "toarray") else n.X[:200, :200]
print(f"  X dtype={n.X.dtype}  sample range=[{xs.min():.2f},{xs.max():.2f}]  integer={np.allclose(xs, np.round(xs))}")

# dataset breakdown
print("\n[dataset breakdown]")
print(n.obs["dataset"].value_counts())

print(f"\n[write] {DST}")
n.write_h5ad(DST, compression="gzip")
print(f"  saved: {os.path.getsize(DST)/1e6:.1f} MB  shape={n.shape}")
print(f"\nelapsed: {time.time()-t0:.1f}s")
