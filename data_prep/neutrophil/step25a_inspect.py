"""Inspect Salcher neutrophil_final.h5ad and own luad_myeloid.h5ad headers.

Goal: decide label-transfer strategy. Print only — no writes.
"""
import os
import sys
import anndata as ad
import numpy as np

SALCHER = "${DATA_ROOT}/High-resolution/neutrophil_final.h5ad"
MYELOID = "${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad"


def banner(msg):
    print("\n" + "=" * 60)
    print(msg)
    print("=" * 60)


def inspect(path, full=False):
    banner(path)
    print(f"size: {os.path.getsize(path) / 1e9:.2f} GB")
    a = ad.read_h5ad(path, backed="r")
    print(f"shape: {a.shape}")
    print(f"obs cols ({len(a.obs.columns)}): {list(a.obs.columns)}")
    print(f"obsm keys: {list(a.obsm.keys())}")
    print(f"layers: {list(a.layers.keys())}")
    print(f"uns keys: {list(a.uns.keys())}")
    if "X_scANVI" in a.obsm:
        x = a.obsm["X_scANVI"]
        print(f"X_scANVI shape={x.shape} dtype={x.dtype} nan={np.isnan(np.asarray(x)).any()}")
    # head of obs
    print("\nobs.head():")
    print(a.obs.head(3).to_string())
    # candidate label columns
    for col in a.obs.columns:
        s = a.obs[col]
        if s.dtype.name in ("category", "object"):
            uniq = s.astype(str).unique()
            if 1 < len(uniq) <= 30:
                print(f"  {col}: {sorted(uniq.tolist())[:25]}")
    return a


s = inspect(SALCHER)
m = inspect(MYELOID)

banner("Salcher barcode samples")
print(s.obs.index[:5].tolist())
banner("Myeloid barcode samples")
print(m.obs.index[:5].tolist())

banner("var_names overlap")
sv = set(s.var_names.astype(str))
mv = set(m.var_names.astype(str))
print(f"Salcher genes: {len(sv)}")
print(f"Myeloid genes: {len(mv)}")
print(f"intersection: {len(sv & mv)}")

# look for scANVI / scVI model files alongside Salcher h5ad
banner("model files near Salcher")
import glob
for p in glob.glob("${DATA_ROOT}/High-resolution/**", recursive=True):
    bn = os.path.basename(p)
    if any(k in bn.lower() for k in ("model", ".pt", "scvi", "scanvi", "_history")):
        print(p)
