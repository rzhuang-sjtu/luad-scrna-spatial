"""Strip problematic uns entries (log1p with None) from cellchat h5ads.

zellkonverter via basilisk anndata 0.10.2 cannot decode encoding_type='null'
which scanpy writes for adata.uns['log1p']={'base': None}.
"""
import anndata as ad
import os, time
from pathlib import Path

t0 = time.time()
DATA = Path("${PROJECT_ROOT}/data/processed")
files = ["cellchat_input_all.h5ad",
         "cellchat_input_Normal.h5ad",
         "cellchat_input_Tumor.h5ad",
         "cellchat_input_Metastasis.h5ad"]

for f in files:
    p = DATA / f
    print(f"[fix] {p}")
    a = ad.read_h5ad(p)
    # nuke uns entirely (we don't need any uns for CellChat)
    keys = list(a.uns.keys())
    for k in keys:
        del a.uns[k]
    # also drop layers (raw counts already moved earlier; keep .X = lognorm)
    a.layers.clear()
    # drop .raw
    if a.raw is not None:
        a.raw = None
    print(f"  shape={a.shape}, uns={list(a.uns.keys())}, layers={list(a.layers.keys())}")
    a.write_h5ad(p, compression="gzip")
    print(f"  rewritten {os.path.getsize(p)/1e6:.0f} MB")

print(f"\nelapsed: {time.time()-t0:.1f}s")
