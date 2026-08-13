"""Step 2.0 + 2.1: verify raw counts in luad_merged_raw.h5ad and presence of key Fig7 genes."""
import os, numpy as np, scipy.sparse as sp, scanpy as sc, pandas as pd

PROC = os.path.expanduser("~/luad/data/processed")
RAW   = os.path.join(PROC, "luad_merged_raw.h5ad")
ANNO  = os.path.join(PROC, "luad_merged_annotated.h5ad")

print(f"== {RAW} ==")
a = sc.read_h5ad(RAW, backed="r")
print(f"shape: {a.shape}")
print(f"obs cols: {a.obs.columns.tolist()}")
print(f"var cols: {a.var.columns.tolist()}")
print(f"first 5 var_names: {a.var_names[:5].tolist()}")

# Sample 3000 random cells, fetch X, check if integer
rng = np.random.default_rng(0)
idx = sorted(rng.choice(a.n_obs, size=min(3000, a.n_obs), replace=False))
X_sample = a.X[idx]
if sp.issparse(X_sample):
    data = X_sample.data
else:
    data = np.asarray(X_sample).ravel()
nz = data[data != 0]
print(f"X type: {type(X_sample).__name__}, dtype: {data.dtype}")
print(f"nonzero values: n={nz.size}, min={nz.min()}, max={nz.max()}, mean={nz.mean():.3f}")
is_int = np.all(np.equal(np.mod(nz, 1), 0))
print(f"all values integer? {is_int}")
print(f"first 20 nonzero values: {nz[:20]}")
a.file.close()

# Compare obs barcodes vs annotated
print(f"\n== {ANNO} (header check) ==")
b = sc.read_h5ad(ANNO, backed="r")
print(f"annotated shape: {b.shape}")
same_n   = (a.shape[0] == b.shape[0])
same_var = (a.shape[1] == b.shape[1])
# Compare first 10 obs and var names without copying full arrays
print(f"obs/var dims match annotated: obs={same_n}, var={same_var}")
b.file.close()

# Re-open RAW to also probe key genes
a = sc.read_h5ad(RAW, backed="r")
KEY_GENES = [
    # Fig 7 ligands / receptors / pathways
    "OSM","OSMR","LIFR","IL6ST",                       # OSM pathway
    "IL1A","IL1B","IL1R1","IL1R2","IL1RAP","IL1RN",   # IL1 pathway
    "SPP1","CD44",                                     # SPP1
    "MMP9","MMP2","MMP1","MMP7","MMP14",               # MMPs
    "CXCL8","CXCL1","CXCL2","CXCL5","CXCR1","CXCR2",   # neutro chemokines
    "PLAU","PLAUR","PLAT",                             # PLAU
    "TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2",         # TGFb (sanity)
    "VEGFA","VEGFB",                                    # angio
    "IFNG","STAT1","JAK1","JAK2","STAT3",              # IFN/JAK-STAT
    "NFKB1","NFKBIA","RELA","RELB","REL",              # NF-kB
    # Tumor program TFs
    "ATF3","FOSB","JUN","JUNB","FOS","ELF3","HIF1A",
    # Macro markers
    "C1QA","C1QB","C1QC","FCN1","FOLR2","MARCO",
    # Neutro core
    "S100A8","S100A9","FCGR3B","CSF3R","ELANE","MPO","PADI4",
    # Malignant epithelial
    "EPCAM","KRT7","KRT8","KRT18","KRT19","SFTPC","SFTPB","NAPSA","KRT5","TP63",
    "LAMC2","VIM","CDH1","CDH2","SNAI1","SNAI2","ZEB1","ZEB2",
]
present = [g for g in KEY_GENES if g in a.var_names]
missing = [g for g in KEY_GENES if g not in a.var_names]
print(f"\n== key gene presence: {len(present)}/{len(KEY_GENES)} present ==")
print(f"PRESENT ({len(present)}): {' '.join(present)}")
print(f"\nMISSING ({len(missing)}): {' '.join(missing)}")

a.file.close()
print("\nDONE")
