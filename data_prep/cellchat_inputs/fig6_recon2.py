"""Pre-CellChat recon: data normalization + label completeness."""
import anndata as ad, pandas as pd, numpy as np, os

# 1) X normalization
print("=" * 60)
print("1) merged_annotated .X check")
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad", backed="r")
print(f"  shape: {a.shape}; layers: {list(a.layers.keys())}")
xs = a.X[:200, :200]
xs_a = xs.toarray() if hasattr(xs, "toarray") else xs
print(f"  X dtype={a.X.dtype}, sample range=[{xs_a.min():.3f}, {xs_a.max():.3f}]")
print(f"  integer? {np.allclose(xs_a, np.round(xs_a))}")
print(f"  median nonzero: {np.median(xs_a[xs_a>0]) if (xs_a>0).any() else 'na'}")

print("\n2) malignant h5ad: dominant_MP availability")
m = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad", backed="r")
print(f"  {m.shape}; dominant_MP counts:")
print(m.obs['dominant_MP'].value_counts().to_string())
print(f"  X range sample: ", end="")
xs2 = m.X[:200, :200]
xs2_a = xs2.toarray() if hasattr(xs2, "toarray") else xs2
print(f"[{xs2_a.min():.3f}, {xs2_a.max():.3f}]  integer? {np.allclose(xs2_a, np.round(xs2_a))}")

print("\n3) Macro_SPP1 label — fig4 refined metadata")
fig4_meta = "${WORK_ROOT}/luad_figures/fig4/panel_major_type_metadata.csv.gz"
if os.path.exists(fig4_meta):
    df = pd.read_csv(fig4_meta, nrows=5)
    print(f"  cols: {list(df.columns)[:25]}")
    full = pd.read_csv(fig4_meta, usecols=[c for c in df.columns if 'subtype' in c or 'major' in c or 'refined' in c.lower() or 'cell_id' in c.lower()])
    print(f"  full shape: {full.shape}")
    for col in full.columns:
        if col == 'cell_id' or col == 'barcode': continue
        vc = full[col].value_counts().head(20)
        print(f"  {col}:")
        print(vc.to_string())

print("\n4) neutrophil annotated obs barcodes sample")
neu = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad", backed="r")
print(f"  {neu.shape}, barcode head: {neu.obs.index[:3].tolist()}")
print(f"  neu_subtype dist: {neu.obs['neu_subtype'].value_counts().head().to_string()}")

print("\n5) merged barcode index format")
print(f"  merged barcode head: {a.obs.index[:3].tolist()}")
print(f"  malignant barcode head: {m.obs.index[:3].tolist()}")

print("\n6) tissue_type buckets in merged")
print(a.obs['tissue_type'].value_counts().to_string())
