"""Recon: assess what's needed for Figure 6 communication panels.

Check:
  1. step12 LIANA outputs — what's in step12_liana_all.csv?
  2. luad_tnk.h5ad / luad_myeloid.h5ad / luad_malignant_scored.h5ad — subtype labels
  3. luad_merged_annotated.h5ad — celltype_coarse for global types
  4. Are endothelial/fibroblast already labeled separately?
"""
import pandas as pd
import anndata as ad
import os

print("=" * 60)
print("1) step12 LIANA — what's there?")
print("=" * 60)
for f in ["${PROJECT_ROOT}/results/step12_liana_all.csv",
          "${PROJECT_ROOT}/results/step12_liana_mp_comparison.csv",
          "${PROJECT_ROOT}/results/step12_liana_mp3_neutrophil.csv"]:
    if os.path.exists(f):
        df = pd.read_csv(f, nrows=3)
        n = sum(1 for _ in open(f)) - 1
        print(f"\n{os.path.basename(f)}  ({os.path.getsize(f)/1e6:.1f} MB, ~{n} rows)")
        print(f"  cols: {list(df.columns)}")
        print(df.head(2).to_string())
        if 'source' in df.columns:
            full = pd.read_csv(f, usecols=['source','target'])
            print(f"  unique sources: {full['source'].unique()[:30]}")
            print(f"  unique targets: {full['target'].unique()[:30]}")

print("\n" + "=" * 60)
print("2) Per-compartment h5ad obs labels")
print("=" * 60)
for path in ["${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad",
             "${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad",
             "${PROJECT_ROOT}/data/processed/luad_tnk.h5ad",
             "${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad",
             "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad"]:
    if not os.path.exists(path): continue
    a = ad.read_h5ad(path, backed="r")
    print(f"\n{os.path.basename(path)}: {a.shape}")
    candidates = [c for c in a.obs.columns if any(k in c.lower() for k in
                  ['celltype','subtype','cluster','major','dominant','annot'])]
    for c in candidates:
        try:
            uniq = a.obs[c].astype(str).unique()
            if 1 < len(uniq) < 50:
                print(f"  {c}: {sorted(uniq.tolist())}")
        except Exception:
            pass
    # also check tissue_type
    if "tissue_type" in a.obs.columns:
        print(f"  tissue_type: {sorted(a.obs['tissue_type'].astype(str).unique().tolist())}")

print("\n" + "=" * 60)
print("3) endothelial / fibroblast — separate h5ads?")
print("=" * 60)
import glob
for d in ["${PROJECT_ROOT}/data/processed/", "${WORK_ROOT}/"]:
    for f in glob.glob(d + "*.h5ad"):
        bn = os.path.basename(f).lower()
        if any(k in bn for k in ['endo','fibro','stroma','epith']):
            print(f"  {f}  ({os.path.getsize(f)/1e9:.2f} GB)")

# Also check celltype_coarse from luad_merged_annotated for global-level labels
print("\n" + "=" * 60)
print("4) cell type coarse breakdown in luad_merged_annotated")
print("=" * 60)
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad", backed="r")
for col in ["celltype_coarse", "celltype_marker", "celltype_celltypist", "celltype_ct_coarse",
            "celltype_marker", "celltype_original_mapped"]:
    if col in a.obs.columns:
        print(f"\n  {col}:")
        print(a.obs[col].value_counts().head(20).to_string())
