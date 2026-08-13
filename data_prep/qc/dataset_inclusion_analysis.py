#!/usr/bin/env python3
# Flattened code cells from dataset_inclusion_analysis.ipynb

import scanpy as sc
import pandas as pd

paths = {
    'GSE164789': r"${DATA_ROOT}/cleaned/GSE164789_LUAD_clean.h5ad",
    'GSE189357': r"${DATA_ROOT}/cleaned/GSE189357_LUAD_clean.h5ad",
}

for name, path in paths.items():
    adata = sc.read_h5ad(path, backed='r')
    print(f"\n{'='*60}")
    print(f"{name}: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  layers: {list(adata.layers.keys())}")
    print(f"  obs columns: {adata.obs.columns.tolist()}")

    # Required field check
    for col in ['sample_id', 'patient_id', 'tissue_type', 'dataset']:
        if col in adata.obs.columns:
            print(f"  {col}: {adata.obs[col].unique().tolist()[:10]}")
        else:
            print(f"  {col}: MISSING!")

    if 'n_genes_by_counts' in adata.obs.columns:
        print(f"  median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    if 'total_counts' in adata.obs.columns:
        print(f"  median UMI/cell: {adata.obs['total_counts'].median():.0f}")

    adata.file.close()


adata = sc.read_h5ad(
    r"${DATA_ROOT}/cleaned/GSE189357_LUAD_clean.h5ad",
    backed='r'
)
print("stage:", adata.obs['stage'].value_counts().to_dict())
print("chemo:", adata.obs['chemotherapy'].value_counts().to_dict())
print("samples per patient:")
print(adata.obs.groupby('patient_id')['sample_id'].nunique())
adata.file.close()


adata = sc.read_h5ad(
    r"${DATA_ROOT}/cleaned/GSE253013_clean.h5ad",
    backed='r'
)
print(f"GSE253013: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"layers: {list(adata.layers.keys())}")
print(f"obs columns: {adata.obs.columns.tolist()}")
for col in ['sample_id', 'patient_id', 'tissue_type', 'dataset']:
    if col in adata.obs.columns:
        print(f"{col}: {adata.obs[col].value_counts().to_dict()}")
    else:
        print(f"{col}: MISSING!")
if 'n_genes_by_counts' in adata.obs.columns:
    print(f"median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"median UMI/cell: {adata.obs['total_counts'].median():.0f}")
adata.file.close()


adata = sc.read_h5ad(
    r"${DATA_ROOT}/cleaned/GSE253013_clean.h5ad",
    backed='r'
)
print("cd45_sorted distribution:")
print(adata.obs['cd45_sorted'].value_counts())
print("\nBy patient_id x cd45_sorted:")
print(adata.obs.groupby(['patient_id','cd45_sorted']).size().unstack(fill_value=0))
print("\nExisting cell_type annotations:")
print(adata.obs['cell_type'].value_counts().head(20))
adata.file.close()


import scanpy as sc
import pandas as pd

paths = {
    'GSE131907': r"${DATA_ROOT}/cleaned/GSE131907_clean.h5ad",
    'GSE123902': r"${DATA_ROOT}/cleaned/GSE123902_clean.h5ad",
    'GSE148071': r"${DATA_ROOT}/cleaned/GSE148071_LUAD_clean.h5ad",
    'GSE143423': r"${DATA_ROOT}/cleaned/GSE143423_LUAD_clean.h5ad",
    'GSE149655': r"${DATA_ROOT}/cleaned/GSE149655/GSE149655_LUAD_clean.h5ad",
}

for name, path in paths.items():
    adata = sc.read_h5ad(path, backed='r')
    print(f"\n{'='*60}")
    print(f"{name}: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  layers: {list(adata.layers.keys())}")

    for col in ['sample_id', 'patient_id', 'tissue_type', 'dataset']:
        if col in adata.obs.columns:
            vals = adata.obs[col].value_counts()
            print(f"  {col}: {len(vals)} unique | {vals.to_dict()}" if len(vals) <= 10
                  else f"  {col}: {len(vals)} unique | top5: {vals.head(5).to_dict()}")
        else:
            print(f"  {col}: MISSING!")

    if 'n_genes_by_counts' in adata.obs.columns:
        print(f"  median genes: {adata.obs['n_genes_by_counts'].median():.0f}")
        print(f"  median UMI: {adata.obs['total_counts'].median():.0f}")

    for col in ['cell_type', 'celltype', 'cell_type_major']:
        if col in adata.obs.columns:
            print(f"  {col}: {adata.obs[col].value_counts().head(10).to_dict()}")

    adata.file.close()


import scanpy as sc

# === Fix 1: GSE148071 add lognorm ===
path_148071 = r"${DATA_ROOT}/cleaned/GSE148071_LUAD_clean.h5ad"
adata = sc.read_h5ad(path_148071)
adata.layers['lognorm'] = adata.X.copy()  # provisional copy, may be replaced below
import numpy as np
max_val = adata.X[:100].toarray().max() if hasattr(adata.X, 'toarray') else adata.X[:100].max()
print(f"GSE148071 X max (first 100 cells): {max_val}")
# If X holds raw counts, regenerate lognorm
if max_val > 50:
    print("X is raw counts, regenerating lognorm...")
    adata.layers['lognorm'] = None
    temp = adata.copy()
    sc.pp.normalize_total(temp, target_sum=1e4)
    sc.pp.log1p(temp)
    adata.layers['lognorm'] = temp.X.copy()
    del temp
    print("lognorm layer added")
else:
    print("X already normalized, using as lognorm")
adata.write(path_148071)
print("GSE148071 saved")
del adata

# === Fix 2: GSE149655 obs_names dedup ===
path_149655 = r"${DATA_ROOT}/cleaned/GSE149655/GSE149655_LUAD_clean.h5ad"
adata = sc.read_h5ad(path_149655)
print(f"\nGSE149655 obs_names unique: {adata.obs_names.is_unique}")
adata.obs_names_make_unique()
print(f"After fix: {adata.obs_names.is_unique}")
adata.write(path_149655)
print("GSE149655 saved")
del adata


import scanpy as sc
import numpy as np

# === Fix 1: GSE148071 add lognorm ===
path_148071 = r"${DATA_ROOT}/cleaned/GSE148071_LUAD_clean.h5ad"
adata = sc.read_h5ad(path_148071)

max_val = adata.X[:100].toarray().max() if hasattr(adata.X, 'toarray') else adata.X[:100].max()
print(f"GSE148071 X max (first 100 cells): {max_val}")

if max_val > 50:
    print("X is raw counts, generating lognorm...")
    temp = adata.copy()
    sc.pp.normalize_total(temp, target_sum=1e4)
    sc.pp.log1p(temp)
    adata.layers['lognorm'] = temp.X.copy()
    del temp
else:
    print("X already normalized, using as lognorm")
    adata.layers['lognorm'] = adata.X.copy()

print(f"layers: {list(adata.layers.keys())}")
adata.write(path_148071)
print("GSE148071 saved")
del adata

# === Fix 2: GSE149655 obs_names dedup ===
path_149655 = r"${DATA_ROOT}/cleaned/GSE149655/GSE149655_LUAD_clean.h5ad"
adata = sc.read_h5ad(path_149655)
print(f"\nGSE149655 obs_names unique: {adata.obs_names.is_unique}")
adata.obs_names_make_unique()
print(f"After fix: {adata.obs_names.is_unique}")
adata.write(path_149655)
print("GSE149655 saved")
del adata


adata = sc.read_h5ad(r"${DATA_ROOT}/cleaned/GSE148071_LUAD_clean.h5ad", backed='r')
print(f"layers: {list(adata.layers.keys())}")
adata.file.close()


import scanpy as sc

paths = {
    'GSE164789': r"${DATA_ROOT}/cleaned/GSE164789_LUAD_clean.h5ad",
    'GSE253013': r"${DATA_ROOT}/cleaned/GSE253013_clean.h5ad",
    'GSE131907': r"${DATA_ROOT}/cleaned/GSE131907_clean.h5ad",
    'GSE189357': r"${DATA_ROOT}/cleaned/GSE189357_LUAD_clean.h5ad",
    'GSE148071': r"${DATA_ROOT}/cleaned/GSE148071_LUAD_clean.h5ad",
    'GSE143423': r"${DATA_ROOT}/cleaned/GSE143423_LUAD_clean.h5ad",
}

gene_sets = {}
for name, path in paths.items():
    adata = sc.read_h5ad(path, backed='r')
    gene_sets[name] = set(adata.var_names)
    print(f"{name}: {len(gene_sets[name])} genes")
    adata.file.close()

# Intersection
common = set.intersection(*gene_sets.values())
print(f"\n6-dataset intersection: {len(common)} genes")

# Drop one at a time to find the bottleneck
for exclude in gene_sets:
    rest = {k: v for k, v in gene_sets.items() if k != exclude}
    n = len(set.intersection(*rest.values()))
    print(f"After excluding {exclude}: {n} genes")
