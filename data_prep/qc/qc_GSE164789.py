#!/usr/bin/env python3
# Flattened code cells from qc_GSE164789.ipynb

import tarfile, os

# List first 50 entries inside RAW.tar
with tarfile.open(r"${DATA_ROOT}/GSE164789/GSE164789_RAW.tar", "r") as tar:
    names = tar.getnames()
    print(f"Total files: {len(names)}")
    for n in names[:50]:
        print(n)


"""
GSE164789 - Precursor LUAD scRNA-seq cleaning pipeline
62 scRNA samples (excludes 16 TCR samples)
"""
import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import anndata as ad
import tarfile
import os
import tempfile
import shutil
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2

# Step 1: extract RAW.tar and read all scRNA samples
raw_tar = r"${DATA_ROOT}/GSE164789/GSE164789_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE164789/raw_files"
os.makedirs(extract_dir, exist_ok=True)

print("Extracting RAW.tar...")
with tarfile.open(raw_tar, "r") as tar:
    _safe_extract(tar, extract_dir)

# Collect scRNA sample names (drop TCR)
all_files = os.listdir(extract_dir)
sample_names = sorted(set(
    f.split('.')[0] for f in all_files if f.endswith('.matrix.mtx.gz')
))
sample_names = [s for s in sample_names if 'TCR' not in s.upper()]
print(f"scRNA samples: {len(sample_names)}")

# Step 2: read & merge per sample
from scipy.io import mmread
import gzip


def _safe_extract(tar, dest, members=None):
    """Extract without letting a member escape dest.

    filter="data" is the one-line answer but only exists from Python 3.12;
    the released environment pins 3.11, so fall back to checking each member.
    """
    import os
    root = os.path.realpath(dest)
    try:
        tar.extractall(dest, members=members, filter="data")
        return
    except TypeError:
        pass
    picked = list(tar.getmembers() if members is None else members)
    for m in picked:
        target = os.path.realpath(os.path.join(root, m.name))
        if target != root and not target.startswith(root + os.sep):
            raise ValueError(f"unsafe path in archive: {m.name}")
        if m.issym() or m.islnk() or m.isdev():
            raise ValueError(f"unsupported member type in archive: {m.name}")
    tar.extractall(dest, members=picked)


adatas = []
for sname in sample_names:
    sample_id = sname.split('_', 1)[1]

    mtx_file = os.path.join(extract_dir, f"{sname}.matrix.mtx.gz")
    bar_file = os.path.join(extract_dir, f"{sname}.barcodes.tsv.gz")
    gene_file = os.path.join(extract_dir, f"{sname}.genes.tsv.gz")

    mat = mmread(mtx_file).T.tocsr()  # cells x genes
    barcodes = pd.read_csv(bar_file, header=None, sep='\t')[0].values
    genes = pd.read_csv(gene_file, header=None, sep='\t')

    adata = sc.AnnData(X=mat)
    adata.obs_names = barcodes.astype(str)

    # Gene names: column 2 is symbol if two columns
    if genes.shape[1] >= 2:
        adata.var_names = genes[1].values.astype(str)
        adata.var['gene_ids'] = genes[0].values
    else:
        adata.var_names = genes[0].values.astype(str)

    adata.var_names_make_unique()

    # Parse metadata from sample id
    parts = sample_id.rsplit('-', 1)
    patient_id = parts[0]
    suffix = parts[1] if len(parts) > 1 else 'Unknown'

    if suffix.startswith('A') or suffix == 'KA':
        tissue_type = 'Normal'
    elif suffix.startswith('T') or suffix == 'KT':
        tissue_type = 'Tumor'
    else:
        tissue_type = 'Unknown'

    adata.obs['sample_id'] = sample_id
    adata.obs['patient_id'] = patient_id
    adata.obs['tissue_type'] = tissue_type
    adata.obs['dataset'] = 'GSE164789'

    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

    adatas.append(adata)
    print(f"  {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes, tissue={tissue_type}")

adata = ad.concat(adatas, join='inner')
adata.obs_names_make_unique()
print(f"\nMerged: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Patients: {adata.obs['patient_id'].nunique()}")
print(adata.obs['tissue_type'].value_counts())

# Step 3: QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                           percent_top=None, inplace=True)

sc.pp.filter_genes(adata, min_cells=3)
print(f"After gene filter: {adata.n_vars} genes")

# Step 4: per-sample MAD QC
def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return median - nmads * mad * 1.4826, median + nmads * mad * 1.4826

qc_records = []
cells_before = adata.n_obs

keep_mask = np.ones(adata.n_obs, dtype=bool)
for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata.obs.loc[idx]
    n_before = idx.sum()

    gl, gu = mad_filter(sub['n_genes_by_counts'], nmads=5)
    cl, cu = mad_filter(sub['total_counts'], nmads=5)

    # Hard lower bounds
    gl = max(gl, 200)
    cl = max(cl, 500)

    mask = (
        (sub['n_genes_by_counts'] >= gl) & (sub['n_genes_by_counts'] <= gu) &
        (sub['total_counts'] >= cl) & (sub['total_counts'] <= cu) &
        (sub['pct_counts_mt'] < 25) &
        (sub['pct_counts_hb'] < 5)
    )
    keep_mask[idx.values] = mask.values

    n_after = mask.sum()
    qc_records.append({
        'sample_id': sample,
        'before': n_before,
        'after': n_after,
        'removed_pct': f"{(1 - n_after/n_before)*100:.1f}%"
    })

adata = adata[keep_mask].copy()
qc_df = pd.DataFrame(qc_records)
print(f"\nQC filter: {cells_before} -> {adata.n_obs} cells")
print(qc_df.to_string(index=False))

# Step 5: per-sample Scrublet
doublet_records = []
for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata[idx].copy()

    if sub.n_obs < 50:
        adata.obs.loc[idx, 'doublet_score'] = 0
        adata.obs.loc[idx, 'predicted_doublet'] = False
        doublet_records.append({'sample': sample, 'n_cells': sub.n_obs,
                                'n_doublets': 0, 'rate': '0.0%'})
        continue

    try:
        scrub = scr.Scrublet(sub.X, expected_doublet_rate=0.06)
        scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                              min_gene_variability_pctl=85)
        adata.obs.loc[idx, 'doublet_score'] = scores
        adata.obs.loc[idx, 'predicted_doublet'] = preds
        n_d = preds.sum()
        doublet_records.append({'sample': sample, 'n_cells': sub.n_obs,
                                'n_doublets': n_d,
                                'rate': f"{n_d/sub.n_obs*100:.1f}%"})
    except Exception as e:
        print(f"  Scrublet failed for {sample}: {e}")
        adata.obs.loc[idx, 'doublet_score'] = 0
        adata.obs.loc[idx, 'predicted_doublet'] = False
        doublet_records.append({'sample': sample, 'n_cells': sub.n_obs,
                                'n_doublets': 0, 'rate': 'failed'})

doublet_df = pd.DataFrame(doublet_records)
print(f"\nDoublet detection results:")
print(doublet_df.to_string(index=False))

n_before_doublet = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"\nAfter doublet removal: {n_before_doublet} -> {adata.n_obs} cells")

# Step 6: save
# Ensure X is raw counts
adata.layers['counts'] = adata.X.copy()

# Log-normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

output_path = r"${DATA_ROOT}/GSE164789/GSE164789_LUAD_cleaned.h5ad"
adata.write_h5ad(output_path)
print(f"\nSaved to: {output_path}")
print(f"Final: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Patients: {adata.obs['patient_id'].nunique()}")
print(f"Samples: {adata.obs['sample_id'].nunique()}")
print(adata.obs['tissue_type'].value_counts())


# Step 5 (alt): Doublet detection (per sample)
import scrublet as scr

doublet_results = []
for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata[idx].copy()

    if sub.n_obs < 50:
        print(f"  {sample}: skip (only {sub.n_obs} cells)")
        adata.obs.loc[idx, 'doublet_score'] = 0.0
        adata.obs.loc[idx, 'predicted_doublet'] = False
        continue

    counts = sub.layers['counts'] if 'counts' in sub.layers else sub.X
    scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
    scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                          min_gene_variability_pctl=85,
                                          n_prin_comps=30)

    adata.obs.loc[idx, 'doublet_score'] = scores
    adata.obs.loc[idx, 'predicted_doublet'] = preds

    n_doublet = preds.sum()
    pct = n_doublet / len(preds) * 100
    doublet_results.append({'sample': sample, 'total': len(preds),
                           'doublets': n_doublet, 'pct': f"{pct:.1f}%"})
    print(f"  {sample}: {len(preds)} cells, {n_doublet} doublets ({pct:.1f}%)")

doublet_df = pd.DataFrame(doublet_results)
print(doublet_df.to_string(index=False))

# Cast to bool to be safe
n_before_doublet = adata.n_obs
mask = adata.obs['predicted_doublet'].astype(bool).fillna(False)
adata = adata[~mask].copy()
print(f"\nDoublet removal: {n_before_doublet} -> {adata.n_obs} cells")

# Step 6: save
if 'counts' not in adata.layers:
    adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

output_path = r"${DATA_ROOT}/GSE164789/GSE164789_LUAD_cleaned.h5ad"
adata.write(output_path)
print(f"\nSaved to: {output_path}")
print(f"Final: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Patients: {adata.obs['patient_id'].nunique()}")
print(f"Samples: {adata.obs['sample_id'].nunique()}")
print(adata.obs['tissue_type'].value_counts())
print(f"\nlayers: {list(adata.layers.keys())}")
print(f"obs columns: {list(adata.obs.columns)}")


# Cast types and save
adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)
adata.obs['doublet_score'] = adata.obs['doublet_score'].astype(float)
for col in adata.obs.columns:
    if adata.obs[col].dtype == object:
        adata.obs[col] = adata.obs[col].astype(str)

output_path = r"${DATA_ROOT}/GSE164789/GSE164789_LUAD_cleaned.h5ad"
adata.write(output_path)
print(f"Saved to: {output_path}")
print(f"Final: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Patients: {adata.obs['patient_id'].nunique()}")
print(f"Samples: {adata.obs['sample_id'].nunique()}")
print(adata.obs['tissue_type'].value_counts())
print(f"layers: {list(adata.layers.keys())}")