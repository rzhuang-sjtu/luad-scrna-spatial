#!/usr/bin/env python3
# Flattened code cells from qc_GSE308103.ipynb

import tarfile
import os

tar_path = r"${DATA_ROOT}/GSE308103/GSE308103_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE308103/RAW"
os.makedirs(extract_dir, exist_ok=True)

# List filenames first, do not extract yet
with tarfile.open(tar_path, 'r') as tar:
    names = tar.getnames()
    print(f"Total files: {len(names)}")
    for n in names[:20]:
        print(n)
    print("...")
    from collections import Counter
    exts = Counter(os.path.splitext(n)[1] for n in names)
    print(f"\nExtension counts: {exts}")


import gzip
import os

# Extract tar
import tarfile
tar_path = r"${DATA_ROOT}/GSE308103/GSE308103_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE308103/RAW"
os.makedirs(extract_dir, exist_ok=True)

with tarfile.open(tar_path, 'r') as tar:
    _safe_extract(tar, extract_dir)

# Peek at first 10 lines of one file
f = os.path.join(extract_dir, "GSM9237901_P3_Normal.raw_counts.mtx.txt.gz")
with gzip.open(f, 'rt') as fh:
    for i, line in enumerate(fh):
        print(repr(line[:200]))
        if i >= 10:
            break


import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import scrublet as scr
import os
import gzip
import re
from scipy import sparse


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


# Part 1: read all samples + build metadata
raw_dir = r"${DATA_ROOT}/GSE308103/RAW"
files = sorted([f for f in os.listdir(raw_dir) if f.endswith('.gz')])
print(f"{len(files)} files")

def parse_filename(fname):
    """Parse metadata from name e.g. GSM9237901_P3_Normal.raw_counts.mtx.txt.gz"""
    base = fname.replace('.raw_counts.mtx.txt.gz', '')
    parts = base.split('_')
    gsm = parts[0]
    patient = parts[1]  # P3, P4, ...
    tissue_raw = '_'.join(parts[2:])  # Normal, AAH, AIS, MIA, LUAD, Normal1, LUAD1...

    # Strip trailing digits to get tissue_type, keep raw form for sample_id
    tissue_type = re.sub(r'\d+$', '', tissue_raw)  # Normal1 -> Normal
    sample_id = f"{patient}_{tissue_raw}"  # P3_Normal, P4_Normal1
    patient_id = patient

    return {
        'gsm': gsm,
        'sample_id': sample_id,
        'patient_id': patient_id,
        'tissue_type': tissue_type,
        'dataset': 'GSE308103'
    }

# Test parser
for f in files[:5]:
    print(parse_filename(f))


import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import os
import gzip
import re
from scipy import sparse
from multiprocessing import Pool, cpu_count
import time

raw_dir = r"${DATA_ROOT}/GSE308103/RAW"
files = sorted([f for f in os.listdir(raw_dir) if f.endswith('.gz')])

def parse_filename(fname):
    base = fname.replace('.raw_counts.mtx.txt.gz', '')
    parts = base.split('_')
    gsm, patient = parts[0], parts[1]
    tissue_raw = '_'.join(parts[2:])
    tissue_type = re.sub(r'\d+$', '', tissue_raw)
    return {
        'gsm': gsm, 'sample_id': f"{patient}_{tissue_raw}",
        'patient_id': patient, 'tissue_type': tissue_type, 'dataset': 'GSE308103'
    }

def read_one_file(args):
    """Read a single file (multiprocessing worker)"""
    idx, fname, raw_dir = args
    meta = parse_filename(fname)
    filepath = os.path.join(raw_dir, fname)

    t0 = time.time()
    with gzip.open(filepath, 'rt') as fh:
        barcodes = fh.readline().strip().split('\t')
        genes = []
        rows = []
        for line in fh:
            parts = line.strip().split('\t')
            genes.append(parts[0])
            rows.append(np.array(parts[1:], dtype=np.int32))

    mat = sparse.csr_matrix(np.vstack(rows).T)  # cell x gene
    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes)
    )
    for k, v in meta.items():
        adata.obs[k] = v

    elapsed = time.time() - t0
    print(f"[{idx+1}/75] {meta['sample_id']}: {adata.n_obs} cells, {elapsed:.0f}s", flush=True)
    return adata

# Multiprocessing read (8 workers, avoid OOM)
if __name__ == '__main__':
    args_list = [(i, f, raw_dir) for i, f in enumerate(files)]

    n_workers = min(8, cpu_count())
    print(f"Using {n_workers} worker processes...")

    with Pool(n_workers) as pool:
        adatas = pool.map(read_one_file, args_list)

    print("\nMerging...")
    adata = ad.concat(adatas, join='inner', merge='same')
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    print(f"Merged: {adata.n_obs} cells, {adata.n_vars} genes")
    print(adata.obs['tissue_type'].value_counts())
    print(f"Patients: {adata.obs['patient_id'].nunique()}")

    out_path = r"${DATA_ROOT}/GSE308103/GSE308103_raw_merged.h5ad"
    adata.write_h5ad(out_path)
    print(f"Saved: {out_path}")


import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import matplotlib.pyplot as plt

# Part 3: QC
adata = sc.read_h5ad(r"${DATA_ROOT}/GSE308103/GSE308103_raw_merged.h5ad")
adata.obs_names_make_unique()
print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

sc.pp.filter_genes(adata, min_cells=3)
print(f"After gene filter: {adata.n_vars} genes")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, inplace=True)

# snRNA-seq MT distribution (FFPE nuclear RNA: MT should be very low)
print(f"\nMT% distribution:")
print(adata.obs['pct_counts_mt'].describe())
print(f"\nBasic QC distribution:")
print(adata.obs[['n_genes_by_counts','total_counts']].describe())


# Part 3 cont.: MAD adaptive QC + doublet removal

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return median - nmads * mad * 1.4826, median + nmads * mad * 1.4826

# MAD filter per sample
qc_results = []
keep_mask = np.zeros(adata.n_obs, dtype=bool)

for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata.obs.loc[idx]
    n_before = idx.sum()

    ng_lo, ng_hi = mad_filter(sub['n_genes_by_counts'], nmads=5)
    nc_lo, nc_hi = mad_filter(sub['total_counts'], nmads=5)

    # Hard lower bounds
    ng_lo = max(ng_lo, 200)
    nc_lo = max(nc_lo, 300)

    # snRNA-seq: tight MT (< 5%) - nuclear RNA has near-zero mito; > 5% indicates lysed cells
    mask = (
        (sub['n_genes_by_counts'] >= ng_lo) & (sub['n_genes_by_counts'] <= ng_hi) &
        (sub['total_counts'] >= nc_lo) & (sub['total_counts'] <= nc_hi) &
        (sub['pct_counts_mt'] < 5) &
        (sub['pct_counts_hb'] < 5)
    )
    keep_mask[np.where(idx)[0]] = mask.values
    n_after = mask.sum()

    qc_results.append({
        'sample_id': sample, 'tissue_type': sub['tissue_type'].iloc[0],
        'before': n_before, 'after': n_after,
        'removed_pct': f"{(1 - n_after/n_before)*100:.1f}%"
    })

qc_df = pd.DataFrame(qc_results)
print(f"Pre-QC: {adata.n_obs} cells")
adata = adata[keep_mask].copy()
print(f"Post-QC: {adata.n_obs} cells")
print(f"\nRemaining tissue types:")
print(adata.obs['tissue_type'].value_counts())
print(f"\nQC summary:")
print(qc_df.to_string(index=False))

# Part 4: Scrublet doublet removal (per sample)
print("\n\n========== Doublet Detection ==========")
doublet_results = []

for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata[idx].copy()
    n_cells = sub.n_obs

    if n_cells < 100:
        print(f"[SKIP] {sample}: {n_cells} cells, too few")
        adata.obs.loc[idx, 'doublet_score'] = 0
        adata.obs.loc[idx, 'predicted_doublet'] = False
        continue

    try:
        scrub = scr.Scrublet(sub.X, expected_doublet_rate=0.06)
        scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                              min_gene_variability_pctl=85,
                                              n_prin_comps=min(30, n_cells-1))
        adata.obs.loc[idx, 'doublet_score'] = scores
        adata.obs.loc[idx, 'predicted_doublet'] = preds
        n_db = preds.sum()
        print(f"{sample}: {n_cells} cells, {n_db} doublets ({n_db/n_cells*100:.1f}%)")
        doublet_results.append({'sample_id': sample, 'cells': n_cells,
                                'doublets': n_db, 'rate': f"{n_db/n_cells*100:.1f}%"})
    except Exception as e:
        print(f"[ERROR] {sample}: {e}")
        adata.obs.loc[idx, 'doublet_score'] = 0
        adata.obs.loc[idx, 'predicted_doublet'] = False

n_before_db = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"\nDoublet removal: {n_before_db} -> {adata.n_obs} cells")
print("Tissue types:")
print(adata.obs['tissue_type'].value_counts())

out_path = r"${DATA_ROOT}/GSE308103/GSE308103_qc_doublet.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved: {out_path}")


# Drop doublets
n_before = adata.n_obs
adata = adata[adata.obs['predicted_doublet'] == False].copy()
print(f"Doublet removal: {n_before} -> {adata.n_obs}, removed {n_before - adata.n_obs}")
print(adata.obs['tissue_type'].value_counts())

# Save raw counts + normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Stage info from tissue_type
stage_map = {
    'Normal': 'Normal',
    'AAH': 'Preinvasive',
    'AIS': 'Preinvasive',
    'MIA': 'Preinvasive',
    'LUAD': 'Invasive'
}
adata.obs['stage'] = adata.obs['tissue_type'].map(stage_map)
adata.obs['chemotherapy'] = 'No'  # treatment-naive

print(f"\n{'='*50}")
print(f"FINAL: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Tissue types:\n{adata.obs['tissue_type'].value_counts()}")
print(f"obs cols: {list(adata.obs.columns)}")
print(f"layers: {list(adata.layers.keys())}")

out_path = r"${DATA_ROOT}/GSE308103/GSE308103_LUAD_clean.h5ad"
adata.write_h5ad(out_path)
print(f"Saved: {out_path}")


import scanpy as sc
import numpy as np
import scrublet as scr

adata = sc.read_h5ad(r"${DATA_ROOT}/GSE308103/GSE308103_raw_merged.h5ad")
adata.obs_names_make_unique()
print(f"Loaded: {adata.n_obs} cells")

# ---- QC ----
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, inplace=True)

def mad_filter(s, nmads=5):
    med = np.median(s)
    mad = np.median(np.abs(s - med))
    return med - nmads*mad*1.4826, med + nmads*mad*1.4826

lo_g, hi_g = mad_filter(adata.obs['n_genes_by_counts'])
lo_c, hi_c = mad_filter(adata.obs['total_counts'])
keep = ((adata.obs['n_genes_by_counts'] > max(200, lo_g)) &
        (adata.obs['n_genes_by_counts'] < hi_g) &
        (adata.obs['total_counts'] > max(500, lo_c)) &
        (adata.obs['total_counts'] < hi_c) &
        (adata.obs['pct_counts_mt'] < 25) &
        (adata.obs['pct_counts_hb'] < 5))
print(f"QC: {adata.n_obs} -> {keep.sum()}")
adata = adata[keep].copy()

# ---- Doublet ----
adata.obs['doublet_score'] = 0.0
adata.obs['predicted_doublet'] = False
for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata[idx].copy()
    print(f"Scrublet {sample}: {sub.n_obs} cells...", end=' ', flush=True)
    scrub = scr.Scrublet(sub.X, expected_doublet_rate=0.06)
    scores, preds = scrub.scrub_doublets(verbose=False)
    adata.obs.loc[idx, 'doublet_score'] = scores
    adata.obs.loc[idx, 'predicted_doublet'] = preds
    print(f"{preds.sum()} doublets ({preds.sum()/sub.n_obs*100:.1f}%)")

n_before = adata.n_obs
adata = adata[adata.obs['predicted_doublet'] == False].copy()
print(f"Doublet removal: {n_before} -> {adata.n_obs}")

# ---- Normalize + Save ----
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

adata.obs['stage'] = adata.obs['tissue_type'].map({
    'Normal':'Normal','AAH':'Preinvasive','AIS':'Preinvasive',
    'MIA':'Preinvasive','LUAD':'Invasive'})
adata.obs['chemotherapy'] = 'No'
adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)

print(f"FINAL: {adata.n_obs} cells, {adata.n_vars} genes")
print(adata.obs['tissue_type'].value_counts())

adata.write_h5ad(r"${DATA_ROOT}/GSE308103/GSE308103_LUAD_clean.h5ad")
print("Saved")