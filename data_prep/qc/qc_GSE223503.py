#!/usr/bin/env python3
# Flattened code cells from qc_GSE223503.ipynb

import scanpy as sc
import tarfile
import os

tar_path = r"${DATA_ROOT}/GSE223503/GSE223503_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE223503/temp_check"
os.makedirs(extract_dir, exist_ok=True)

with tarfile.open(tar_path, 'r') as tar:
    h5_files = [m for m in tar.getmembers() if m.name.endswith('.h5')]
    print(f"h5 file count: {len(h5_files)}")
    tar.extract(h5_files[0], extract_dir)

h5_path = os.path.join(extract_dir, h5_files[0].name)
adata = sc.read_10x_h5(h5_path)

print(f"Shape: {adata.shape}")
print(f"X dtype: {adata.X.dtype}")
print(f"first 5 genes: {adata.var_names[:5].tolist()}")
print(f"MT gene count: {adata.var_names.str.startswith('MT-').sum()}")
print(f"var cols: {adata.var.columns.tolist()}")
print(f"is integer: {(adata.X.data == adata.X.data.astype(int)).all()}")
print(f"mean genes per cell: {(adata.X > 0).sum(1).mean():.0f}")
print(f"mean counts per cell: {adata.X.sum(1).mean():.0f}")


import scanpy as sc
import tarfile
import os

tar_path = r"${DATA_ROOT}/GSE223503/GSE223503_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE223503/temp_check"

# Extract one sample's mtx triplet
with tarfile.open(tar_path, 'r') as tar:
    targets = [m for m in tar.getmembers()
               if 'PA001' in m.name and ('barcodes' in m.name or 'features' in m.name or 'matrix.mtx' in m.name)]
    for t in targets:
        print(f"extract: {t.name}")
        tar.extract(t, extract_dir)

# Read mtx
adata = sc.read_mtx(os.path.join(extract_dir, [t.name for t in targets if 'matrix' in t.name][0])).T

import pandas as pd


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

barcodes = pd.read_csv(os.path.join(extract_dir, [t.name for t in targets if 'barcodes' in t.name][0]),
                       header=None, compression='gzip')
features = pd.read_csv(os.path.join(extract_dir, [t.name for t in targets if 'features' in t.name][0]),
                       header=None, sep='\t', compression='gzip')

print(f"MTX shape: {adata.shape}")
print(f"Barcodes: {barcodes.shape}, first 3: {barcodes[0][:3].tolist()}")
print(f"Features: {features.shape}, ncols: {features.shape[1]}")
print(f"Features head:\n{features.head()}")
print(f"is integer: {(adata.X.data == adata.X.data.astype(int)).all()}")
print(f"mean genes per cell: {(adata.X > 0).sum(1).mean():.0f}")
print(f"mean counts per cell: {adata.X.sum(1).mean():.0f}")


import scanpy as sc
import pandas as pd
import numpy as np
import tarfile
import os
import scrublet as scr
import warnings
warnings.filterwarnings('ignore')

# 1. Extract all mtx triplets from tar (skip h5 and csv)
tar_path = r"${DATA_ROOT}/GSE223503/GSE223503_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE223503/extracted"
os.makedirs(extract_dir, exist_ok=True)

with tarfile.open(tar_path, 'r') as tar:
    members = [m for m in tar.getmembers()
               if ('barcodes' in m.name or 'features' in m.name or 'matrix.mtx' in m.name)]
    print(f"Extracting {len(members)} files...")
    _safe_extract(tar, extract_dir, members=members)

# Sample list
all_files = os.listdir(extract_dir)
sample_ids = sorted(set('_'.join(f.split('_')[1:-1]).replace('_sn', '') for f in all_files if 'matrix' in f))
print(f"Samples: {len(sample_ids)}")
for s in sample_ids:
    print(f"  {s}")
# 2. Build Paper ID -> filename -> metadata mapping
supp_table = pd.read_excel(r"${DATA_ROOT}/GSE223503/41591_2025_3530_MOESM3_ESM.xlsx", skiprows=2)
# Keep only rows with snRNA-seq
supp_snrna = supp_table[supp_table['snRNA-seq chemistry'].notna()].copy()

# Normalise Paper ID (strip dots etc.)
supp_snrna['paper_id_clean'] = supp_snrna['Paper ID'].str.replace('.', 'dot', regex=False).str.strip()
print(f"Supp table snRNA samples: {len(supp_snrna)}")
print(supp_snrna[['Paper ID', 'Sample type(s)', 'Other mutation status ']].to_string())


# 3. Build Paper ID -> filename map
def filename_to_paper_id(fname):
    """NSCLC_PA001 -> PA001, NSCLC_STK_5dot1 -> STK_5.1"""
    pid = fname.replace('NSCLC_', '')
    pid = pid.replace('dot', '.')
    return pid

file_to_paper = {s: filename_to_paper_id(s) for s in sample_ids}
paper_to_file = {v: k for k, v in file_to_paper.items()}

# Verify mapping
supp_snrna['file_name'] = supp_snrna['paper_id_clean'].map(
    {filename_to_paper_id(s).replace('.', 'dot'): s for s in sample_ids}
)
mapped = 0
for pid in supp_snrna['Paper ID']:
    pid_clean = pid.strip().replace('.', 'dot')
    match = [s for s in sample_ids if pid_clean in s or pid.strip().replace('.', 'dot') in s]
    if match:
        mapped += 1
    else:
        print(f"unmatched: {pid}")
print(f"Matched: {mapped}/{len(supp_snrna)}")


# 4. Read all 42 samples + QC + doublet removal
from scipy.io import mmread
import gc

# Build metadata dict from supp table
meta_dict = {}
for _, row in supp_snrna.iterrows():
    pid = row['Paper ID'].strip()
    tissue = 'Brain_met' if 'BRAIN' in str(row['Sample type(s)']) else 'Tumor'
    meta_dict[pid] = {
        'tissue_type': tissue,
        'age': row['Age at resection of profiled specimen'],
        'sex': row['Sex'],
        'mutation': str(row['Other mutation status ']).strip() if pd.notna(row['Other mutation status ']) else 'Unknown',
        'stk11_status': str(row['STK mutation status ']).strip() if pd.notna(row['STK mutation status ']) else 'WT',
        'treatment': str(row['Subsequent treatment received']).strip() if pd.notna(row['Subsequent treatment received']) else 'Unknown',
    }

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    lower = median - nmads * mad * 1.4826
    upper = median + nmads * mad * 1.4826
    return lower, upper

results = []
adatas = []

for sample_name in sample_ids:
    paper_id = filename_to_paper_id(sample_name)
    print(f"\n{'='*50}")
    print(f"Processing: {sample_name} (Paper ID: {paper_id})")

    # --- Read mtx triplet ---
    prefix = [f for f in os.listdir(extract_dir) if sample_name in f and 'matrix' in f][0]
    gsm = prefix.split('_')[0]

    mtx_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_matrix.mtx.gz")
    bar_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_barcodes.tsv.gz")
    feat_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_features.tsv.gz")

    adata = sc.read_mtx(mtx_file).T
    barcodes = pd.read_csv(bar_file, header=None, compression='gzip')
    features = pd.read_csv(feat_file, header=None, sep='\t', compression='gzip')

    adata.obs_names = barcodes[0].values
    adata.var_names = features[1].values  # gene symbol
    adata.var['gene_ids'] = features[0].values  # Ensembl ID
    adata.var['feature_types'] = features[2].values
    # --- Make var names unique ---
    adata.var_names_make_unique()
    n_before = adata.n_obs

    # --- QC metrics (snRNA: do not filter by MT) ---
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, inplace=True)

    # --- Gene filter ---
    sc.pp.filter_genes(adata, min_cells=3)

    # --- Cell filter (MAD adaptive) ---
    n_genes_lower, n_genes_upper = mad_filter(adata.obs['n_genes_by_counts'], nmads=5)
    n_counts_lower, n_counts_upper = mad_filter(adata.obs['total_counts'], nmads=5)

    # snRNA: relax lower bounds
    n_genes_lower = max(n_genes_lower, 150)
    n_counts_lower = max(n_counts_lower, 300)

    keep = (
        (adata.obs['n_genes_by_counts'] > n_genes_lower) &
        (adata.obs['n_genes_by_counts'] < n_genes_upper) &
        (adata.obs['total_counts'] > n_counts_lower) &
        (adata.obs['total_counts'] < n_counts_upper) &
        (adata.obs['pct_counts_hb'] < 5)
        # snRNA: skip pct_counts_mt
    )
    adata = adata[keep].copy()
    n_after_qc = adata.n_obs

    # --- Doublet detection ---
    try:
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3)
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        n_doublets = predicted_doublets.sum()
        adata = adata[~adata.obs['predicted_doublet']].copy()
    except Exception as e:
        print(f"  Scrublet failed: {e}, skip doublet detection")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        n_doublets = 0

    n_final = adata.n_obs

    # --- Add metadata ---
    adata.obs['sample_id'] = sample_name
    adata.obs['patient_id'] = paper_id
    adata.obs['dataset'] = 'GSE223503'
    adata.obs['seq_type'] = 'snRNA'

    meta = meta_dict.get(paper_id, {})
    adata.obs['tissue_type'] = meta.get('tissue_type', 'Unknown')
    adata.obs['age'] = meta.get('age', np.nan)
    adata.obs['sex'] = meta.get('sex', 'Unknown')
    adata.obs['mutation'] = meta.get('mutation', 'Unknown')
    adata.obs['stk11_status'] = meta.get('stk11_status', 'WT')
    adata.obs['treatment'] = meta.get('treatment', 'Unknown')

    # --- Save raw counts ---
    adata.layers['counts'] = adata.X.copy()

    adatas.append(adata)

    result = {
        'sample': sample_name, 'paper_id': paper_id,
        'before': n_before, 'after_qc': n_after_qc,
        'doublets': n_doublets, 'final': n_final,
        'genes_lower': f"{n_genes_lower:.0f}", 'genes_upper': f"{n_genes_upper:.0f}",
    }
    results.append(result)
    print(f"  {n_before} -> QC:{n_after_qc} -> -doublet:{n_final} (doublets:{n_doublets})")

    del adata
    gc.collect()

# --- Summary ---
df_results = pd.DataFrame(results)
print(f"\n{'='*50}")
print(f"Summary:")
print(df_results.to_string(index=False))
print(f"\nTotal cells: {df_results['final'].sum()}")


# 4b. Resume processing from checkpoint
import tarfile

processed = {a.obs['sample_id'].iloc[0] for a in adatas}
print(f"Processed: {len(processed)} samples; continuing with the rest...\n")

for sample_name in sample_ids:
    if sample_name in processed:
        continue

    paper_id = filename_to_paper_id(sample_name)
    print(f"\n{'='*50}")
    print(f"Processing: {sample_name} (Paper ID: {paper_id})")

    # --- Read ---
    prefix = [f for f in os.listdir(extract_dir) if sample_name in f and 'matrix' in f][0]
    gsm = prefix.split('_')[0]

    mtx_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_matrix.mtx.gz")
    bar_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_barcodes.tsv.gz")
    feat_file = os.path.join(extract_dir, f"{gsm}_{sample_name}_sn_features.tsv.gz")

    barcodes = pd.read_csv(bar_file, header=None, compression='gzip')
    features = pd.read_csv(feat_file, header=None, sep='\t', compression='gzip')

    adata = sc.read_mtx(mtx_file).T
    adata.var_names = features[1].values
    adata.var['gene_ids'] = features[0].values
    adata.var['feature_types'] = features[2].values

    if adata.n_obs == len(barcodes):
        adata.obs_names = barcodes[0].values
    else:
        # mtx unfiltered; fall back to h5 + barcodes filter
        print(f"  MTX unfiltered ({adata.n_obs}); switching to h5 + barcodes filter ({len(barcodes)})")
        del adata
        gc.collect()
        h5_name = f"{gsm}_{sample_name}_sn_raw_feature_bc_matrix.h5"
        h5_path = os.path.join(extract_dir, h5_name)
        if not os.path.exists(h5_path):
            with tarfile.open(tar_path, 'r') as tar:
                member = [m for m in tar.getmembers() if m.name == h5_name][0]
                tar.extract(member, extract_dir)
        adata = sc.read_10x_h5(h5_path)
        keep_bc = set(barcodes[0].values)
        adata = adata[adata.obs_names.isin(keep_bc)].copy()
        print(f"  After h5 filter: {adata.n_obs} cells")

    # --- Make var names unique ---
    adata.var_names_make_unique()
    n_before = adata.n_obs

    # --- QC metrics ---
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, inplace=True)

    # --- Gene filter ---
    sc.pp.filter_genes(adata, min_cells=3)

    # --- Cell filter (MAD adaptive) ---
    n_genes_lower, n_genes_upper = mad_filter(adata.obs['n_genes_by_counts'], nmads=5)
    n_counts_lower, n_counts_upper = mad_filter(adata.obs['total_counts'], nmads=5)
    n_genes_lower = max(n_genes_lower, 150)
    n_counts_lower = max(n_counts_lower, 300)

    keep = (
        (adata.obs['n_genes_by_counts'] > n_genes_lower) &
        (adata.obs['n_genes_by_counts'] < n_genes_upper) &
        (adata.obs['total_counts'] > n_counts_lower) &
        (adata.obs['total_counts'] < n_counts_upper) &
        (adata.obs['pct_counts_hb'] < 5)
    )
    adata = adata[keep].copy()
    n_after_qc = adata.n_obs

    # --- Doublet detection ---
    try:
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3)
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        n_doublets = predicted_doublets.sum()
        adata = adata[~adata.obs['predicted_doublet']].copy()
    except Exception as e:
        print(f"  Scrublet failed: {e}")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        n_doublets = 0

    n_final = adata.n_obs

    # --- metadata ---
    adata.obs['sample_id'] = sample_name
    adata.obs['patient_id'] = paper_id
    adata.obs['dataset'] = 'GSE223503'
    adata.obs['seq_type'] = 'snRNA'

    meta = meta_dict.get(paper_id, {})
    adata.obs['tissue_type'] = meta.get('tissue_type', 'Unknown')
    adata.obs['age'] = meta.get('age', np.nan)
    adata.obs['sex'] = meta.get('sex', 'Unknown')
    adata.obs['mutation'] = meta.get('mutation', 'Unknown')
    adata.obs['stk11_status'] = meta.get('stk11_status', 'WT')
    adata.obs['treatment'] = meta.get('treatment', 'Unknown')

    adata.layers['counts'] = adata.X.copy()
    adatas.append(adata)

    results.append({
        'sample': sample_name, 'paper_id': paper_id,
        'before': n_before, 'after_qc': n_after_qc,
        'doublets': n_doublets, 'final': n_final,
        'genes_lower': f"{n_genes_lower:.0f}", 'genes_upper': f"{n_genes_upper:.0f}",
    })
    print(f"  {n_before} -> QC:{n_after_qc} -> -doublet:{n_final} (doublets:{n_doublets})")

    del adata
    gc.collect()

# --- Summary ---
df_results = pd.DataFrame(results)
print(f"\n{'='*50}")
print(f"Summary:")
print(df_results.to_string(index=False))
print(f"\nTotal cells: {df_results['final'].sum()}")


# Repair PA072: manual empty-droplet filter
pa072_prefix = [f for f in os.listdir(extract_dir) if 'PA072' in f and 'matrix' in f][0]
gsm_072 = pa072_prefix.split('_')[0]
h5_name = f"{gsm_072}_NSCLC_PA072_sn_raw_feature_bc_matrix.h5"
h5_path = os.path.join(extract_dir, h5_name)

if not os.path.exists(h5_path):
    with tarfile.open(tar_path, 'r') as tar:
        member = [m for m in tar.getmembers() if 'PA072' in m.name and m.name.endswith('.h5')][0]
        tar.extract(member, extract_dir)

adata = sc.read_10x_h5(h5_path)
print(f"PA072 raw: {adata.shape}")

# Look at gene-count distribution to set empty-droplet threshold
gene_counts = np.array((adata.X > 0).sum(axis=1)).flatten()
print(f"gene-count quantiles: p50={np.median(gene_counts):.0f}, p90={np.percentile(gene_counts, 90):.0f}, p99={np.percentile(gene_counts, 99):.0f}")
print(f"gene>200: {(gene_counts > 200).sum()}")
print(f"gene>500: {(gene_counts > 500).sum()}")
print(f"gene>1000: {(gene_counts > 1000).sum()}")


for i, a in enumerate(adatas):
    if a.n_obs == 0 or (a.n_obs > 0 and a.obs['sample_id'].iloc[0] == 'NSCLC_PA072'):
        adatas[i] = adata
        print(f"Replaced index {i} (PA072)")
        break
else:
    adatas.append(adata)
    print("Appended PA072")


# Merge + save
import anndata as ad

adata_merged = ad.concat(adatas, join='inner')
adata_merged.var_names_make_unique()
adata_merged.obs_names_make_unique()

print(f"Merged: {adata_merged.shape}")
print(f"\ntissue_type:\n{adata_merged.obs['tissue_type'].value_counts()}")
print(f"\npatients: {adata_merged.obs['patient_id'].nunique()}")
print(f"samples: {adata_merged.obs['sample_id'].nunique()}")

x = adata_merged.layers['counts']
print(f"\ncounts integer: {(x.data == x.data.astype(int)).all()}")
print(f"layers: {list(adata_merged.layers.keys())}")

out_path = r"${DATA_ROOT}/GSE223503/GSE223503_LUAD_clean.h5ad"
adata_merged.write(out_path)
print(f"\nSaved: {out_path}")