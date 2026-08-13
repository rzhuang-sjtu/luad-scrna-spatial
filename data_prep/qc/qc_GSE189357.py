#!/usr/bin/env python3
# Flattened code cells from qc_GSE189357.ipynb

# GSE189357 cleaning
# LUAD: AIS (TD1-3) -> MIA (TD4-6) -> IAC (TD7-9)

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os
import warnings
warnings.filterwarnings('ignore')

# 1. List directory contents (after extraction)
data_dir = r"${DATA_ROOT}/GSE189357"
for f in sorted(os.listdir(data_dir)):
    print(f)


import tarfile
import os

data_dir = r"${DATA_ROOT}/GSE189357"
tar_path = os.path.join(data_dir, "GSE189357_RAW.tar")

# Extract
with tarfile.open(tar_path, 'r') as tar:
    _safe_extract(tar, data_dir)

for f in sorted(os.listdir(data_dir)):
    print(f)


import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os
import shutil
import warnings


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

warnings.filterwarnings('ignore')

data_dir = r"${DATA_ROOT}/GSE189357"
output_dir = r"${DATA_ROOT}/GSE189357/processed"
os.makedirs(output_dir, exist_ok=True)

# 1. Reorganize files into per-sample folders (scanpy needs the standard layout)
samples = {
    'TD1': {'gsm': 'GSM5699777', 'stage': 'AIS',  'patient_id': 'P1'},
    'TD2': {'gsm': 'GSM5699778', 'stage': 'AIS',  'patient_id': 'P2'},
    'TD3': {'gsm': 'GSM5699779', 'stage': 'AIS',  'patient_id': 'P3'},
    'TD4': {'gsm': 'GSM5699780', 'stage': 'MIA',  'patient_id': 'P4'},
    'TD5': {'gsm': 'GSM5699781', 'stage': 'MIA',  'patient_id': 'P5'},
    'TD6': {'gsm': 'GSM5699782', 'stage': 'MIA',  'patient_id': 'P6'},
    'TD7': {'gsm': 'GSM5699783', 'stage': 'IAC',  'patient_id': 'P7'},
    'TD8': {'gsm': 'GSM5699784', 'stage': 'IAC',  'patient_id': 'P8'},
    'TD9': {'gsm': 'GSM5699785', 'stage': 'IAC',  'patient_id': 'P9'},
}

for sample, info in samples.items():
    sample_dir = os.path.join(data_dir, sample)
    os.makedirs(sample_dir, exist_ok=True)
    gsm = info['gsm']
    # Rename to standard 10X filenames
    for src_suffix, dst_name in [
        (f'{gsm}_{sample}_barcodes.tsv.gz', 'barcodes.tsv.gz'),
        (f'{gsm}_{sample}_features.tsv.gz', 'features.tsv.gz'),
        (f'{gsm}_{sample}_matrix.mtx.gz',   'matrix.mtx.gz'),
    ]:
        src = os.path.join(data_dir, src_suffix)
        dst = os.path.join(sample_dir, dst_name)
        if not os.path.exists(dst):
            shutil.copy2(src, dst)

print("File reorg done")


# 2. Read all samples
adatas = {}
for sample, info in samples.items():
    sample_dir = os.path.join(data_dir, sample)
    adata = sc.read_10x_mtx(sample_dir, var_names='gene_symbols', cache=False)
    adata.obs['sample_id'] = sample
    adata.obs['patient_id'] = info['patient_id']
    adata.obs['dataset'] = 'GSE189357'
    adata.obs['stage'] = info['stage']
    adata.obs['tissue_type'] = 'Tumor'  # all resected tumor tissue
    adata.obs['chemotherapy'] = 'treatment_naive'
    adata.var_names_make_unique()
    adatas[sample] = adata
    print(f"{sample} ({info['stage']}): {adata.n_obs} cells, {adata.n_vars} genes")

adata = sc.concat(list(adatas.values()), merge='same')
adata.obs_names_make_unique()
print(f"\nMerged: {adata.n_obs} cells, {adata.n_vars} genes")


# 3. QC metrics
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'],
                           percent_top=None, inplace=True)

print("\n=== Cells per sample (pre-QC) ===")
print(adata.obs['sample_id'].value_counts().sort_index())


# 4. MAD adaptive QC
def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    lower = median - nmads * mad * 1.4826
    upper = median + nmads * mad * 1.4826
    return lower, upper

# Per-sample adaptive thresholds
qc_mask = pd.Series(True, index=adata.obs_names)

print("\n=== Per-sample MAD thresholds ===")
for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    sub = adata.obs.loc[idx]

    ng_lo, ng_hi = mad_filter(sub['n_genes_by_counts'], nmads=5)
    nc_lo, nc_hi = mad_filter(sub['total_counts'], nmads=5)

    # Hard lower bounds
    ng_lo = max(ng_lo, 200)
    nc_lo = max(nc_lo, 500)

    mask = (
        (sub['n_genes_by_counts'] >= ng_lo) &
        (sub['n_genes_by_counts'] <= ng_hi) &
        (sub['total_counts'] >= nc_lo) &
        (sub['total_counts'] <= nc_hi) &
        (sub['pct_counts_mt'] <= 25) &
        (sub['pct_counts_hb'] <= 5)
    )
    qc_mask.loc[idx] = mask.values

    n_before = idx.sum()
    n_after = mask.sum()
    print(f"{sample}: genes[{ng_lo:.0f}-{ng_hi:.0f}], counts[{nc_lo:.0f}-{nc_hi:.0f}], "
          f"kept {n_after}/{n_before} ({n_after/n_before*100:.1f}%)")

adata = adata[qc_mask].copy()
print(f"\nAfter QC: {adata.n_obs} cells")


# 5. Per-sample Scrublet
print("\n=== Doublet detection ===")
adata.obs['doublet_score'] = 0.0
adata.obs['predicted_doublet'] = False

for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    adata_sub = adata[idx].copy()

    scrub = scr.Scrublet(adata_sub.X, expected_doublet_rate=0.06)
    scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                          min_gene_variability_pctl=85)

    adata.obs.loc[idx, 'doublet_score'] = scores
    adata.obs.loc[idx, 'predicted_doublet'] = preds

    n_doublet = preds.sum()
    print(f"{sample}: {n_doublet}/{idx.sum()} doublets ({n_doublet/idx.sum()*100:.1f}%)")

adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"\nAfter doublet removal: {adata.n_obs} cells")


# 6. Save
adata.X = adata.X.copy()
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

save_path = os.path.join(output_dir, "GSE189357_LUAD_cleaned.h5ad")
adata.write(save_path)
print(f"\nSaved: {save_path}")
print(f"Final: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"\nobs cols: {list(adata.obs.columns)}")
print(f"layers: {list(adata.layers.keys())}")

summary = adata.obs.groupby('sample_id').size().reset_index(name='final_cells')
summary['stage'] = summary['sample_id'].map({k: v['stage'] for k, v in samples.items()})
print("\n=== Final cells per sample ===")
print(summary.to_string(index=False))