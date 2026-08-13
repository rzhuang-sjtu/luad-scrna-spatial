#!/usr/bin/env python3
# Flattened code cells from qc_GSE149655.ipynb

import tarfile
import os

# List RAW.tar contents
tar_path = r"${DATA_ROOT}/GSE149655/GSE149655_RAW.tar"
with tarfile.open(tar_path, 'r') as tar:
    for member in tar.getmembers():
        print(member.name)


import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import anndata as ad
import tarfile
import os
import warnings
warnings.filterwarnings('ignore')

# GSE149655 cleaning pipeline

# Step 0: extract RAW.tar
tar_path = r"${DATA_ROOT}/GSE149655/GSE149655_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE149655/RAW"
os.makedirs(extract_dir, exist_ok=True)

with tarfile.open(tar_path, 'r') as tar:
    _safe_extract(tar, extract_dir)

for root, dirs, files in os.walk(extract_dir):
    for f in files:
        print(os.path.join(root, f))


import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import anndata as ad
from scipy import sparse
from scipy.io import mmread
import gzip
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

raw_dir = r"${DATA_ROOT}/GSE149655/RAW"

# Sample table
samples = {
    'GSM4506698': {
        'mtx': 'GSM4506698_Case2_nor_matrix.mtx.gz',
        'barcodes': 'GSM4506698_Case2_nor_barcodes.tsv.gz',
        'features': 'GSM4506698_Case2_nor_features.tsv.gz',
        'sample_id': 'Case2_Normal', 'patient_id': 'Case2', 'tissue_type': 'Normal'
    },
    'GSM4506699': {
        'mtx': 'GSM4506699_Case2_ade1_matrix.mtx.gz',
        'barcodes': 'GSM4506699_Case2_ade1_barcodes.tsv.gz',
        'features': 'GSM4506699_Case2_ade1_features.tsv.gz',
        'sample_id': 'Case2_Tumor', 'patient_id': 'Case2', 'tissue_type': 'Tumor'
    },
    'GSM4506700': {
        'mtx': 'GSM4506700_Case3_nor_matrix.mtx.gz',
        'barcodes': 'GSM4506700_Case3_nor_barcodes.tsv.gz',
        'features': 'GSM4506700_Case3_nor_feature.tsv.gz',  # note: singular 'feature'
        'sample_id': 'Case3_Normal', 'patient_id': 'Case3', 'tissue_type': 'Normal'
    },
    'GSM4506701': {
        'mtx': 'GSM4506701_Case3_ade1_matrix.tsv.gz',  # note: tsv not mtx
        'barcodes': 'GSM4506701_Case3_ade1_barcodes.tsv.gz',
        'features': 'GSM4506701_Case3_ade1_features.tsv.gz',
        'sample_id': 'Case3_Tumor', 'patient_id': 'Case3', 'tissue_type': 'Tumor'
    },
}

# Step 0: peek at GSM4506701's matrix.tsv.gz format
import os
for gsm, info in samples.items():
    mtx_path = os.path.join(raw_dir, info['mtx'])
    with gzip.open(mtx_path, 'rt') as f:
        first_lines = [f.readline() for _ in range(5)]
    print(f"\n{gsm} ({info['mtx']}) first 5 lines:")
    for line in first_lines:
        print(repr(line.strip()))


import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import anndata as ad
from scipy.io import mmread
from scipy import sparse
import gzip
import os
import warnings
warnings.filterwarnings('ignore')

raw_dir = r"${DATA_ROOT}/GSE149655/RAW"

samples = {
    'GSM4506698': {
        'mtx': 'GSM4506698_Case2_nor_matrix.mtx.gz',
        'barcodes': 'GSM4506698_Case2_nor_barcodes.tsv.gz',
        'features': 'GSM4506698_Case2_nor_features.tsv.gz',
        'sample_id': 'Case2_Normal', 'patient_id': 'Case2', 'tissue_type': 'Normal'
    },
    'GSM4506699': {
        'mtx': 'GSM4506699_Case2_ade1_matrix.mtx.gz',
        'barcodes': 'GSM4506699_Case2_ade1_barcodes.tsv.gz',
        'features': 'GSM4506699_Case2_ade1_features.tsv.gz',
        'sample_id': 'Case2_Tumor', 'patient_id': 'Case2', 'tissue_type': 'Tumor'
    },
    'GSM4506700': {
        'mtx': 'GSM4506700_Case3_nor_matrix.mtx.gz',
        'barcodes': 'GSM4506700_Case3_nor_barcodes.tsv.gz',
        'features': 'GSM4506700_Case3_nor_feature.tsv.gz',
        'sample_id': 'Case3_Normal', 'patient_id': 'Case3', 'tissue_type': 'Normal'
    },
    'GSM4506701': {
        'mtx': 'GSM4506701_Case3_ade1_matrix.tsv.gz',
        'barcodes': 'GSM4506701_Case3_ade1_barcodes.tsv.gz',
        'features': 'GSM4506701_Case3_ade1_features.tsv.gz',
        'sample_id': 'Case3_Tumor', 'patient_id': 'Case3', 'tissue_type': 'Tumor'
    },
}

# Step 1: read
adatas = []
for gsm, info in samples.items():
    mtx = mmread(os.path.join(raw_dir, info['mtx'])).T.tocsr()  # genes x cells -> cells x genes
    barcodes = pd.read_csv(os.path.join(raw_dir, info['barcodes']), header=None, sep='\t')
    features = pd.read_csv(os.path.join(raw_dir, info['features']), header=None, sep='\t')

    adata = ad.AnnData(X=mtx)
    adata.obs_names = barcodes[0].values
    adata.var_names = features[1].values  # column 2 = gene symbol
    adata.var['gene_ids'] = features[0].values  # column 1 = Ensembl ID
    adata.var_names_make_unique()

    adata.obs['sample_id'] = info['sample_id']
    adata.obs['patient_id'] = info['patient_id']
    adata.obs['tissue_type'] = info['tissue_type']
    adata.obs['dataset'] = 'GSE149655'
    adata.obs['stage'] = 'Early'
    adata.obs['chemotherapy'] = 'No'
    adata.obs['kras_status'] = 'Mutant'
    adata.obs['cd45_filtered'] = True

    print(f"{gsm} ({info['sample_id']}): {adata.shape[0]} cells x {adata.shape[1]} genes")
    adatas.append(adata)

adata = ad.concat(adatas, join='inner')
adata.var_names_make_unique()
print(f"\nMerged: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Step 2: QC
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                           percent_top=None, inplace=True)

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    lower = median - nmads * mad * 1.4826
    upper = median + nmads * mad * 1.4826
    return lower, upper

n_genes_lower, n_genes_upper = mad_filter(adata.obs['n_genes_by_counts'], nmads=5)
n_counts_lower, n_counts_upper = mad_filter(adata.obs['total_counts'], nmads=5)
n_genes_lower = max(n_genes_lower, 200)
n_counts_lower = max(n_counts_lower, 500)

print(f"n_genes cutoff: [{n_genes_lower:.0f}, {n_genes_upper:.0f}]")
print(f"n_counts cutoff: [{n_counts_lower:.0f}, {n_counts_upper:.0f}]")

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')

pre_filter = adata.n_obs
adata = adata[
    (adata.obs['n_genes_by_counts'] >= n_genes_lower) &
    (adata.obs['n_genes_by_counts'] <= n_genes_upper) &
    (adata.obs['total_counts'] >= n_counts_lower) &
    (adata.obs['total_counts'] <= n_counts_upper) &
    (adata.obs['pct_counts_mt'] < 25) &
    (adata.obs['pct_counts_hb'] < 5)
].copy()
post_filter = adata.n_obs
print(f"QC filter: {pre_filter} -> {post_filter} (removed {pre_filter-post_filter}, {(pre_filter-post_filter)/pre_filter*100:.1f}%)")
print(adata.obs['sample_id'].value_counts())

# Step 3: Scrublet
adata.obs['doublet_score'] = np.nan
adata.obs['predicted_doublet'] = False

for sample in adata.obs['sample_id'].unique():
    mask = adata.obs['sample_id'] == sample
    adata_sub = adata[mask].copy()
    try:
        scrub = scr.Scrublet(adata_sub.X, expected_doublet_rate=0.06)
        scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85)
        adata.obs.loc[mask, 'doublet_score'] = scores
        adata.obs.loc[mask, 'predicted_doublet'] = preds
        print(f"{sample}: {adata_sub.n_obs} cells, {preds.sum()} doublets ({preds.sum()/adata_sub.n_obs*100:.1f}%)")
    except Exception as e:
        print(f"{sample}: Scrublet error - {e}")

pre_d = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"Doublet filter: {pre_d} -> {adata.n_obs} (removed {pre_d - adata.n_obs})")

# Step 4: normalize + save
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Step 5: final check
print("\n" + "="*60)
print("GSE149655 cleaning done")
print("="*60)
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"obs cols: {list(adata.obs.columns)}")
print(f"layers: {list(adata.layers.keys())}")
print(f"\nCells per sample:\n{adata.obs['sample_id'].value_counts()}")
print(f"\nMT genes: {adata.var['mt'].sum()}, Ribo genes: {adata.var['ribo'].sum()}")

# Verify counts integer
test = adata.layers['counts'][:10,:10]
if sparse.issparse(test):
    test = test.toarray()
print(f"counts is integer: {np.allclose(test % 1, 0)}")

output_path = r"${DATA_ROOT}/GSE149655/GSE149655_LUAD_clean.h5ad"
adata.write(output_path)
print(f"Saved: {output_path}")