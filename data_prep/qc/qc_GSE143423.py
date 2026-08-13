#!/usr/bin/env python3
# Flattened code cells from qc_GSE143423.ipynb

import pandas as pd
import gzip

base = r"${DATA_ROOT}/GSE143423"

# Inspect metadata structure
meta = pd.read_csv(f"{base}/GSE143423_lbm_scRNAseq_metadata.csv.gz")
print("metadata shape:", meta.shape)
print(meta.head(10))
print("\ncolumns:", meta.columns.tolist())

# Peek at counts matrix (first 5 rows / 5 cols)
counts_peek = pd.read_csv(f"{base}/GSE143423_lbm_scRNAseq_gene_expression_counts.csv.gz",
                          index_col=0, nrows=5)
print("\ncounts shape (5 rows):", counts_peek.shape)
print(counts_peek.iloc[:5, :5])


import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import anndata as ad
from scipy.sparse import csr_matrix

base = r"${DATA_ROOT}/GSE143423"

# LUAD lbm subset only, skip TNBC
print("Reading counts matrix...")
counts_df = pd.read_csv(f"{base}/GSE143423_lbm_scRNAseq_gene_expression_counts.csv.gz", index_col=0)
meta = pd.read_csv(f"{base}/GSE143423_lbm_scRNAseq_metadata.csv.gz")

print(f"Counts: {counts_df.shape[0]} genes x {counts_df.shape[1]} cells")
print(f"Metadata: {meta.shape[0]} cells")
print(f"Samples: {meta['samples'].value_counts().to_dict()}")

# Transpose: genes x cells -> cells x genes
adata = ad.AnnData(
    X=csr_matrix(counts_df.values.T),
    obs=pd.DataFrame(index=counts_df.columns),
    var=pd.DataFrame(index=counts_df.index)
)
print(f"AnnData: {adata.shape}")

meta = meta.set_index('cell.id')
meta = meta.loc[adata.obs_names]

# Sample -> patient mapping (3 NSCLC patients, 1 brain-met sample each)
adata.obs['sample_id'] = meta['samples'].values
adata.obs['patient_id'] = adata.obs['sample_id'].map({
    'lbm1': 'GSE143423_P1',
    'lbm2': 'GSE143423_P2',
    'lbm3': 'GSE143423_P3'
})
adata.obs['dataset'] = 'GSE143423'
adata.obs['tissue_type'] = 'Brain_met'
adata.obs['stage'] = 'IV'           # brain met = stage IV
adata.obs['chemotherapy'] = 'No'    # treatment-naive

print("\nSample distribution:")
print(adata.obs['sample_id'].value_counts())

sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                           percent_top=None, inplace=True)

print(f"\nMT genes found: {adata.var['mt'].sum()}")
print(f"Before QC: {adata.n_obs} cells, {adata.n_vars} genes")

# Pre-filter visualisation
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, groupby='sample_id')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    lower = median - nmads * mad * 1.4826
    upper = median + nmads * mad * 1.4826
    return lower, upper

n_genes_lo, n_genes_hi = mad_filter(adata.obs['n_genes_by_counts'], nmads=5)
n_counts_lo, n_counts_hi = mad_filter(adata.obs['total_counts'], nmads=5)

# Hard lower bounds
n_genes_lo = max(n_genes_lo, 200)
n_counts_lo = max(n_counts_lo, 500)

print(f"\nQC thresholds:")
print(f"  n_genes: [{n_genes_lo:.0f}, {n_genes_hi:.0f}]")
print(f"  n_counts: [{n_counts_lo:.0f}, {n_counts_hi:.0f}]")
print(f"  pct_mt: < 25%")

cell_mask = (
    (adata.obs['n_genes_by_counts'] >= n_genes_lo) &
    (adata.obs['n_genes_by_counts'] <= n_genes_hi) &
    (adata.obs['total_counts'] >= n_counts_lo) &
    (adata.obs['total_counts'] <= n_counts_hi) &
    (adata.obs['pct_counts_mt'] < 25)
)

n_before = adata.n_obs
adata = adata[cell_mask].copy()
print(f"After QC: {adata.n_obs} cells ({n_before - adata.n_obs} removed, {(n_before-adata.n_obs)/n_before*100:.1f}%)")

adata.obs['doublet_score'] = np.nan
adata.obs['predicted_doublet'] = False

for sample in adata.obs['sample_id'].unique():
    idx = adata.obs['sample_id'] == sample
    adata_sub = adata[idx].copy()
    print(f"\nScrublet on {sample}: {adata_sub.n_obs} cells")

    scrub = scr.Scrublet(adata_sub.X, expected_doublet_rate=0.06)
    scores, preds = scrub.scrub_doublets()

    adata.obs.loc[idx, 'doublet_score'] = scores
    adata.obs.loc[idx, 'predicted_doublet'] = preds
    print(f"  Doublets detected: {preds.sum()} ({preds.sum()/len(preds)*100:.1f}%)")

n_before_dbl = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"\nAfter doublet removal: {adata.n_obs} cells ({n_before_dbl - adata.n_obs} removed)")

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

print("\n" + "="*50)
print("FINAL SUMMARY")
print("="*50)
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"\nSample distribution:")
print(adata.obs['sample_id'].value_counts())
print(f"\nobs columns: {adata.obs.columns.tolist()}")
print(f"\nlayers: {list(adata.layers.keys())}")
print(f"\nMetadata preview:")
print(adata.obs[['sample_id','patient_id','dataset','tissue_type','stage','chemotherapy']].drop_duplicates())

out_path = f"{base}/GSE143423_LUAD_clean.h5ad"
adata.write(out_path)
print(f"\nSaved to: {out_path}")


import scanpy as sc
import numpy as np

paths = {
    'GSE131907': r"${DATA_ROOT}/GSE131907/GSE131907_clean.h5ad",
    'GSE143423': r"${DATA_ROOT}/GSE143423/GSE143423_LUAD_clean.h5ad",
}

required_obs = ['sample_id','patient_id','tissue_type','dataset','stage','chemotherapy']

for name, path in paths.items():
    adata = sc.read_h5ad(path)
    print(f"\n{'='*50}")
    print(f"{name}: {adata.shape}")
    print(f"layers: {list(adata.layers.keys())}")

    if hasattr(adata.X, 'data'):
        x_int = np.allclose(adata.X.data[:1000] % 1, 0)
    else:
        x_int = np.allclose(adata.X[:100] % 1, 0)
    print(f"X is integer: {x_int}  X max: {adata.X.max():.4f}")

    if 'counts' in adata.layers:
        d = adata.layers['counts']
        c_int = np.allclose(d.data[:1000] % 1, 0) if hasattr(d,'data') else np.allclose(d[:100] % 1, 0)
        print(f"counts layer is integer: {c_int}")

    missing = [c for c in required_obs if c not in adata.obs.columns]
    print(f"missing obs columns: {missing if missing else 'none'}")
    print(f"obs columns: {list(adata.obs.columns)}")

    if 'tissue_type' in adata.obs.columns:
        print(f"tissue_type:\n{adata.obs['tissue_type'].value_counts()}")
    if 'dataset' in adata.obs.columns:
        print(f"dataset: {adata.obs['dataset'].unique()}")
