#!/usr/bin/env python3
# Flattened code cells from qc_GSE131907.ipynb

import pandas as pd
import gzip

ann_path = r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
ann = pd.read_csv(ann_path, sep='\t', index_col=0)
print(ann.shape)
print(ann.columns.tolist())
print(ann.head(3))


import gzip

mat_path = r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
with gzip.open(mat_path, 'rt') as f:
    header = f.readline()
    row1 = f.readline()

cells = header.strip().split('\t')
print(f"n_columns (cells + 1): {len(cells)}")
print(f"first 5 column names: {cells[:5]}")
print(f"first row gene + 5 values: {row1.strip().split(chr(9))[:6]}")


import pandas as pd, gzip
ann = pd.read_csv(r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep='\t', index_col=0)
print(ann['Sample_Origin'].value_counts())
print(ann['Sample'].value_counts().head(20))


import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scrublet as scr
import gzip

mat_path = r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
ann_path = r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
out_path = r"${DATA_ROOT}/GSE131907/GSE131907_clean.h5ad"

# Annotation
ann = pd.read_csv(ann_path, sep='\t', index_col=0)

origin_map = {
    'nLung':  'Normal_Lung',
    'tLung':  'Primary_Tumor',
    'nLN':    'LN_Normal',
    'mLN':    'LN_Met',
    'mBrain': 'Brain_Met',
    'PE':     'Pleural_Effusion',
    'tL/B':   'Advanced_Tumor',
}
ann['tissue_type'] = ann['Sample_Origin'].map(origin_map).fillna(ann['Sample_Origin'])

ann['patient_id'] = ann['Sample'].str.extract(r'(\d+)$')[0].str.zfill(2).apply(lambda x: f"P{x}")

def assign_stage_chemo(sample):
    if any(sample.startswith(p) for p in ['LUNG_N', 'LUNG_T', 'LN_', 'NS_']):
        return 'Resectable', 'No'
    elif any(sample.startswith(p) for p in ['EBUS_', 'BRONCHO_', 'EFFUSION_']):
        return 'Advanced', 'Unknown'
    return 'Unknown', 'Unknown'

ann[['stage', 'chemotherapy']] = ann['Sample'].apply(
    lambda s: pd.Series(assign_stage_chemo(s))
)
ann['sample_id'] = ann['Sample']
ann['dataset']   = 'GSE131907'

print(ann[['tissue_type', 'stage', 'chemotherapy']].value_counts())

# Read matrix in chunks (genes x cells)
CHUNK_SIZE = 2000
chunks_data = []
gene_names  = []

print("Reading matrix...")
with gzip.open(mat_path, 'rt') as f:
    header     = f.readline().strip().split('\t')
    cell_names = header[1:]
    batch_rows = []

    for i, line in enumerate(f):
        parts = line.rstrip('\n').split('\t')
        gene_names.append(parts[0])
        batch_rows.append(np.array(parts[1:], dtype=np.float32))

        if len(batch_rows) == CHUNK_SIZE:
            chunks_data.append(sp.csr_matrix(np.vstack(batch_rows)))
            batch_rows = []
            print(f"  {i+1} genes processed...", end='\r')

    if batch_rows:
        chunks_data.append(sp.csr_matrix(np.vstack(batch_rows)))

X = sp.vstack(chunks_data).T.tocsr()
print(f"\nMatrix (cells x genes): {X.shape}")

# Align annotation (Index = Barcode_Sample)
obs_aligned = ann.reindex(cell_names)
missing = obs_aligned.isnull().any(axis=1).sum()
print(f"Annotation missing: {missing} cells")
if missing > 0:
    print("Sample of missing:", [c for c in cell_names[:20] if c not in ann.index])

adata = ad.AnnData(
    X   = X,
    obs = obs_aligned,
    var = pd.DataFrame(index=gene_names)
)
adata.obs_names = cell_names
adata.var_names = gene_names
adata.layers['counts'] = adata.X.copy()
print(adata)

# QC
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt']   = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb']   = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'],
                           percent_top=None, inplace=True)

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad    = np.median(np.abs(series - median))
    return median - nmads * mad * 1.4826, median + nmads * mad * 1.4826

n_genes_lo,  n_genes_hi  = mad_filter(adata.obs['n_genes_by_counts'])
n_counts_lo, n_counts_hi = mad_filter(adata.obs['total_counts'])
n_genes_lo  = max(n_genes_lo,  200)
n_counts_lo = max(n_counts_lo, 500)

print(f"n_genes  cutoff: [{n_genes_lo:.0f}, {n_genes_hi:.0f}]")
print(f"n_counts cutoff: [{n_counts_lo:.0f}, {n_counts_hi:.0f}]")

cells_before = adata.n_obs
adata = adata[
    (adata.obs['n_genes_by_counts']  >= n_genes_lo)  &
    (adata.obs['n_genes_by_counts']  <= n_genes_hi)  &
    (adata.obs['total_counts']       >= n_counts_lo) &
    (adata.obs['total_counts']       <= n_counts_hi) &
    (adata.obs['pct_counts_mt']      <  25)
].copy()
print(f"QC: {cells_before} -> {adata.n_obs} (removed {cells_before - adata.n_obs})")

# Scrublet per sample
adata.obs['doublet_score']     = 0.0
adata.obs['predicted_doublet'] = False

for sample in sorted(adata.obs['sample_id'].unique()):
    mask   = adata.obs['sample_id'] == sample
    n_cells = mask.sum()
    if n_cells < 50:
        print(f"  Skip {sample}: {n_cells} cells")
        continue
    try:
        scrub = scr.Scrublet(adata[mask].layers['counts'], expected_doublet_rate=0.06)
        scores, doublets = scrub.scrub_doublets(verbose=False)
        adata.obs.loc[mask, 'doublet_score']     = scores
        adata.obs.loc[mask, 'predicted_doublet'] = doublets
        print(f"  {sample}: {n_cells} cells, {doublets.sum()} doublets ({doublets.mean()*100:.1f}%)")
    except Exception as e:
        print(f"  {sample} error: {e}")

cells_before = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"Doublet removal: {cells_before} -> {adata.n_obs}")

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Field check
required = ['sample_id','patient_id','tissue_type','dataset',
            'stage','chemotherapy','n_genes_by_counts',
            'total_counts','pct_counts_mt','doublet_score']
for col in required:
    print(f"  {'OK' if col in adata.obs.columns else 'MISSING'}  {col}")

adata.write_h5ad(out_path, compression='gzip')
print(f"\nSaved: {out_path}")
print(adata)


import scanpy as sc

adata = sc.read_h5ad(r"${DATA_ROOT}/GSE131907/GSE131907_clean.h5ad")

missing_mask = adata.obs['tissue_type'].isnull()
print(f"Missing annotation cells: {missing_mask.sum()}")
print("\nFirst 10 cell names with missing annotation:")
print(adata.obs.loc[missing_mask].index[:10].tolist())

import pandas as pd, gzip
ann = pd.read_csv(r"${DATA_ROOT}/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz", sep='\t', index_col=0)
print("\nAnnotation index sample:", ann.index[:3].tolist())
print("Matrix cell_name sample (missing):", adata.obs.loc[missing_mask].index[:3].tolist())
print("Matrix cell_name sample (annotated):", adata.obs.loc[~missing_mask].index[:3].tolist())
