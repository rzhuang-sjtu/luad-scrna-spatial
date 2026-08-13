#!/usr/bin/env python3
# Flattened code cells from qc_GSE253013.ipynb

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
import scrublet as scr
import anndata as ad
import warnings
warnings.filterwarnings('ignore')

base = "${DATA_ROOT}/GSE253013/"

print("Reading matrix...")
X = mmread(f"{base}matrix.mtx").T.tocsr()  # cells x genes
meta = pd.read_csv(f"{base}metadata.csv", index_col=0)
genes = pd.read_csv(f"{base}genes.csv")

adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes['gene'].values))
adata.var_names_make_unique()
print(f"Raw: {adata.shape[0]} cells, {adata.shape[1]} genes")

adata.obs['sample_id'] = adata.obs['library_id'].astype(str)
adata.obs['patient_id'] = adata.obs['SampleGroup'].astype(str)
adata.obs['dataset'] = 'GSE253013'

# Unify tissue_type naming
tissue_map = {'T': 'Tumor', 'NAT': 'Normal'}
adata.obs['tissue_type_clean'] = adata.obs['tissue_type'].map(tissue_map)

# stage
adata.obs['stage'] = adata.obs['AJCC_Stage_8th_Edition'].astype(str)
adata.obs['chemotherapy'] = 'treatment_naive'

# Flag MRC008 ANT as CD45+ sorted
adata.obs['cd45_sorted'] = False
mask_cd45 = (adata.obs['patient_id'] == 'MRC008') & (adata.obs['tissue_type'] == 'NAT')
adata.obs.loc[mask_cd45, 'cd45_sorted'] = True
print(f"CD45+ sorted cells (MRC008 ANT): {mask_cd45.sum()}")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.var['hb'] = adata.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                           percent_top=None, inplace=True)

print("\nQC summary before filtering:")
for col in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_hb']:
    print(f"  {col}: median={adata.obs[col].median():.1f}, "
          f"mean={adata.obs[col].mean():.1f}, "
          f"max={adata.obs[col].max():.1f}")

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return median - nmads * mad * 1.4826, median + nmads * mad * 1.4826

qc_results = []
keep_mask = np.ones(adata.n_obs, dtype=bool)

for patient in adata.obs['patient_id'].unique():
    idx = adata.obs['patient_id'] == patient
    sub = adata.obs.loc[idx]
    n_before = idx.sum()

    # MAD thresholds
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

    keep_mask[idx.values] &= mask.values
    n_after = mask.sum()

    qc_results.append({
        'patient': patient,
        'before': n_before,
        'after': n_after,
        'removed_pct': f"{(1 - n_after/n_before)*100:.1f}%",
        'ng_range': f"{ng_lo:.0f}-{ng_hi:.0f}",
        'nc_range': f"{nc_lo:.0f}-{nc_hi:.0f}"
    })

qc_df = pd.DataFrame(qc_results)
print("\nQC per patient:")
print(qc_df.to_string(index=False))

adata = adata[keep_mask].copy()
print(f"\nAfter QC: {adata.shape[0]} cells")

sc.pp.filter_genes(adata, min_cells=3)
print(f"After gene filter: {adata.shape[1]} genes")

print("\nRunning Scrublet per sample...")
adata.obs['doublet_score'] = 0.0
adata.obs['predicted_doublet'] = False

sample_list = adata.obs['sample_id'].unique()
for i, sample in enumerate(sample_list):
    idx = adata.obs['sample_id'] == sample
    n_cells = idx.sum()
    if n_cells < 100:
        print(f"  {sample}: skipped ({n_cells} cells)")
        continue

    sub_X = adata[idx].X.copy()
    try:
        scrub = scr.Scrublet(sub_X, expected_doublet_rate=0.06)
        scores, preds = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                              min_gene_variability_pctl=85,
                                              verbose=False)
        adata.obs.loc[idx, 'doublet_score'] = scores
        adata.obs.loc[idx, 'predicted_doublet'] = preds
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  [{i+1}/{len(sample_list)}] {sample}: {preds.sum()}/{n_cells} doublets "
                  f"({preds.sum()/n_cells*100:.1f}%)")
    except Exception as e:
        print(f"  {sample}: Scrublet failed - {e}")

n_doublets = adata.obs['predicted_doublet'].sum()
print(f"\nTotal doublets: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")

adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"After doublet removal: {adata.shape[0]} cells")

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Tidy obs columns
keep_cols = ['sample_id', 'patient_id', 'dataset', 'tissue_type_clean',
             'stage', 'chemotherapy', 'cd45_sorted',
             'PatientID', 'Age', 'Sex', 'Histology', 'Site',
             'AJCC_TNM_8th_Edition', 'SmokingStatus', 'cell_type',
             'n_genes_by_counts', 'total_counts', 'pct_counts_mt',
             'pct_counts_ribo', 'pct_counts_hb', 'doublet_score']
adata.obs = adata.obs[[c for c in keep_cols if c in adata.obs.columns]]
adata.obs.rename(columns={'tissue_type_clean': 'tissue_type'}, inplace=True)

out_path = f"{base}GSE253013_cleaned.h5ad"
adata.write(out_path)
print(f"\nSaved: {out_path}")
print(f"Final: {adata.shape[0]} cells, {adata.shape[1]} genes")
print(f"Layers: {list(adata.layers.keys())}")
print(f"\nTissue distribution:")
print(adata.obs['tissue_type'].value_counts())
print(f"\nPatient distribution:")
print(adata.obs['patient_id'].value_counts())
