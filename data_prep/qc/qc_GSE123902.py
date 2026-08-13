#!/usr/bin/env python3
# Flattened code cells from qc_GSE123902.ipynb

import pandas as pd
# Peek at file head to confirm orientation
df = pd.read_csv(r"${DATA_ROOT}/GSE123902/extracted/GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz",
                 nrows=3, index_col=0)
print(df.shape)
print(df.index[:5].tolist())
print(df.columns[:5].tolist())


df2 = pd.read_csv(r"${DATA_ROOT}/GSE123902/extracted/GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz",
                  nrows=5, index_col=0)
print(df2.iloc[:3, :3])
print(df2.dtypes[:3])
print("has_floats:", (df2.values % 1 != 0).any())


import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scrublet as scr
from scipy.sparse import csr_matrix
import os

file_meta = {
    "GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv.gz":  ("LX653", "LX653", "Tumor"),
    "GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv.gz":  ("LX661", "LX661", "Tumor"),
    "GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz":      ("LX666", "LX666", "Metastasis"),
    "GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv.gz":  ("LX675_T", "LX675", "Tumor"),
    "GSM3516666_MSK_LX675_NORMAL_dense.csv.gz":          ("LX675_N", "LX675", "Normal"),
    "GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv.gz":  ("LX676", "LX676", "Tumor"),
    "GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz":     ("LX255B", "LX255B", "Metastasis"),
    "GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv.gz":  ("LX679", "LX679", "Tumor"),
    "GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv.gz":  ("LX680", "LX680", "Tumor"),
    "GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz":      ("LX681", "LX681", "Metastasis"),
    "GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv.gz":  ("LX682_T", "LX682", "Tumor"),
    "GSM3516673_MSK_LX682_NORMAL_dense.csv.gz":          ("LX682_N", "LX682", "Normal"),
    "GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv.gz":  ("LX684_T", "LX684", "Tumor"),
    "GSM3516675_MSK_LX684_NORMAL_dense.csv.gz":          ("LX684_N", "LX684", "Normal"),
    "GSM3516676_MSK_LX685_NORMAL_dense.csv.gz":          ("LX685", "LX685", "Normal"),
    "GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz":      ("LX699", "LX699", "Metastasis"),
    "GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz":      ("LX701", "LX701", "Metastasis"),
}

data_dir = r"${DATA_ROOT}/GSE123902/extracted"
out_dir  = r"${DATA_ROOT}/GSE123902/qc_out"
os.makedirs(out_dir, exist_ok=True)

def mad_bounds(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return median - nmads * mad * 1.4826, median + nmads * mad * 1.4826

qc_summary = []
adatas = []

for fname, (sample_id, patient_id, tissue_type) in file_meta.items():
    print(f"\n{'='*50}\nProcessing: {sample_id}")

    df = pd.read_csv(os.path.join(data_dir, fname), index_col=0)
    df.columns = df.columns.str.upper()

    adata = ad.AnnData(X=csr_matrix(df.values, dtype=np.float32),
                       obs=pd.DataFrame(index=df.index.astype(str)),
                       var=pd.DataFrame(index=df.columns))
    adata.obs['sample_id']   = sample_id
    adata.obs['patient_id']  = patient_id
    adata.obs['tissue_type'] = tissue_type
    adata.obs['dataset']     = 'GSE123902'

    n_before = adata.n_obs

    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt']   = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb']   = adata.var_names.isin(['HBA1','HBA2','HBB','HBM'])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'],
                               percent_top=None, inplace=True)

    _, n_genes_upper = mad_bounds(adata.obs['n_genes_by_counts'])
    _, n_counts_upper = mad_bounds(adata.obs['total_counts'])

    mask = (
        (adata.obs['n_genes_by_counts'] >= 200) &
        (adata.obs['n_genes_by_counts'] <= n_genes_upper) &
        (adata.obs['total_counts'] >= 500) &
        (adata.obs['total_counts'] <= n_counts_upper) &
        (adata.obs['pct_counts_mt'] < 25) &
        (adata.obs['pct_counts_hb'] < 5)
    )
    adata = adata[mask].copy()
    n_after_qc = adata.n_obs

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
    scores, predicted = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                              n_prin_comps=30)
    adata.obs['doublet_score']     = scores
    adata.obs['predicted_doublet'] = predicted
    adata = adata[~predicted].copy()
    n_after_doublet = adata.n_obs

    adata.layers['counts'] = adata.X.copy()

    print(f"  {n_before} -> after_QC {n_after_qc} -> after_doublet {n_after_doublet}")
    print(f"  n_genes cutoff: 200..{n_genes_upper:.0f}, n_counts cutoff: 500..{n_counts_upper:.0f}")

    qc_summary.append({
        'sample_id': sample_id,
        'tissue_type': tissue_type,
        'n_before': n_before,
        'n_after_qc': n_after_qc,
        'n_after_doublet': n_after_doublet,
        'pct_removed': round((1 - n_after_doublet/n_before)*100, 1)
    })

    adatas.append(adata)

df_summary = pd.DataFrame(qc_summary)
print("\n", df_summary.to_string())
df_summary.to_csv(os.path.join(out_dir, 'qc_summary.csv'), index=False)

print("\nMerging...")
adata_merged = ad.concat(adatas, join='inner', merge='same')
adata_merged.var_names_make_unique()
print(f"Merged: {adata_merged.n_obs} cells x {adata_merged.n_vars} genes")

out_path = os.path.join(out_dir, 'GSE123902_qc.h5ad')
adata_merged.write_h5ad(out_path)
print(f"Saved: {out_path}")


import scanpy as sc
adata = sc.read_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
adata.obs_names_make_unique()
adata.write_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
print("obs_names made unique")


import scanpy as sc
adata = sc.read_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
print(adata)
print("\n--- obs head ---")
print(adata.obs.head())
print("\n--- var head ---")
print(adata.var.head())
print("\n--- layers ---")
print(list(adata.layers.keys()))


import scanpy as sc
import numpy as np

adata = sc.read_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")

# Confirm .X is integer raw counts
print("X is integer:", np.allclose(adata.X.data % 1, 0))

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

adata.write_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
print("done, layers:", list(adata.layers.keys()))


import urllib.request

gsm_primary = {
    'GSM3516662': 'LX653',
    'GSM3516663': 'LX661',
    'GSM3516665': 'LX675_T',
    'GSM3516666': 'LX675_N',
    'GSM3516667': 'LX676',
    'GSM3516669': 'LX679',
    'GSM3516670': 'LX680',
    'GSM3516672': 'LX682_T',
    'GSM3516673': 'LX682_N',
    'GSM3516674': 'LX684_T',
    'GSM3516675': 'LX684_N',
    'GSM3516676': 'LX685',
}

for gsm, sample in gsm_primary.items():
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full"
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            text = r.read().decode()
        print(f"\n=== {gsm} ({sample}) ===")
        for line in text.split('\n'):
            if any(k in line.lower() for k in ['stage', 'chemo', 'characteristics', 'source']):
                print(line.strip())
    except Exception as e:
        print(f"{gsm} failed: {e}")


import scanpy as sc

adata = sc.read_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")

meta = {
    'LX653':   ('Tumor',       '1A',  'No'),
    'LX661':   ('Tumor',       '1A',  'No'),
    'LX666':   ('Bone_met',    'IV',  'Yes'),
    'LX675_T': ('Tumor',       'IV',  'No'),
    'LX675_N': ('Normal',      'NA',  'No'),
    'LX676':   ('Tumor',       '1A',  'No'),
    'LX255B':  ('Brain_met',   'IV',  'Yes'),
    'LX679':   ('Tumor',       'IIA', 'Yes'),
    'LX680':   ('Tumor',       'IB',  'No'),
    'LX681':   ('Brain_met',   'IV',  'Yes'),
    'LX682_T': ('Tumor',       'IB',  'No'),
    'LX682_N': ('Normal',      'NA',  'No'),
    'LX684_T': ('Tumor',       'IA',  'No'),
    'LX684_N': ('Normal',      'NA',  'No'),
    'LX685':   ('Normal',      'NA',  'No'),
    'LX699':   ('Adrenal_met', 'IV',  'Yes'),
    'LX701':   ('Brain_met',   'IV',  'Yes'),
}

new_tissue = ['Bone_met', 'Brain_met', 'Adrenal_met']
adata.obs['tissue_type'] = adata.obs['tissue_type'].cat.add_categories(
    [c for c in new_tissue if c not in adata.obs['tissue_type'].cat.categories])

adata.obs['stage'] = 'Unknown'
adata.obs['chemotherapy'] = 'Unknown'

for sample, (tissue, stage, chemo) in meta.items():
    mask = adata.obs['sample_id'] == sample
    adata.obs.loc[mask, 'tissue_type']  = tissue
    adata.obs.loc[mask, 'stage']        = stage
    adata.obs.loc[mask, 'chemotherapy'] = chemo

adata.obs['tissue_type'] = adata.obs['tissue_type'].cat.remove_unused_categories()

print(adata.obs[['sample_id','tissue_type','stage','chemotherapy']].drop_duplicates())
adata.write_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
print("saved")


import scanpy as sc
import numpy as np

adata = sc.read_h5ad(r"${DATA_ROOT}/GSE123902/qc_out/GSE123902_qc.h5ad")
print("layers:", list(adata.layers.keys()))
print("X is integer:", np.allclose(adata.X.data % 1, 0))
