#!/usr/bin/env python3
# Flattened code cells from qc_GSE148071.ipynb

import tarfile
import gzip
import io
import pandas as pd
import numpy as np
import os


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


tar_path = r"${DATA_ROOT}/GSE148071/GSE148071_RAW.tar"
output_dir = r"${DATA_ROOT}/GSE148071"
os.makedirs(output_dir, exist_ok=True)

# Marker genes (paper Fig. S2)
luad_markers = ['NAPSA', 'NKX2-1']
lusc_markers = ['TP63', 'KRT5']
cancer_marker = 'EPCAM'
normal_epithelial = ['SFTPC', 'SCGB1A1', 'CLDN18', 'KRT5']

print("Processing 42 samples in GSE148071_RAW.tar...\n")

luad_samples = []
lusc_samples = []
other_samples = []

with tarfile.open(tar_path, 'r') as tar:
    for member in tar.getmembers():
        if not member.name.endswith('_exp.txt.gz'):
            continue

        gsm = member.name.split('_')[0]   # GSM4453576
        p_id = member.name.split('_')[1].replace('.txt.gz', '')  # P1

        f = tar.extractfile(member)
        with gzip.GzipFile(fileobj=f) as gz:
            content = gz.read()

        # Print head of first file only, to confirm format
        if len(luad_samples) + len(lusc_samples) + len(other_samples) == 0:
            print("File head preview:")
            print(content.decode('utf-8')[:1000])
            print("-" * 60)

        df = pd.read_csv(io.BytesIO(content), sep='\t', index_col=0)

        # Coerce numeric, drop empty columns
        df = df.apply(pd.to_numeric, errors='coerce')
        df = df.dropna(how='all', axis=1)

        # Cancer cell selection (genes are row index)
        if cancer_marker not in df.index:
            print(f"{gsm} ({p_id}) no EPCAM, skip")
            continue

        # Transpose to cells x genes
        df = df.T

        cancer_mask = (df[cancer_marker] > 0) & (df[normal_epithelial].mean(axis=1) < 0.5)
        cancer_cells = df[cancer_mask]

        if len(cancer_cells) < 10:
            print(f"{gsm} ({p_id}) too few cancer cells ({len(cancer_cells)}), skip")
            continue

        # Marker percent (paper method)
        luad_pct = (cancer_cells[luad_markers] > 0).mean().mean() * 100
        lusc_pct = (cancer_cells[lusc_markers] > 0).mean().mean() * 100

        result = f"{gsm} ({p_id}) | LUAD: {luad_pct:.1f}% | LUSC: {lusc_pct:.1f}%"

        if luad_pct > 5 and lusc_pct < 5:
            luad_samples.append((gsm, p_id))
            print(f"OK {result} -> LUAD")
        elif lusc_pct > 5 and luad_pct < 5:
            lusc_samples.append((gsm, p_id))
            print(f"OK {result} -> LUSC")
        else:
            other_samples.append((gsm, p_id))
            print(f"   {result} -> other/mixed")

print("\n" + "="*70)
print(f"LUAD samples ({len(luad_samples)}, matches paper):")
for gsm, p in sorted(luad_samples, key=lambda x: int(x[1][1:])):
    print(f"  {gsm} ({p})")

print(f"\nLUSC samples ({len(lusc_samples)}):")
for gsm, p in sorted(lusc_samples, key=lambda x: int(x[1][1:])):
    print(f"  {gsm} ({p})")

pd.DataFrame(luad_samples, columns=['GSM', 'P_ID']).to_csv(
    os.path.join(output_dir, 'LUAD_samples_list.csv'), index=False)

print(f"\nSaved: {output_dir}\\LUAD_samples_list.csv")


import tarfile
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import warnings
warnings.filterwarnings('ignore')

tar_path = r"${DATA_ROOT}/GSE148071/GSE148071_RAW.tar"
extract_dir = r"${DATA_ROOT}/GSE148071/raw_extracted"

os.makedirs(extract_dir, exist_ok=True)
with tarfile.open(tar_path, 'r') as tar:
    _safe_extract(tar, extract_dir)
    members = tar.getnames()
print(f"Extracted {len(members)} files")
for m in members[:20]:
    print(m)

for root, dirs, files in os.walk(extract_dir):
    for f in sorted(files)[:10]:
        fpath = os.path.join(root, f)
        print(f"{f}  ({os.path.getsize(fpath)} bytes)")
        with open(fpath, 'r', errors='ignore') as fh:
            for i, line in enumerate(fh):
                if i < 3:
                    print(f"  {line.rstrip()[:120]}")
                else:
                    break
    break


import gzip
import os

extract_dir = r"${DATA_ROOT}/GSE148071/raw_extracted"
f = os.path.join(extract_dir, "GSM4453576_P1_exp.txt.gz")

with gzip.open(f, 'rt') as fh:
    for i, line in enumerate(fh):
        if i < 5:
            cols = line.rstrip().split('\t')
            print(f"row {i}: {len(cols)} cols, first 5: {cols[:5]}")
        else:
            break

with gzip.open(f, 'rt') as fh:
    total = sum(1 for _ in fh)
print(f"total rows: {total}")


import gzip
import os
import pandas as pd
import numpy as np

extract_dir = r"${DATA_ROOT}/GSE148071/raw_extracted"

# LUAD markers: NAPSA, NKX2-1 (TTF-1)
# LUSC markers: KRT5, TP63, DSG3
# Extra: EPCAM (epithelial), KRT7 (adeno), KRT14 (squamous basal)
markers = ['NAPSA', 'NKX2-1', 'KRT5', 'TP63', 'DSG3', 'EPCAM', 'KRT7', 'KRT14', 'SFTPC', 'SFTPA1']

results = []

for fname in sorted(os.listdir(extract_dir)):
    if not fname.endswith('_exp.txt.gz'):
        continue
    sample = fname.split('_')[1]
    fpath = os.path.join(extract_dir, fname)

    marker_data = {}
    n_cells = 0

    with gzip.open(fpath, 'rt') as fh:
        header = fh.readline().rstrip().split('\t')
        n_cells = len(header)

        for line in fh:
            parts = line.rstrip().split('\t')
            gene = parts[0]
            if gene in markers:
                vals = np.array(parts[1:], dtype=float)
                marker_data[gene] = vals

    row = {'Sample': sample, 'n_cells': n_cells}
    for g in markers:
        if g in marker_data:
            arr = marker_data[g]
            row[f'{g}_pct'] = np.mean(arr > 0) * 100
            row[f'{g}_mean'] = np.mean(arr)
        else:
            row[f'{g}_pct'] = 0.0
            row[f'{g}_mean'] = 0.0

    luad_score = np.mean([row.get('NAPSA_pct', 0), row.get('NKX2-1_pct', 0)]) / 100
    lusc_score = np.mean([row.get('KRT5_pct', 0), row.get('TP63_pct', 0), row.get('DSG3_pct', 0)]) / 100

    if luad_score < 0.05 and lusc_score < 0.05:
        subtype = 'NSCLC_NOS'
    elif luad_score >= lusc_score:
        subtype = 'LUAD'
    else:
        subtype = 'LUSC'

    row['LUAD_score'] = round(luad_score, 4)
    row['LUSC_score'] = round(lusc_score, 4)
    row['Subtype'] = subtype

    results.append(row)
    print(f"{sample}: {n_cells} cells | LUAD={luad_score:.3f} LUSC={lusc_score:.3f} -> {subtype}")

df = pd.DataFrame(results)

print("\n========== Subtype counts ==========")
print(df['Subtype'].value_counts())
print(f"\nLUAD samples: {df[df.Subtype=='LUAD']['Sample'].tolist()}")
print(f"LUSC samples: {df[df.Subtype=='LUSC']['Sample'].tolist()}")
print(f"NSCLC_NOS samples: {df[df.Subtype=='NSCLC_NOS']['Sample'].tolist()}")

out_path = r"${DATA_ROOT}/GSE148071/sample_subtype_classification.csv"
df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

print("\n========== Key marker expression % ==========")
cols_show = ['Sample', 'n_cells', 'NAPSA_pct', 'NKX2-1_pct', 'KRT5_pct', 'TP63_pct', 'DSG3_pct', 'Subtype']
existing_cols = [c for c in cols_show if c in df.columns]
print(df[existing_cols].to_string(index=False))


import gzip
import os
import numpy as np
import pandas as pd

extract_dir = r"${DATA_ROOT}/GSE148071/raw_extracted"
markers = ['NAPSA', 'NKX2-1', 'KRT5', 'TP63', 'DSG3', 'EPCAM']

results = []

for fname in sorted(os.listdir(extract_dir)):
    if not fname.endswith('_exp.txt.gz'):
        continue
    sample = fname.split('_')[1]
    fpath = os.path.join(extract_dir, fname)

    marker_data = {}
    n_cells = 0

    with gzip.open(fpath, 'rt') as fh:
        header = fh.readline().rstrip().split('\t')
        n_cells = len(header)
        for line in fh:
            parts = line.rstrip().split('\t')
            gene = parts[0]
            if gene in markers:
                marker_data[gene] = np.array(parts[1:], dtype=float)

    # Filter epithelial / cancer cells with EPCAM>0
    if 'EPCAM' in marker_data:
        epcam_mask = marker_data['EPCAM'] > 0
    else:
        epcam_mask = np.ones(n_cells, dtype=bool)  # fallback: all cells

    n_epcam = epcam_mask.sum()

    row = {'Sample': sample, 'n_cells': n_cells, 'n_EPCAM+': int(n_epcam)}

    # Marker percent within EPCAM+ cells
    for g in ['NAPSA', 'NKX2-1', 'KRT5', 'TP63', 'DSG3']:
        if g in marker_data and n_epcam > 0:
            row[f'{g}_pct_epi'] = np.mean(marker_data[g][epcam_mask] > 0) * 100
        else:
            row[f'{g}_pct_epi'] = 0.0

    # Score on EPCAM+ cells
    if n_epcam > 0:
        luad_score = np.mean([row['NAPSA_pct_epi'], row['NKX2-1_pct_epi']]) / 100
        lusc_score = np.mean([row['KRT5_pct_epi'], row['TP63_pct_epi'], row['DSG3_pct_epi']]) / 100
    else:
        luad_score = 0
        lusc_score = 0

    # Paper rule: both <5% -> NSCLC_NOS, otherwise higher score wins
    if luad_score < 0.05 and lusc_score < 0.05:
        subtype = 'NSCLC_NOS'
    elif luad_score >= lusc_score:
        subtype = 'LUAD'
    else:
        subtype = 'LUSC'

    row['LUAD_score'] = round(luad_score, 4)
    row['LUSC_score'] = round(lusc_score, 4)
    row['Subtype'] = subtype
    results.append(row)
    print(f"{sample}: {n_cells} cells, {n_epcam} EPCAM+ | LUAD={luad_score:.3f} LUSC={lusc_score:.3f} -> {subtype}")

df = pd.DataFrame(results)
print("\n========== Subtype counts ==========")
print(df['Subtype'].value_counts())
print(f"\nLUAD: {sorted(df[df.Subtype=='LUAD']['Sample'].tolist(), key=lambda x: int(x[1:]))}")
print(f"LUSC: {sorted(df[df.Subtype=='LUSC']['Sample'].tolist(), key=lambda x: int(x[1:]))}")
print(f"NSCLC_NOS: {sorted(df[df.Subtype=='NSCLC_NOS']['Sample'].tolist(), key=lambda x: int(x[1:]))}")

out_path = r"${DATA_ROOT}/GSE148071/sample_subtype_classification_v2.csv"
df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")


import scanpy as sc
import scrublet as scr
import anndata as ad
import pandas as pd
import numpy as np
import gzip
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore')

extract_dir = r"${DATA_ROOT}/GSE148071/raw_extracted"
output_dir = r"${DATA_ROOT}/GSE148071/processed"
fig_dir = r"${DATA_ROOT}/GSE148071/processed/figures"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

# 18 LUAD samples
luad_samples = ['P2','P5','P8','P9','P11','P12','P13','P16','P20','P21',
                'P24','P28','P29','P32','P33','P35','P38','P39']

gsm_map = {f'P{i}': f'GSM{4453575+i}' for i in range(1, 43)}

print("=" * 60)
print("Step 1: read LUAD samples, build AnnData")
print("=" * 60)

adatas = []
for sample in luad_samples:
    gsm = gsm_map[sample]
    fname = f"{gsm}_{sample}_exp.txt.gz"
    fpath = os.path.join(extract_dir, fname)

    if not os.path.exists(fpath):
        print(f"  [WARN] {fname} not found, skip")
        continue

    # rows = genes, cols = cells
    df = pd.read_csv(fpath, sep='\t', index_col=0, compression='gzip')

    # Transpose to cells x genes, sparsify
    adata = ad.AnnData(X=csr_matrix(df.values.T.astype(np.float32)),
                       obs=pd.DataFrame(index=df.columns),
                       var=pd.DataFrame(index=df.index))

    adata.obs['sample_id'] = sample
    adata.obs['patient_id'] = sample
    adata.obs['dataset'] = 'GSE148071'
    adata.obs['tissue_type'] = 'Tumor'  # all biopsy tumor tissue
    adata.obs['cancer_type'] = 'LUAD'
    adata.obs['stage'] = 'III/IV'  # paper: advanced NSCLC
    adata.obs['platform'] = 'GEXSCOPE'  # Singleron platform

    # Prefix barcodes with sample to avoid collisions
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]

    adatas.append(adata)
    print(f"  {sample} ({gsm}): {adata.n_obs} cells x {adata.n_vars} genes")

print(f"\nLoaded {len(adatas)} samples")

print("\n" + "=" * 60)
print("Step 2: per-sample QC")
print("=" * 60)

def mad_filter(series, nmads=5):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    lower = max(median - nmads * mad * 1.4826, 0)
    upper = median + nmads * mad * 1.4826
    return lower, upper

qc_summary = []

for i, adata in enumerate(adatas):
    sample = adata.obs['sample_id'].iloc[0]
    n_before = adata.n_obs

    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'],
                               percent_top=None, inplace=True)

    ng_lo, ng_hi = mad_filter(adata.obs['n_genes_by_counts'], nmads=5)
    nc_lo, nc_hi = mad_filter(adata.obs['total_counts'], nmads=5)

    # Hard lower bounds
    ng_lo = max(ng_lo, 200)
    nc_lo = max(nc_lo, 500)

    keep = (
        (adata.obs['n_genes_by_counts'] >= ng_lo) &
        (adata.obs['n_genes_by_counts'] <= ng_hi) &
        (adata.obs['total_counts'] >= nc_lo) &
        (adata.obs['total_counts'] <= nc_hi) &
        (adata.obs['pct_counts_mt'] < 25) &
        (adata.obs['pct_counts_hb'] < 5)
    )

    n_after = keep.sum()
    adata_filtered = adata[keep].copy()
    adatas[i] = adata_filtered

    qc_summary.append({
        'sample': sample,
        'cells_before': n_before,
        'genes_after_filter': adata_filtered.n_vars,
        'ng_lo': round(ng_lo), 'ng_hi': round(ng_hi),
        'nc_lo': round(nc_lo), 'nc_hi': round(nc_hi),
        'cells_after_qc': n_after,
        'pct_removed': round((1 - n_after / n_before) * 100, 1)
    })

    print(f"  {sample}: {n_before} -> {n_after} cells "
          f"(-{round((1-n_after/n_before)*100,1)}%) | "
          f"genes={adata_filtered.n_vars} | "
          f"nGene=[{round(ng_lo)},{round(ng_hi)}] mt<25%")

df_qc = pd.DataFrame(qc_summary)
df_qc.to_csv(os.path.join(output_dir, 'qc_summary_GSE148071.csv'), index=False)
print(f"\nQC summary saved")

# QC violins per sample
print("\nGenerating QC plots...")
fig, axes = plt.subplots(3, 1, figsize=(max(len(adatas)*0.8, 12), 12))

adata_plot = ad.concat(adatas)
for ax, metric in zip(axes, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
    sc.pl.violin(adata_plot, metric, groupby='sample_id', ax=ax, show=False, rotation=45)
    ax.set_title(metric)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'qc_violin_post_filter.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  QC violin plot saved")

print("\n" + "=" * 60)
print("Step 3: per-sample Scrublet")
print("=" * 60)

doublet_summary = []

for i, adata in enumerate(adatas):
    sample = adata.obs['sample_id'].iloc[0]
    n_before = adata.n_obs

    if n_before < 50:
        print(f"  {sample}: too few cells ({n_before}), skip doublet detection")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False
        continue

    try:
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85,
            n_prin_comps=min(30, n_before - 1)
        )
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets

        n_doublet = predicted_doublets.sum()
        adatas[i] = adata[~adata.obs['predicted_doublet']].copy()
        n_after = adatas[i].n_obs

        doublet_summary.append({
            'sample': sample, 'cells_before': n_before,
            'doublets_detected': n_doublet,
            'cells_after': n_after,
            'doublet_rate': round(n_doublet / n_before * 100, 1)
        })
        print(f"  {sample}: {n_before} -> {n_after} cells | "
              f"doublets={n_doublet} ({round(n_doublet/n_before*100,1)}%)")
    except Exception as e:
        print(f"  {sample}: Scrublet failed ({e}), skip")
        adata.obs['doublet_score'] = 0.0
        adata.obs['predicted_doublet'] = False

df_doublet = pd.DataFrame(doublet_summary)
df_doublet.to_csv(os.path.join(output_dir, 'doublet_summary_GSE148071.csv'), index=False)

print("\n" + "=" * 60)
print("Merge and save")
print("=" * 60)

adata_merged = ad.concat(adatas)
adata_merged.var_names_make_unique()

# Recompute var flags after merge
adata_merged.var['mt'] = adata_merged.var_names.str.startswith('MT-')
adata_merged.var['ribo'] = adata_merged.var_names.str.startswith(('RPS', 'RPL'))
adata_merged.var['hb'] = adata_merged.var_names.isin(['HBA1', 'HBA2', 'HBB', 'HBM'])

adata_merged.layers['counts'] = adata_merged.X.copy()

total_before = df_qc['cells_before'].sum()
total_after = adata_merged.n_obs
print(f"Total cells before: {total_before}")
print(f"After QC + doublet: {total_after} ({round(total_after/total_before*100,1)}% retained)")
print(f"Genes: {adata_merged.n_vars}")
print(f"Samples: {adata_merged.obs['sample_id'].nunique()}")

out_path = os.path.join(output_dir, 'GSE148071_LUAD_clean.h5ad')
adata_merged.write(out_path)
print(f"\nSaved: {out_path}")
print(f"File size: {os.path.getsize(out_path)/1024/1024:.1f} MB")

print("\n========== Final cells per sample ==========")
print(adata_merged.obs['sample_id'].value_counts().sort_index().to_string())


import scanpy as sc
import anndata as ad
import os

# Re-merge with outer join (fill missing genes with 0)
adata_merged = ad.concat(adatas, join='outer')
adata_merged.var_names_make_unique()

# Replace NaN in sparse matrix with 0
import numpy as np
from scipy.sparse import issparse
if issparse(adata_merged.X):
    adata_merged.X.data = np.nan_to_num(adata_merged.X.data)
else:
    adata_merged.X = np.nan_to_num(adata_merged.X)

# Re-filter genes after merge: at least 0.1% of total cells
min_cells = max(int(adata_merged.n_obs * 0.001), 3)
sc.pp.filter_genes(adata_merged, min_cells=min_cells)

adata_merged.var['mt'] = adata_merged.var_names.str.startswith('MT-')
adata_merged.var['ribo'] = adata_merged.var_names.str.startswith(('RPS', 'RPL'))
adata_merged.var['hb'] = adata_merged.var_names.isin(['HBA1','HBA2','HBB','HBM'])

adata_merged.layers['counts'] = adata_merged.X.copy()

print(f"Merged: {adata_merged.n_obs} cells x {adata_merged.n_vars} genes")
print(f"MT genes: {adata_merged.var['mt'].sum()}")

output_dir = r"${DATA_ROOT}/GSE148071/processed"
out_path = os.path.join(output_dir, 'GSE148071_LUAD_clean.h5ad')
adata_merged.write(out_path)
print(f"Saved: {out_path} ({os.path.getsize(out_path)/1024/1024:.1f} MB)")