"""Step 10e: Compute DEG + GO enrichment for the two refined-only subtypes
Macro_FOLR2 and Macro_SPP1, then append to fig4/myeloid_go_enrichment.csv.

Inputs:
  ~/luad/data/processed/luad_myeloid.h5ad
  ~/luad/results/step10c_obs_labels.csv.gz   (refined labels from 10c)
  ${WORK_ROOT}/luad_figures/fig4/myeloid_go_enrichment.csv  (existing 11-sub)

Outputs (replaced/appended):
  ~/luad/results/step10e_folr2_spp1_markers.csv
  ~/luad/results/step10_myeloid_go_enrichment.csv     (full 13-sub top10)
  ${WORK_ROOT}/luad_figures/fig4/myeloid_go_enrichment.csv
"""
from __future__ import annotations
import os, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp

H5  = Path.home() / 'luad/data/processed/luad_myeloid.h5ad'
LBL = Path.home() / 'luad/results/step10c_obs_labels.csv.gz'
RES = Path.home() / 'luad/results'
FIG = Path('${WORK_ROOT}/luad_figures/fig4')
GENE_SETS = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
TARGET = ['Macro_FOLR2', 'Macro_SPP1']

def log(m): print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)

def main():
    t0 = time.time()
    log(f'load {H5}')
    a = sc.read_h5ad(H5)
    log(f'  shape={a.shape}')

    log(f'load refined labels {LBL}')
    lbl = pd.read_csv(LBL, index_col=0)
    a.obs['myeloid_subtype_refined'] = (
        lbl['myeloid_subtype_refined'].reindex(a.obs.index).astype('category')
    )
    counts = a.obs['myeloid_subtype_refined'].value_counts()
    log(f'  refined counts:\n{counts.to_string()}')

    for t in TARGET:
        if t not in counts.index:
            log(f'WARN {t} not present; skip'); 
    targets = [t for t in TARGET if t in counts.index]
    if not targets:
        log('no targets present, abort'); return

    log(f'rank_genes_groups for {targets}')
    sc.tl.rank_genes_groups(a, groupby='myeloid_subtype_refined',
                              groups=targets, reference='rest',
                              method='wilcoxon', n_genes=2000, pts=True)

    # extract per subtype
    deg_rows = []
    for t in targets:
        df = sc.get.rank_genes_groups_df(a, group=t)
        df['subtype'] = t
        df = df.rename(columns={'names':'gene','logfoldchanges':'logFC',
                                 'pvals_adj':'pval_adj','scores':'score'})
        deg_rows.append(df[['subtype','gene','score','logFC','pval_adj']])
    deg = pd.concat(deg_rows, ignore_index=True)
    deg.to_csv(RES/'step10e_folr2_spp1_markers.csv', index=False)
    log(f'  DEG rows: {len(deg)}; saved')

    # enrichr per subtype
    enr_rows = []
    for t, d in deg.groupby('subtype'):
        tight = d[(d['logFC'] > 0.5) & (d['pval_adj'] < 0.05)]
        if len(tight) >= 10:
            genes = tight['gene'].dropna().astype(str).tolist()
            src = 'logFC>0.5 & padj<0.05'
        else:
            genes = d.sort_values('score', ascending=False).head(100)['gene'].astype(str).tolist()
            src = 'top 100 by score'
        log(f'  {t}: {len(genes)} genes ({src})')
        if len(genes) < 5:
            continue
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=GENE_SETS,
                             organism='human', outdir=None, no_plot=True)
            for _, row in enr.results.iterrows():
                enr_rows.append({
                    'subtype': t,
                    'gene_set': row['Gene_set'],
                    'term': row['Term'],
                    'overlap': row['Overlap'],
                    'p_value': row['P-value'],
                    'adj_p_value': row['Adjusted P-value'],
                    'odds_ratio': row.get('Odds Ratio', None),
                    'combined_score': row.get('Combined Score', None),
                    'genes': row['Genes'],
                    'source_filter': src,
                })
        except Exception as e:
            log(f'    enrichr fail {t}: {e}')
    new_enr = pd.DataFrame(enr_rows)
    log(f'  new enr rows: {len(new_enr)}')

    # merge with existing
    existing = pd.read_csv(FIG/'myeloid_go_enrichment.csv')
    log(f'  existing rows: {len(existing)} (subtypes: {sorted(existing.subtype.unique())})')
    combined = pd.concat([existing, new_enr], ignore_index=True)
    # de-dup by (subtype, gene_set, term) keeping first (existing pre-computed wins for orig)
    combined = combined.drop_duplicates(subset=['subtype','gene_set','term'])
    log(f'  combined rows: {len(combined)} (subtypes: {sorted(combined.subtype.unique())})')

    # full and top10 versions
    combined.to_csv(RES/'step10_myeloid_go_enrichment_full.csv', index=False)
    top10 = (combined.sort_values('combined_score', ascending=False, na_position='last')
             .groupby(['subtype','gene_set']).head(10).reset_index(drop=True))
    top10.to_csv(RES/'step10_myeloid_go_enrichment.csv', index=False)
    top10.to_csv(FIG/'myeloid_go_enrichment.csv', index=False)
    log(f'  fig4 csv updated. {top10.subtype.nunique()} subtypes covered')
    log(f'DONE in {time.time()-t0:.1f}s')

if __name__ == '__main__':
    main()
