#!/usr/bin/env python
"""
Fig.1 data export
Export CSV files for R plotting from luad_integrated.h5ad

Input: ~/luad/data/processed/luad_integrated.h5ad (read-only)
Output: ${WORK_ROOT}/luad_figures/fig1/*.csv

Usage: conda activate scst && python data_prep/plot_data_exports/export_fig1_data.py
"""

import gc
import os
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

IN_PATH = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
OUT_DIR = Path("${WORK_ROOT}/luad_figures/fig1")

# ── Fig 1B: dot plot markers ──
# Aligned to HCC Fig.1B in the source paper; LUAD version
DOTPLOT_MARKERS = {
    "Epithelial":  ["EPCAM", "KRT18", "KRT19"],
    "Epithelial_prolif": ["TOP2A", "MKI67"],
    "Endothelial": ["PECAM1", "CDH5"],
    "Fibroblast":  ["COL1A1", "COL1A2", "ACTA2"],
    "T_NK":        ["CD3D", "CD3E", "NKG7", "KLRD1"],
    "B":           ["CD79A", "MS4A1"],
    "Plasma":      ["MZB1", "IGHG1"],
    "Myeloid":     ["LYZ", "CD68", "CD163", "CD14"],
    "Mast":        ["KIT", "TPSAB1", "CPA3"],
}

# Markers missing from the inner join; take from obsm['marker_counts']
# Known missing: VWF, NKX2-1 (TTF1), etc.
OBSM_FALLBACK_MARKERS = ["VWF", "NKX2-1", "SFTPC"]


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    log("=" * 60)
    log("Fig.1 data export")
    log("=" * 60)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log(f"Reading {IN_PATH}")
    adata = sc.read_h5ad(IN_PATH, backed="r")  # backed mode saves memory
    log(f"  shape: {adata.shape}")

    # Extract obs (full load into memory)
    obs = adata.obs.to_frame() if hasattr(adata.obs, "to_frame") else adata.obs.copy()
    obs = pd.DataFrame(obs)

    # UMAP
    umap = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    )

    # 1. cell_metadata.csv.gz — obs + UMAP (shared by Fig 1C/D/E/F)
    log("Export cell_metadata.csv.gz")
    meta_cols = [
        "dataset", "sample_id", "patient_id", "tissue_type",
        "celltype_coarse", "celltype_confidence",
        "celltype_celltypist", "celltype_ct_coarse",
        "leiden_1.0", "leiden_2.0",
    ]
    # Keep columns that exist
    meta_cols = [c for c in meta_cols if c in obs.columns]
    meta = obs[meta_cols].copy()
    meta["UMAP_1"] = umap["UMAP_1"].values
    meta["UMAP_2"] = umap["UMAP_2"].values
    meta.to_csv(OUT_DIR / "cell_metadata.csv.gz", compression="gzip")
    log(f"  {len(meta):,} cells × {len(meta.columns)} cols")

    # 2. dotplot_markers.csv — Fig 1B
    log("Export dotplot_markers.csv")

    all_markers = []
    for genes in DOTPLOT_MARKERS.values():
        all_markers.extend(genes)
    all_markers = list(dict.fromkeys(all_markers))  # dedupe keeping order

    # Check which are in var_names
    var_set = set(adata.var_names)
    in_var = [g for g in all_markers if g in var_set]
    not_in_var = [g for g in all_markers if g not in var_set]
    log(f"In var_names: {len(in_var)}/{len(all_markers)}")
    if not_in_var:
        log(f"Missing: {not_in_var}")

    # Expression from lognorm layer
    celltype = obs["celltype_coarse"].values
    ct_unique = sorted(obs["celltype_coarse"].unique())

    rows = []
    for gene in in_var:
        gene_idx = list(adata.var_names).index(gene)
        # backed mode requires column-wise access
        expr = adata.X[:, gene_idx]
        if sp.issparse(expr):
            expr = np.asarray(expr.todense()).flatten()
        else:
            expr = np.asarray(expr).flatten()

        # Also take lognorm version for mean expression
        expr_ln = adata.layers["lognorm"][:, gene_idx]
        if sp.issparse(expr_ln):
            expr_ln = np.asarray(expr_ln.todense()).flatten()
        else:
            expr_ln = np.asarray(expr_ln).flatten()

        for ct in ct_unique:
            mask = celltype == ct
            n_cells = mask.sum()
            if n_cells == 0:
                continue
            vals_raw = expr[mask]
            vals_ln = expr_ln[mask]
            frac = (vals_raw > 0).mean()
            mean_expr = vals_ln.mean()
            rows.append({
                "gene": gene,
                "celltype": ct,
                "frac_expressed": round(float(frac), 4),
                "mean_expression": round(float(mean_expr), 4),
                "n_cells": int(n_cells),
            })

    # Fill missing markers from obsm['marker_counts']
    if "marker_counts" in adata.obsm:
        marker_df = adata.obsm["marker_counts"]
        if isinstance(marker_df, pd.DataFrame):
            obsm_available = [g for g in not_in_var + OBSM_FALLBACK_MARKERS
                              if g in marker_df.columns]
            if obsm_available:
                log(f"Supplemented from obsm['marker_counts']: {obsm_available}")
                for gene in obsm_available:
                    vals = marker_df[gene].values.astype(np.float32)
                    # marker_counts are raw counts; log-normalize manually
                    # Row-wise total normalization (approximate: marker subset ≠ all genes)
                    # Simplified: use log1p(raw) as approximation
                    vals_ln = np.log1p(vals)
                    for ct in ct_unique:
                        mask = celltype == ct
                        n_cells = mask.sum()
                        if n_cells == 0:
                            continue
                        v = vals[mask]
                        v_ln = vals_ln[mask]
                        frac = (v > 0).mean()
                        # Handle NaN
                        valid = ~np.isnan(v_ln)
                        mean_expr = v_ln[valid].mean() if valid.any() else 0.0
                        frac_val = float((v[~np.isnan(v)] > 0).mean()) if (~np.isnan(v)).any() else 0.0
                        rows.append({
                            "gene": gene,
                            "celltype": ct,
                            "frac_expressed": round(frac_val, 4),
                            "mean_expression": round(float(mean_expr), 4),
                            "n_cells": int(n_cells),
                        })

    dotplot_df = pd.DataFrame(rows)
    # Add marker_group column
    gene_to_group = {}
    for group, genes in DOTPLOT_MARKERS.items():
        for g in genes:
            gene_to_group[g] = group
    for g in OBSM_FALLBACK_MARKERS:
        if g not in gene_to_group:
            gene_to_group[g] = "Supplementary"
    dotplot_df["marker_group"] = dotplot_df["gene"].map(gene_to_group)
    dotplot_df.to_csv(OUT_DIR / "dotplot_markers.csv", index=False)
    log(f"  {len(dotplot_df)} rows")

    # 3. proportion_by_tissue.csv — Fig 1D + 1E
    log("Export proportion_by_tissue.csv")

    ct_tissue = pd.crosstab(obs["tissue_type"], obs["celltype_coarse"])
    ct_tissue_pct = ct_tissue.div(ct_tissue.sum(axis=1), axis=0) * 100

    # Long format
    prop_rows = []
    for tissue in ct_tissue.index:
        for ct_col in ct_tissue.columns:
            prop_rows.append({
                "tissue_type": tissue,
                "celltype": ct_col,
                "count": int(ct_tissue.loc[tissue, ct_col]),
                "percent": round(float(ct_tissue_pct.loc[tissue, ct_col]), 2),
            })
    prop_df = pd.DataFrame(prop_rows)
    prop_df.to_csv(OUT_DIR / "proportion_by_tissue.csv", index=False)
    log(f"  {len(prop_df)} rows ({ct_tissue.shape[0]} tissues × {ct_tissue.shape[1]} celltypes)")

    # 4. proportion_by_sample.csv — Fig 1F (raw data for correlation matrix)
    log("Export proportion_by_sample.csv")

    ct_sample = pd.crosstab(obs["sample_id"], obs["celltype_coarse"])
    ct_sample_pct = ct_sample.div(ct_sample.sum(axis=1), axis=0)
    ct_sample_pct.to_csv(OUT_DIR / "proportion_by_sample.csv")
    log(f"  {ct_sample_pct.shape[0]} samples × {ct_sample_pct.shape[1]} celltypes")

    # 5. correlation_matrix.csv — Fig 1F
    log("Export correlation_matrix.csv")

    corr = ct_sample_pct.corr(method="pearson")
    corr.to_csv(OUT_DIR / "correlation_matrix.csv")
    log(f"  {corr.shape[0]}×{corr.shape[1]} Pearson correlation")

    # Also export p-value matrix (for R plot stars)
    from scipy.stats import pearsonr
    n_ct = len(corr.columns)
    pval_mat = pd.DataFrame(
        np.ones((n_ct, n_ct)),
        index=corr.index,
        columns=corr.columns,
    )
    for i in range(n_ct):
        for j in range(i + 1, n_ct):
            _, p = pearsonr(ct_sample_pct.iloc[:, i], ct_sample_pct.iloc[:, j])
            pval_mat.iloc[i, j] = p
            pval_mat.iloc[j, i] = p
    pval_mat.to_csv(OUT_DIR / "correlation_pvalues.csv")
    log(f"p-value matrix exported")

    # 6. enrichment_heatmap.csv — Fig 1D (Ro/e log2FC)
    log("Export enrichment_heatmap.csv")

    # Compute observed/expected ratio
    total = ct_tissue.values.sum()
    row_totals = ct_tissue.sum(axis=1).values  # per tissue
    col_totals = ct_tissue.sum(axis=0).values  # per celltype
    expected = np.outer(row_totals, col_totals) / total
    observed = ct_tissue.values.astype(float)

    # log2(O/E) with pseudocount to avoid log(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2_oe = np.log2((observed + 1) / (expected + 1))

    enrichment = pd.DataFrame(
        log2_oe,
        index=ct_tissue.index,
        columns=ct_tissue.columns,
    )
    enrichment.to_csv(OUT_DIR / "enrichment_heatmap.csv")
    log(f"  {enrichment.shape[0]}×{enrichment.shape[1]} log2(O/E)")

    # Significance (chi-squared per cell)
    from scipy.stats import chi2
    chi2_vals = ((observed - expected) ** 2) / (expected + 1)
    pvals_enrichment = 1 - chi2.cdf(chi2_vals, df=1)
    pval_enrich_df = pd.DataFrame(
        pvals_enrichment,
        index=ct_tissue.index,
        columns=ct_tissue.columns,
    )
    pval_enrich_df.to_csv(OUT_DIR / "enrichment_pvalues.csv")
    log(f"enrichment p-value matrix exported")

    # 7. cluster_celltype_summary.csv — for auxiliary checks
    log("Export cluster_celltype_summary.csv")
    if "leiden_1.0" in obs.columns:
        cluster_ct = pd.crosstab(obs["leiden_1.0"], obs["celltype_coarse"])
        cluster_ct.to_csv(OUT_DIR / "cluster_celltype_summary.csv")
        log(f"  {cluster_ct.shape[0]} clusters × {cluster_ct.shape[1]} celltypes")

    # Summary
    log("=" * 60)
    log("Export done; file list:")
    for f in sorted(OUT_DIR.iterdir()):
        size_mb = f.stat().st_size / 1e6
        log(f"  {f.name}: {size_mb:.1f} MB")
    log("=" * 60)
    log("R plotting path: ${WORK_ROOT}/luad_figures/fig1/")


if __name__ == "__main__":
    main()
