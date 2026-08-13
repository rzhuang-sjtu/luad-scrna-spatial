#!/usr/bin/env python3
"""
Step 3a: Extract raw counts for the extended marker panel from the 7 source h5ad files,
attach as obsm['marker_counts'] on the merged AnnData. Avoids genes lost by inner-join.

Usage:
  python 03a_extract_marker_layer.py --dry-run   # hit-rate diagnostics only
  python 03a_extract_marker_layer.py             # full extract + write new h5ad
"""
import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

DATA_DIR = Path("${WORK_ROOT}/数据清洗")
MERGED_PATH = Path("${PROJECT_ROOT}/data/processed/luad_merged_raw.h5ad")
OUT_PATH = Path("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad")

DATASETS = {
    "GSE123902": "GSE123902_clean.h5ad",
    "GSE131907": "GSE131907_clean.h5ad",
    "GSE143423": "GSE143423_LUAD_clean.h5ad",
    "GSE148071": "GSE148071_LUAD_clean.h5ad",
    "GSE164789": "GSE164789_LUAD_clean.h5ad",
    "GSE189357": "GSE189357_LUAD_clean.h5ad",
    "GSE253013": "GSE253013_clean.h5ad",
}

EXTENDED_MARKERS = {
    "Epithelial": ["EPCAM", "KRT18", "KRT19", "KRT8", "CDH1"],
    "Epithelial_AT2": ["SFTPC", "SFTPB", "SFTPD", "SFTPA1", "SFTPA2", "ABCA3", "LAMP3"],
    "Epithelial_AT1": ["AGER", "PDPN", "HOPX", "CAV1"],
    "Epithelial_Club": ["SCGB1A1", "SCGB3A2"],
    "Epithelial_Ciliated": ["FOXJ1", "TPPP3", "PIFO"],
    "Epithelial_LUAD": ["NKX2-1", "NAPSA"],
    "Endothelial": ["VWF", "PECAM1", "CDH5", "CLDN5", "ENG"],
    "Endothelial_Cap": ["CA4", "EDNRB"],
    "Endothelial_Tumor": ["ESM1", "ACKR1", "PTTG1"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "Fibroblast_CAF": ["ACTA2", "POSTN", "FAP", "CD36"],
    "Fibroblast_ap": ["CXCL12", "HLA-DRA"],
    "Myeloid": ["CD68", "CD163", "LYZ", "CD14", "C1QC", "C1QA", "C1QB"],
    "Myeloid_SPP1": ["SPP1", "TREM2"],
    "Myeloid_CXCL9": ["CXCL9", "CXCL10", "CXCL11"],
    "Myeloid_Alveolar": ["MARCO", "FABP4", "PPARG"],
    "Myeloid_DC": ["CLEC9A", "CLEC10A", "CD1C", "LAMP3"],
    "T_NK": ["CD3D", "CD3E", "CD3G", "NKG7", "KLRD1", "GNLY"],
    "T_CD8": ["CD8A", "CD8B"],
    "T_CD4": ["CD4", "IL7R"],
    "T_Treg": ["FOXP3", "CTLA4", "IL2RA"],
    "T_exh": ["HAVCR2", "LAG3", "PDCD1", "TIGIT"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Plasma": ["MZB1", "IGHG1", "IGKC", "JCHAIN", "DERL3", "XBP1"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3"],
}


def _flatten_panel(panel):
    """Flatten keeping first-occurrence order (LAMP3 counted once)"""
    seen, ordered = set(), []
    for genes in panel.values():
        for g in genes:
            if g not in seen:
                ordered.append(g); seen.add(g)
    return ordered


def _is_integer_matrix(x, n=200):
    if sparse.issparse(x):
        s = x[:min(n, x.shape[0])].toarray()
    else:
        s = np.asarray(x[:min(n, x.shape[0])])
    if s.size == 0:
        return False
    return np.allclose(s, s.astype(int))


def _get_raw_counts(adata):
    if _is_integer_matrix(adata.X):
        X = adata.X
    elif "counts" in adata.layers and _is_integer_matrix(adata.layers["counts"]):
        X = adata.layers["counts"]
    else:
        raise ValueError("Integer raw counts not found")
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    return X.tocsc()  # CSC speeds column slicing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    marker_list = _flatten_panel(EXTENDED_MARKERS)
    n_markers = len(marker_list)
    marker_to_col = {g: i for i, g in enumerate(marker_list)}
    print(f"Marker panel: {n_markers} unique genes ({len(EXTENDED_MARKERS)} categories)")

    # Read merged h5ad
    print(f"\nReading merged h5ad: {MERGED_PATH}")
    if args.dry_run:
        merged = sc.read_h5ad(MERGED_PATH, backed="r")
    else:
        merged = sc.read_h5ad(MERGED_PATH)
    print(f"  shape: {merged.shape}")
    n_cells = merged.n_obs

    merged_obs_names = merged.obs_names.to_numpy()
    merged_datasets = merged.obs["dataset"].astype(str).to_numpy()

    if not args.dry_run:
        marker_mat = np.full((n_cells, n_markers), np.nan, dtype=np.float32)

    diag = []

    for ds, f in DATASETS.items():
        ds_idx_in_merged = np.flatnonzero(merged_datasets == ds)
        n_ds_cells = len(ds_idx_in_merged)
        print(f"\n=== {ds} ({f}) ===")
        print(f"Cells for this dataset in merge: {n_ds_cells:,}")

        path = DATA_DIR / f
        adata = sc.read_h5ad(path, backed="r" if args.dry_run else None)
        print(f"Source shape: {adata.shape}")

        present = [g for g in marker_list if g in adata.var_names]
        missing = [g for g in marker_list if g not in adata.var_names]
        print(f"Marker hits: {len(present)}/{n_markers}")
        if missing:
            show = missing[:8]
            more = f"... (+{len(missing)-8} more)" if len(missing) > 8 else ""
            print(f"Missing: {show}{more}")

        diag.append({"dataset": ds, "n_cells": n_ds_cells,
                     "n_hit": len(present), "missing": list(missing)})

        if args.dry_run:
            if hasattr(adata, "file") and adata.file is not None:
                adata.file.close()
            del adata; gc.collect()
            continue

        # Strip prefix to recover original barcode
        prefix_len = len(ds) + 1
        merged_barcodes = np.array([s[prefix_len:] for s in merged_obs_names[ds_idx_in_merged]])

        orig_barcode_to_idx = {b: i for i, b in enumerate(adata.obs_names.astype(str))}
        try:
            orig_rows = np.array([orig_barcode_to_idx[b] for b in merged_barcodes])
        except KeyError as e:
            raise KeyError(f"{ds}: barcode {e} not found in source h5ad; merge naming may be inconsistent")

        counts = _get_raw_counts(adata)
        orig_var_idx = {v: i for i, v in enumerate(adata.var_names.astype(str))}
        present_var_cols = np.array([orig_var_idx[g] for g in present], dtype=np.int64)
        present_result_cols = np.array([marker_to_col[g] for g in present], dtype=np.int64)

        # Slice columns (fast on CSC), then rows
        sub = counts[:, present_var_cols].toarray()[orig_rows].astype(np.float32)
        marker_mat[np.ix_(ds_idx_in_merged, present_result_cols)] = sub

        del adata, counts, sub; gc.collect()

    # Summary
    print(f"\n{'='*60}\nDiagnostic summary\n{'='*60}")
    for d in diag:
        print(f"{d['dataset']}: cells={d['n_cells']:>7,}  marker={d['n_hit']:>3}/{n_markers}  missing={len(d['missing'])}")

    print(f"\nGSE123902 missing markers ({len([d for d in diag if d['dataset']=='GSE123902'][0]['missing'])}):")
    for g in [d for d in diag if d["dataset"] == "GSE123902"][0]["missing"]:
        print(f"  - {g}")

    if args.dry_run:
        print(f"\n[DRY RUN] No write. Estimated obsm['marker_counts'] shape = ({n_cells:,}, {n_markers})")
        return

    merged.obsm["marker_counts"] = pd.DataFrame(
        marker_mat, index=merged.obs_names.copy(), columns=marker_list
    )
    merged.uns["marker_panel"] = {k: list(v) for k, v in EXTENDED_MARKERS.items()}
    merged.uns["marker_panel_gene_order"] = list(merged.obsm["marker_counts"].columns)
    merged.uns["marker_panel_category_order"] = list(EXTENDED_MARKERS.keys())

    assert isinstance(merged.obsm["marker_counts"], pd.DataFrame), "obsm should be a DataFrame"
    assert list(merged.obsm["marker_counts"].columns) == merged.uns["marker_panel_gene_order"]

    nan_frac = np.isnan(marker_mat).sum() / marker_mat.size * 100
    print(f"\nobsm['marker_counts'] shape: {merged.obsm['marker_counts'].shape}")
    print(f"NaN fraction: {nan_frac:.2f}% (mainly missing columns from GSE123902)")
    print(f"uns['marker_panel'] categories: {len(merged.uns['marker_panel'])}")
    print(f"uns['marker_panel_gene_order'] length: {len(merged.uns['marker_panel_gene_order'])}")
    print(f"uns['marker_panel_category_order']: {merged.uns['marker_panel_category_order'][:5]}...")

    print(f"\nWriting {OUT_PATH} (gzip) ...", flush=True)
    merged.write_h5ad(OUT_PATH, compression="gzip")
    size_gb = OUT_PATH.stat().st_size / 1e9
    print(f"Done: {OUT_PATH} ({size_gb:.2f} GB)")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
