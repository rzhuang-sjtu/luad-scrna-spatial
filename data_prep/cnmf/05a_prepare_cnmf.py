"""
Step 5a: Prepare per-patient cNMF inputs.

Filters:
- Genes: drop MT-/RPS/RPL/MTRNR-prefixed + {MALAT1, NEAT1, XIST, KCNQ1OT1, MEG3}
  (mitochondrial / ribosomal / ubiquitous high-expression lncRNAs)
- Cells (scheme C): keep Malignant cells EXCEPT
    * celltype_marker in {T_NK, Myeloid, Mast, Plasma}  (clear immune by marker scoring)
    * doublet_score > 0.25
  Both filters apply as OR (exclude if either is true).

- Drop patients with < 50 cells after filter
- Export per-patient raw-counts h5ad + per-patient HVG
- Write run_order.txt sorted by cell count ascending
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

INPUT_H5AD = "${PROJECT_ROOT}/data/processed/luad_copykat.h5ad"
OUTPUT_ROOT = Path("${PROJECT_ROOT}/data/cnmf_input")
MIN_CELLS = 50
N_TOP_GENES = 3000

GENE_EXCLUDE_PREFIXES = ("MT-", "RPS", "RPL", "MTRNR")
GENE_EXCLUDE_EXACT = {"MALAT1", "NEAT1", "XIST", "KCNQ1OT1", "MEG3"}

IMMUNE_MARKER_LABELS = {"T_NK", "Myeloid", "Mast", "Plasma"}
DOUBLET_SCORE_THRESHOLD = 0.25


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def filter_genes(adata):
    gn = adata.var_names
    mask_pref = np.array([g.startswith(GENE_EXCLUDE_PREFIXES) for g in gn])
    mask_exact = np.array([g in GENE_EXCLUDE_EXACT for g in gn])
    drop = mask_pref | mask_exact
    kept = adata[:, ~drop].copy()
    log(f"Gene filter:")
    log(f"  before: {len(gn):,}")
    log(f"  MT-/RPS/RPL/MTRNR prefixes dropped: {mask_pref.sum():,}")
    log(f"  exact-match lncRNAs dropped: {mask_exact.sum():,}  "
        f"({[g for g in GENE_EXCLUDE_EXACT if g in set(gn)]})")
    log(f"  after:  {kept.shape[1]:,}")
    return kept


def filter_cells(adata):
    obs = adata.obs
    n0 = len(obs)
    mal = obs["malignant"].astype(str) == "Malignant"
    n_mal = int(mal.sum())

    immune_hit = obs["celltype_marker"].astype(str).isin(IMMUNE_MARKER_LABELS)
    doublet_hit = obs["doublet_score"].astype(float) > DOUBLET_SCORE_THRESHOLD

    drop_mask = mal & (immune_hit | doublet_hit)
    keep_mask = mal & ~(immune_hit | doublet_hit)
    n_drop_immune = int((mal & immune_hit & ~doublet_hit).sum())
    n_drop_doublet = int((mal & ~immune_hit & doublet_hit).sum())
    n_drop_both = int((mal & immune_hit & doublet_hit).sum())

    log(f"Cell filter (scheme C):")
    log(f"  Malignant initial: {n_mal:,}")
    log(f"  dropped: immune-marker only      : {n_drop_immune:,}")
    log(f"  dropped: doublet_score only      : {n_drop_doublet:,}")
    log(f"  dropped: immune AND doublet both : {n_drop_both:,}")
    log(f"  dropped total                    : {int(drop_mask.sum()):,}")
    log(f"  kept                             : {int(keep_mask.sum()):,}")

    # Breakdown of immune-marker drops by label
    if (mal & immune_hit).sum() > 0:
        log(f"  immune-marker-hit breakdown:")
        bd = (
            obs.loc[mal & immune_hit, "celltype_marker"]
            .astype(str)
            .value_counts()
        )
        for lab, n in bd.items():
            log(f"    {lab:12s} {int(n):6d}")

    return adata[keep_mask].copy()


def main():
    if OUTPUT_ROOT.exists():
        log(f"OUTPUT_ROOT {OUTPUT_ROOT} already exists; existing per-patient "
            "dirs will not be overwritten (counts.h5ad checked individually).")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    log(f"Reading {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    log(f"Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; cNMF requires raw counts")

    # Replace X with counts, drop lognorm layer to save memory after gene-filter
    adata = filter_genes(adata)

    adata = filter_cells(adata)
    log(f"Final subset: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    adata.obs["_patient_key"] = (
        adata.obs["dataset"].astype(str) + "__"
        + adata.obs["patient_id"].astype(str)
    )

    stats = (
        adata.obs.groupby("_patient_key", observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    meta = (
        adata.obs[["_patient_key", "dataset", "patient_id", "tissue_type"]]
        .drop_duplicates("_patient_key")
        .reset_index(drop=True)
    )
    stats = stats.merge(meta, on="_patient_key", how="left")
    stats = stats.sort_values("n_cells").reset_index(drop=True)

    log("=" * 70)
    log("Per-patient cell counts after filters (sorted ascending):")
    log("=" * 70)
    log(stats.to_string(index=False))
    log("=" * 70)

    qualified = stats[stats["n_cells"] >= MIN_CELLS].copy()
    dropped = stats[stats["n_cells"] < MIN_CELLS].copy()
    log(f"Qualified patients (>= {MIN_CELLS} cells): "
        f"{len(qualified)} / {len(stats)} (dropped {len(dropped)})")
    if len(dropped) > 0:
        log("Dropped patients:")
        log(dropped[["_patient_key", "n_cells"]].to_string(index=False))
    log(f"Total kept cells: {qualified['n_cells'].sum():,} / "
        f"{stats['n_cells'].sum():,}")

    stats.to_csv(OUTPUT_ROOT / "patient_stats.csv", index=False)
    log(f"Wrote {OUTPUT_ROOT / 'patient_stats.csv'}")

    order_path = OUTPUT_ROOT / "run_order.txt"
    with open(order_path, "w") as f:
        f.write("patient_key\tn_cells\tdataset\n")
        for _, r in qualified.iterrows():
            f.write(f"{r['_patient_key']}\t{int(r['n_cells'])}\t{r['dataset']}\n")
    log(f"Wrote {order_path}")

    log("Exporting per-patient raw-counts h5ad ...")
    n_done = 0
    for _, row in qualified.iterrows():
        key = row["_patient_key"]
        pdir = OUTPUT_ROOT / key
        pdir.mkdir(parents=True, exist_ok=True)
        out_h5 = pdir / "counts.h5ad"

        if out_h5.exists():
            n_done += 1
            continue

        sub = adata[adata.obs["_patient_key"] == key]
        X = sub.layers["counts"]
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        else:
            X = X.tocsr()
        if X.dtype != np.float32:
            X = X.astype(np.float32)

        patient_ad = ad.AnnData(
            X=X,
            obs=sub.obs[
                ["dataset", "patient_id", "sample_id", "tissue_type",
                 "celltype_coarse", "celltype_marker", "doublet_score"]
            ].copy(),
            var=sub.var[[]].copy(),
        )
        patient_ad.var_names = sub.var_names.copy()
        patient_ad.obs_names = sub.obs_names.copy()

        try:
            sc.pp.highly_variable_genes(
                patient_ad,
                n_top_genes=min(N_TOP_GENES, patient_ad.shape[1]),
                flavor="seurat_v3",
                inplace=True,
            )
        except Exception as e:
            log(f"  WARN [{key}] HVG failed ({e}); cNMF will re-select internally")

        patient_ad.write_h5ad(out_h5, compression="gzip")
        n_done += 1
        log(f"  [{n_done}/{len(qualified)}] {key}: "
            f"{patient_ad.shape[0]} cells -> {out_h5}")

    log(f"Done. Exported {n_done} patient directories under {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
