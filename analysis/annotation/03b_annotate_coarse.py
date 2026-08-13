#!/usr/bin/env python3
"""
Step 3: coarse cell-type annotation (three paths in parallel + voting)

Input: ${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad
Output: overwrite the same file (new obs columns)
Report: ${PROJECT_ROOT}/results/step3_annotation_summary.md

Usage:
  python 03b_annotate_coarse.py --path-1-only    # marker scoring only
  python 03b_annotate_coarse.py --skip-celltypist  # skip path 2
  python 03b_annotate_coarse.py                  # full three-path run
"""
import argparse
import gc
import sys
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

H5AD_PATH = Path("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad")
RESULTS_DIR = Path("${PROJECT_ROOT}/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CATEGORY_TO_COARSE = {
    # Epithelial lineage
    "Epithelial": "Epithelial",
    "Epithelial_AT2": "Epithelial",
    "Epithelial_AT1": "Epithelial",
    "Epithelial_Club": "Epithelial",
    "Epithelial_Ciliated": "Epithelial",
    "Epithelial_LUAD": "Epithelial",
    # Endothelial lineage
    "Endothelial": "Endothelial",
    "Endothelial_Cap": "Endothelial",
    "Endothelial_Tumor": "Endothelial",
    # Fibroblast lineage
    "Fibroblast": "Fibroblast",
    "Fibroblast_CAF": "Fibroblast",
    "Fibroblast_ap": "Fibroblast",
    # Myeloid lineage
    "Myeloid": "Myeloid",
    "Myeloid_SPP1": "Myeloid",
    "Myeloid_CXCL9": "Myeloid",
    "Myeloid_Alveolar": "Myeloid",
    "Myeloid_DC": "Myeloid",
    # T / NK lineage
    "T_NK": "T_NK",
    "T_CD8": "T_NK",
    "T_CD4": "T_NK",
    "T_Treg": "T_NK",
    "T_exh": "T_NK",
    # B lineage
    "B": "B",
    # Plasma
    "Plasma": "Plasma",
    # Mast
    "Mast": "Mast",
}

COARSE_CATEGORIES = ["Epithelial", "Endothelial", "Fibroblast",
                     "Myeloid", "Mast", "T_NK", "B", "Plasma"]

ORIGINAL_CELLTYPE_MAP = {
    # GSE131907 (Cell_type, Kim 2020 Nat Commun)
    "T lymphocytes": "T_NK",
    "NK cells": "T_NK",
    "B lymphocytes": "B",
    "Epithelial cells": "Epithelial",
    "Myeloid cells": "Myeloid",
    "Fibroblasts": "Fibroblast",
    "MAST cells": "Mast",
    "Endothelial cells": "Endothelial",
    "Oligodendrocytes": "Uncertain",  # contamination from brain-metastasis samples
    "Undetermined": "Uncertain",
    # GSE253013 (cell_type)
    "T cells": "T_NK",
    "Airway Epithelium": "Epithelial",
    "Myeloid": "Myeloid",
    "B cells": "B",
    "Fibroblasts": "Fibroblast",
    "Unknown": "Uncertain",
    "Endothelial": "Endothelial",
    "CD45+": "Uncertain",
    "Granulocytes": "Myeloid",  # map granulocytes to Myeloid
    "Epithelial": "Epithelial",
}

# This model has ~60 fine-grained labels; major classes are covered here
CELLTYPIST_MAP = {
    # Epithelial
    "AT1": "Epithelial", "AT2": "Epithelial",
    "Basal": "Epithelial", "Basal resting": "Epithelial",
    "Club": "Epithelial", "Goblet": "Epithelial",
    "Multiciliated": "Epithelial", "Multiciliated (non-nasal)": "Epithelial",
    "Deuterosomal": "Epithelial", "Ionocyte": "Epithelial",
    "Tuft": "Epithelial", "Mesothelium": "Epithelial",
    "SMG serous": "Epithelial", "SMG mucous": "Epithelial",
    "Neuroendocrine": "Epithelial", "pre-TB secretory": "Epithelial",
    "Suprabasal": "Epithelial",
    # Endothelial
    "EC aerocyte capillary": "Endothelial", "EC general capillary": "Endothelial",
    "EC venous pulmonary": "Endothelial", "EC venous systemic": "Endothelial",
    "EC arterial": "Endothelial", "Lymphatic EC differentiating": "Endothelial",
    "Lymphatic EC mature": "Endothelial", "Lymphatic EC proliferating": "Endothelial",
    # Fibroblast & stroma
    "Alveolar fibroblasts": "Fibroblast", "Adventitial fibroblasts": "Fibroblast",
    "Subpleural fibroblasts": "Fibroblast", "Peribronchial fibroblasts": "Fibroblast",
    "Myofibroblasts": "Fibroblast", "Pericytes": "Fibroblast",
    "Smooth muscle": "Fibroblast", "SM activated stress response": "Fibroblast",
    "Fibromyocyte": "Fibroblast", "Mesothelium": "Fibroblast",
    # Myeloid
    "Alveolar macrophages": "Myeloid", "Alveolar Mph CCL3+": "Myeloid",
    "Alveolar Mph MT-positive": "Myeloid", "Alveolar Mph proliferating": "Myeloid",
    "Classical monocytes": "Myeloid", "Non-classical monocytes": "Myeloid",
    "Monocyte-derived Mph": "Myeloid", "Interstitial Mph perivascular": "Myeloid",
    "DC1": "Myeloid", "DC2": "Myeloid", "Migratory DCs": "Myeloid",
    "Plasmacytoid DCs": "Myeloid",
    "Macrophages": "Myeloid",
    # Mast
    "Mast cells": "Mast",
    # T/NK
    "CD4 T cells": "T_NK", "CD8 T cells": "T_NK",
    "T cells proliferating": "T_NK",
    "NK cells": "T_NK", "NKT cells": "T_NK",
    "Tregs": "T_NK", "ILC3": "T_NK",
    # B
    "B cells": "B", "Naive B cells": "B", "Memory B cells": "B",
    "GC B cells": "B",
    # Plasma
    "Plasma cells": "Plasma", "Plasmablasts": "Plasma",
}


def path1_marker_scoring(adata):
    """
    Score each celltype category on obsm['marker_counts'].
    score_category_i = log1p(nanmean(marker_expr)) - log1p(total_counts / 1e4)
    """
    print("\n[Path 1] Marker scoring ...", flush=True)
    marker_df = adata.obsm["marker_counts"]   # DataFrame (n_cells, 93)
    marker_panel = adata.uns["marker_panel"]

    # Per-cell total counts (from X)
    total_per_cell = np.asarray(adata.X.sum(axis=1)).flatten()
    norm = np.log1p(total_per_cell / 1e4)

    # Score each category
    scores = {}
    for cat, gene_list in marker_panel.items():
        # Keep gene_list entries that are columns in marker_df
        present = [g for g in gene_list if g in marker_df.columns]
        if not present:
            print(f"category '{cat}' has no matching markers; skip")
            continue
        sub = marker_df[present].values  # (n_cells, k)
        mean_expr = np.nanmean(sub, axis=1)  # NaN skipped automatically
        mean_expr = np.nan_to_num(mean_expr, nan=0.0)  # all-NaN → 0
        score = np.log1p(mean_expr) - norm
        scores[cat] = score

    score_df = pd.DataFrame(scores, index=adata.obs_names)
    print(f"  25 category score matrix: {score_df.shape}")

    # Aggregate to 8 major classes (max score across sub-categories)
    coarse_scores = {}
    for coarse in COARSE_CATEGORIES:
        member_cats = [c for c, cg in CATEGORY_TO_COARSE.items() if cg == coarse]
        member_cats = [c for c in member_cats if c in score_df.columns]
        if not member_cats:
            continue
        coarse_scores[coarse] = score_df[member_cats].max(axis=1).values
    coarse_df = pd.DataFrame(coarse_scores, index=adata.obs_names)

    # argmax as final class
    top = coarse_df.idxmax(axis=1)
    top_score = coarse_df.max(axis=1)
    # second-highest score (borderline check)
    second_score = coarse_df.apply(
        lambda row: row.nlargest(2).iloc[1] if len(row) >= 2 else np.nan,
        axis=1
    )

    adata.obs["celltype_marker"] = top.values
    adata.obs["celltype_marker_score"] = top_score.values.astype("float32")
    adata.obs["celltype_marker_second_score"] = second_score.values.astype("float32")

    print(f"celltype_marker distribution:")
    print(adata.obs["celltype_marker"].value_counts())


def path2_celltypist(adata):
    """Run CellTypist Human_Lung_Atlas prediction"""
    print("\n[Path 2] CellTypist Human_Lung_Atlas ...", flush=True)
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        print("celltypist not installed; run pip install celltypist and retry")
        raise

    # Download / load model
    model_name = "Human_Lung_Atlas.pkl"
    try:
        model = models.Model.load(model=model_name)
    except Exception:
        print(f"Downloading model {model_name} (first time)...")
        models.download_models(model=[model_name], force_update=False)
        model = models.Model.load(model=model_name)

    # CellTypist requires log1p(normalize_total) data
    # Build a temporary copy (do not modify main adata)
    print("Preparing normalized copy ...")
    tmp = ad.AnnData(
        X=adata.X.copy(),
        obs=adata.obs[["dataset"]].copy(),
        var=pd.DataFrame(index=adata.var_names),
    )
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    print("Running CellTypist prediction (majority_voting=True; may take 10-15 min)...")
    pred = celltypist.annotate(
        tmp, model=model,
        majority_voting=True,
        mode="best match",
    )

    # predicted_labels is a DataFrame with predicted_labels / majority_voting / conf_score
    labels_df = pred.predicted_labels
    print(f"CellTypist output columns: {list(labels_df.columns)}")

    # Take majority_voting column
    if "majority_voting" in labels_df.columns:
        raw_labels = labels_df["majority_voting"].values
    else:
        raw_labels = labels_df["predicted_labels"].values

    # Confidence
    if "conf_score" in labels_df.columns:
        conf = labels_df["conf_score"].values
    else:
        conf = np.full(adata.n_obs, np.nan)

    # Map to 8 major classes
    mapped = pd.Series(raw_labels).map(CELLTYPIST_MAP).fillna("Uncertain").values

    adata.obs["celltype_celltypist_raw"] = raw_labels
    adata.obs["celltype_celltypist"] = mapped
    adata.obs["celltype_celltypist_prob"] = conf.astype("float32")

    print(f"celltype_celltypist distribution:")
    print(adata.obs["celltype_celltypist"].value_counts())

    # Report unmapped raw labels
    unmapped = set(raw_labels) - set(CELLTYPIST_MAP.keys())
    if unmapped:
        print(f"Unmapped CellTypist raw labels ({len(unmapped)}): {list(unmapped)[:20]}")

    del tmp; gc.collect()


def path3_original_mapping(adata):
    print("\n[Path 3] Align original author annotations ...", flush=True)
    orig = adata.obs["celltype_original"].astype(str)
    mapped = orig.map(ORIGINAL_CELLTYPE_MAP)
    # Empty string or unmapped → NaN
    mapped = mapped.where(orig != "", np.nan)
    mapped = mapped.where(mapped.notna(), np.nan)

    adata.obs["celltype_original_mapped"] = mapped.astype(str).replace("nan", "")

    n_covered = (mapped.notna() & (orig != "")).sum()
    print(f"Path 3 coverage: {n_covered:,}/{adata.n_obs:,} ({100*n_covered/adata.n_obs:.1f}%)")
    print(adata.obs["celltype_original_mapped"].value_counts())


def vote(adata, has_celltypist, has_original):
    print("\n[Voting] Combine three paths ...", flush=True)

    m = adata.obs["celltype_marker"].astype(str).values
    c = adata.obs["celltype_celltypist"].astype(str).values if has_celltypist else np.array([""] * adata.n_obs)
    o = adata.obs["celltype_original_mapped"].astype(str).values if has_original else np.array([""] * adata.n_obs)

    final = np.empty(adata.n_obs, dtype=object)
    conf = np.empty(adata.n_obs, dtype=object)
    n_votes_arr = np.zeros(adata.n_obs, dtype=np.int8)

    for i in range(adata.n_obs):
        votes = []
        for v in (m[i], c[i], o[i]):
            if v and v not in ("", "nan", "Uncertain"):
                votes.append(v)
        n_votes_arr[i] = len(votes)

        if len(votes) == 0:
            final[i] = "Uncertain"
            conf[i] = "none"
            continue

        counter = Counter(votes)
        top, top_count = counter.most_common(1)[0]

        if len(votes) == 1:
            final[i] = top
            conf[i] = "single-source"
        elif top_count == len(votes):
            # Full agreement
            final[i] = top
            conf[i] = "high" if len(votes) == 3 else "medium"
        elif top_count >= 2:
            # Majority
            final[i] = top
            conf[i] = "medium"
        else:
            # All three disagree: prefer marker (path 1)
            final[i] = m[i] if (m[i] and m[i] != "Uncertain") else "Uncertain"
            conf[i] = "low"

    adata.obs["celltype_coarse"] = final
    adata.obs["celltype_confidence"] = conf
    adata.obs["n_votes"] = n_votes_arr

    print(f"\ncelltype_coarse distribution:")
    print(adata.obs["celltype_coarse"].value_counts())
    print(f"\ncelltype_confidence distribution:")
    print(adata.obs["celltype_confidence"].value_counts())


def make_report(adata, has_celltypist, has_original):
    L = ["# Step 3 coarse cell-type annotation report", ""]
    L.append(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- Total cells: {adata.n_obs:,}")
    L.append(f"- Paths run: path1=  path2={'' if has_celltypist else ''}  path3={'' if has_original else ''}")
    L.append("")

    L.append("## Final celltype_coarse distribution")
    for ct, n in adata.obs["celltype_coarse"].value_counts().items():
        pct = 100 * n / adata.n_obs
        L.append(f"- {ct}: {n:,} ({pct:.2f}%)")
    L.append("")

    L.append("## celltype_confidence distribution")
    for c, n in adata.obs["celltype_confidence"].value_counts().items():
        L.append(f"- {c}: {n:,}")
    L.append("")

    L.append("## Three-path agreement matrix (path1 vs path2)")
    if has_celltypist:
        ct = pd.crosstab(adata.obs["celltype_marker"],
                          adata.obs["celltype_celltypist"])
        L.append("```")
        L.append(ct.to_string())
        L.append("```")
        diag = np.trace(ct.values) if ct.shape[0] == ct.shape[1] else np.nan
        total = ct.values.sum()
        if not np.isnan(diag):
            L.append(f"\nDiagonal agreement: {100*diag/total:.2f}%")
    L.append("")

    L.append("## Three-path agreement matrix (path1 vs path3; cells covered by path3 only)")
    if has_original:
        mask = adata.obs["celltype_original_mapped"].astype(str) != ""
        sub = adata.obs.loc[mask]
        ct = pd.crosstab(sub["celltype_marker"], sub["celltype_original_mapped"])
        L.append("```")
        L.append(ct.to_string())
        L.append("```")
    L.append("")

    L.append("## celltype_coarse × dataset crosstab")
    ct = pd.crosstab(adata.obs["dataset"], adata.obs["celltype_coarse"])
    L.append("```")
    L.append(ct.to_string())
    L.append("```")
    L.append("")

    L.append("## celltype_coarse × tissue_type crosstab")
    ct = pd.crosstab(adata.obs["tissue_type"], adata.obs["celltype_coarse"])
    L.append("```")
    L.append(ct.to_string())
    L.append("```")
    L.append("")

    L.append("## celltype_marker_score distribution by dataset (for QC)")
    L.append("```")
    q = adata.obs.groupby("dataset")["celltype_marker_score"].describe()
    L.append(q.to_string())
    L.append("```")

    return "\n".join(L)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-1-only", action="store_true",
                        help="Run path1 only; skip CellTypist and original annotations")
    parser.add_argument("--skip-celltypist", action="store_true",
                        help="Skip CellTypist (saves 10-15 min)")
    args = parser.parse_args()

    print("=" * 60)
    print(f"Step 3 coarse annotation started: {datetime.now()}")
    print("=" * 60)

    print(f"\nReading {H5AD_PATH} ...")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"shape: {adata.shape}")

    # Path 1 always runs
    path1_marker_scoring(adata)

    # Path 2
    has_celltypist = not (args.path_1_only or args.skip_celltypist)
    if has_celltypist:
        try:
            path2_celltypist(adata)
        except Exception as e:
            print(f"CellTypist failed: {e}; continue without path 2")
            has_celltypist = False

    # Path 3
    has_original = not args.path_1_only
    if has_original:
        path3_original_mapping(adata)

    # Voting
    vote(adata, has_celltypist, has_original)

    # Write h5ad first
    print(f"\nWriting {H5AD_PATH} ...")
    adata.write_h5ad(H5AD_PATH, compression="gzip")
    print(f"h5ad updated: {H5AD_PATH}")

    # Then write report
    try:
        report = make_report(adata, has_celltypist, has_original)
        report_path = RESULTS_DIR / "step3_annotation_summary.md"
        report_path.write_text(report, encoding="utf-8")
        print(f"Report: {report_path}")
    except Exception as e:
        print(f"Report generation failed but h5ad saved: {e}", file=sys.stderr)

    print(f"\nDone: {datetime.now()}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
