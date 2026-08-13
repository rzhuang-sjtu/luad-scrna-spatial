#!/usr/bin/env python3
"""
Step 3c: Leiden cluster-based coarse annotation (CellTypist primary + author anchors + marker check)

Input:  ${PROJECT_ROOT}/data/processed/luad_integrated.h5ad  (integrated + clustered)
Output:  overwrite the same file; new obs columns:
         celltype_celltypist       (fine-grained, CellTypist majority_voting)
         celltype_ct_coarse        (8 major classes, mapped from fine labels)
         celltype_original_mapped  (8 major classes, mapped from author labels; GSE131907/GSE253013 only)
         celltype_coarse           (final call, cluster-level)
         celltype_confidence       (high / medium)
       obsm['marker_scores_coarse'] = DataFrame (8 major-class scores)
Report: ${PROJECT_ROOT}/results/step3_annotation_summary.md
"""
import gc
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

IN_PATH = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
RESULTS_DIR = Path("${PROJECT_ROOT}/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
REPORT_PATH = RESULTS_DIR / "step3_annotation_summary.md"

CLUSTER_KEY = "leiden_1.0"

COARSE_CATEGORIES = ["Epithelial", "Endothelial", "Fibroblast",
                     "Myeloid", "Mast", "T_NK", "B", "Plasma"]

# marker category → 8 major classes (for Part3 marker-score aggregation)
CATEGORY_TO_COARSE = {
    "Epithelial": "Epithelial", "Epithelial_AT2": "Epithelial",
    "Epithelial_AT1": "Epithelial", "Epithelial_Club": "Epithelial",
    "Epithelial_Ciliated": "Epithelial", "Epithelial_LUAD": "Epithelial",
    "Endothelial": "Endothelial", "Endothelial_Cap": "Endothelial",
    "Endothelial_Tumor": "Endothelial",
    "Fibroblast": "Fibroblast", "Fibroblast_CAF": "Fibroblast",
    "Fibroblast_ap": "Fibroblast",
    "Myeloid": "Myeloid", "Myeloid_SPP1": "Myeloid",
    "Myeloid_CXCL9": "Myeloid", "Myeloid_Alveolar": "Myeloid",
    "Myeloid_DC": "Myeloid",
    "T_NK": "T_NK", "T_CD8": "T_NK", "T_CD4": "T_NK",
    "T_Treg": "T_NK", "T_exh": "T_NK",
    "B": "B",
    "Plasma": "Plasma",
    "Mast": "Mast",
}

# CellTypist Human_Lung_Atlas fine labels → 8 major classes
CELLTYPIST_MAP = {
    # Epithelial
    "AT0": "Epithelial", "AT1": "Epithelial", "AT2": "Epithelial",
    "AT2 proliferating": "Epithelial",
    "Basal": "Epithelial", "Basal resting": "Epithelial",
    "Suprabasal": "Epithelial",
    "Club": "Epithelial", "Club (non-nasal)": "Epithelial",
    "Club (nasal)": "Epithelial",
    "Goblet": "Epithelial", "Goblet (nasal)": "Epithelial",
    "Goblet (subsegmental)": "Epithelial", "Goblet (bronchial)": "Epithelial",
    "Multiciliated": "Epithelial", "Multiciliated (non-nasal)": "Epithelial",
    "Multiciliated (nasal)": "Epithelial",
    "Deuterosomal": "Epithelial", "Ionocyte": "Epithelial",
    "Tuft": "Epithelial", "Mesothelium": "Epithelial",
    "SMG serous (bronchial)": "Epithelial", "SMG serous (nasal)": "Epithelial",
    "SMG mucous": "Epithelial", "SMG duct": "Epithelial",
    "Neuroendocrine": "Epithelial", "pre-TB secretory": "Epithelial",
    "Hillock-like": "Epithelial",
    # Endothelial
    "EC aerocyte capillary": "Endothelial", "EC general capillary": "Endothelial",
    "EC capillary": "Endothelial",
    "EC venous pulmonary": "Endothelial", "EC venous systemic": "Endothelial",
    "EC venous": "Endothelial",
    "EC arterial": "Endothelial",
    "Lymphatic EC differentiating": "Endothelial",
    "Lymphatic EC mature": "Endothelial",
    "Lymphatic EC proliferating": "Endothelial",
    # Fibroblast / stroma
    "Alveolar fibroblasts": "Fibroblast",
    "Adventitial fibroblasts": "Fibroblast",
    "Subpleural fibroblasts": "Fibroblast",
    "Peribronchial fibroblasts": "Fibroblast",
    "Myofibroblasts": "Fibroblast",
    "Pericytes": "Fibroblast",
    "Smooth muscle": "Fibroblast",
    "Smooth muscle FAM83D+": "Fibroblast",
    "SM activated stress response": "Fibroblast",
    "Fibromyocyte": "Fibroblast",
    "Chondrocytes": "Fibroblast",
    # Myeloid
    "Alveolar macrophages": "Myeloid",
    "Alveolar Mph CCL3+": "Myeloid", "Alveolar Mph MT-positive": "Myeloid",
    "Alveolar Mph proliferating": "Myeloid",
    "Classical monocytes": "Myeloid", "Non-classical monocytes": "Myeloid",
    "Monocyte-derived Mph": "Myeloid",
    "Interstitial Mph perivascular": "Myeloid",
    "Macrophages": "Myeloid",
    "DC1": "Myeloid", "DC2": "Myeloid", "Migratory DCs": "Myeloid",
    "Plasmacytoid DCs": "Myeloid",
    # Mast
    "Mast cells": "Mast",
    # T / NK
    "CD4 T cells": "T_NK", "CD8 T cells": "T_NK",
    "T cells proliferating": "T_NK",
    "NK cells": "T_NK", "NKT cells": "T_NK",
    "Tregs": "T_NK", "ILC3": "T_NK", "ILC": "T_NK",
    # B
    "B cells": "B", "Naive B cells": "B", "Memory B cells": "B",
    "GC B cells": "B",
    # Plasma
    "Plasma cells": "Plasma", "Plasmablasts": "Plasma",
}

# Author celltype_original → 8 major classes
ORIGINAL_CELLTYPE_MAP = {
    # GSE131907 (Kim 2020)
    "T lymphocytes": "T_NK",
    "NK cells": "T_NK",
    "B lymphocytes": "B",
    "Epithelial cells": "Epithelial",
    "Myeloid cells": "Myeloid",
    "Fibroblasts": "Fibroblast",
    "MAST cells": "Mast",
    "Mast cells": "Mast",
    "Endothelial cells": "Endothelial",
    "Oligodendrocytes": "Uncertain",
    "Undetermined": "Uncertain",
    # GSE253013
    "T cells": "T_NK",
    "Airway Epithelium": "Epithelial",
    "Epithelial": "Epithelial",
    "Myeloid": "Myeloid",
    "B cells": "B",
    "Endothelial": "Endothelial",
    "Granulocytes": "Myeloid",
    "Unknown": "Uncertain",
    "CD45+": "Uncertain",
}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def part1_celltypist(adata):
    log("=" * 60)
    log("Part 1: CellTypist Human_Lung_Atlas")
    log("=" * 60)
    import celltypist
    from celltypist import models

    model_name = "Human_Lung_Atlas.pkl"
    try:
        model = models.Model.load(model=model_name)
        log(f"Loaded model {model_name}")
    except Exception:
        log(f"Downloading model {model_name} (first time)")
        models.download_models(model=[model_name], force_update=False)
        model = models.Model.load(model=model_name)

    log("Building temporary AnnData (X = lognorm)")
    tmp = ad.AnnData(
        X=adata.layers["lognorm"].copy(),
        obs=adata.obs[[CLUSTER_KEY]].copy(),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )

    log(f"  celltypist.annotate (majority_voting=True, over_clustering='{CLUSTER_KEY}')")
    pred = celltypist.annotate(
        tmp, model=model, majority_voting=True, over_clustering=CLUSTER_KEY,
    )
    labels_df = pred.predicted_labels
    log(f"Output columns: {list(labels_df.columns)}")

    # Fine labels (prefer majority_voting)
    if "majority_voting" in labels_df.columns:
        raw_labels = labels_df["majority_voting"].astype(str).values
    else:
        raw_labels = labels_df["predicted_labels"].astype(str).values

    # Print all unique labels for review
    uniq = sorted(set(raw_labels))
    log(f"Unique CellTypist fine labels: {len(uniq)}")
    for u in uniq:
        n = (raw_labels == u).sum()
        mapped = CELLTYPIST_MAP.get(u, "Uncertain")
        log(f"    {u!r} → {mapped}  ({n:,} cells)")

    unmapped = set(uniq) - set(CELLTYPIST_MAP.keys())
    if unmapped:
        log(f"Labels not in the mapping table ({len(unmapped)}) will be marked Uncertain:")
        for u in sorted(unmapped):
            log(f"    - {u!r}")

    adata.obs["celltype_celltypist"] = pd.Categorical(raw_labels)
    ct_coarse = np.array([CELLTYPIST_MAP.get(u, "Uncertain") for u in raw_labels])
    adata.obs["celltype_ct_coarse"] = pd.Categorical(
        ct_coarse, categories=COARSE_CATEGORIES + ["Uncertain"]
    )

    log("celltype_ct_coarse distribution:")
    for ct, n in adata.obs["celltype_ct_coarse"].value_counts().items():
        log(f"    {ct}: {n:,}")

    del tmp, pred
    gc.collect()


def part2_original(adata):
    log("=" * 60)
    log("Part 2: Map author annotations")
    log("=" * 60)
    orig = adata.obs["celltype_original"].astype(str).fillna("")
    orig = orig.where(~orig.isin(["", "nan", "None", "NaN"]), "")

    # Print original label distribution per dataset
    for ds in adata.obs["dataset"].unique():
        mask = adata.obs["dataset"] == ds
        sub = orig[mask]
        sub_non_empty = sub[sub != ""]
        log(f"{ds}: annotated {len(sub_non_empty):,}/{mask.sum():,}")
        if len(sub_non_empty) > 0:
            for lab, n in sub_non_empty.value_counts().items():
                tgt = ORIGINAL_CELLTYPE_MAP.get(lab, "?unmapped?")
                log(f"    {lab!r} → {tgt}  ({n:,})")

    mapped_vals = []
    for v in orig.values:
        if v == "":
            mapped_vals.append("")
        else:
            mapped_vals.append(ORIGINAL_CELLTYPE_MAP.get(v, "Uncertain"))
    adata.obs["celltype_original_mapped"] = pd.Categorical(
        mapped_vals,
        categories=[""] + COARSE_CATEGORIES + ["Uncertain"],
    )

    n_covered = (adata.obs["celltype_original_mapped"].astype(str) != "").sum()
    log(f"celltype_original_mapped coverage: {n_covered:,}/{adata.n_obs:,}"
        f"({100*n_covered/adata.n_obs:.1f}%)")
    vc = adata.obs["celltype_original_mapped"].value_counts()
    for ct, n in vc.items():
        tag = "(no annotation)" if ct == "" else ""
        log(f"    {ct!r} {tag}: {n:,}")


def part3_markers(adata):
    log("=" * 60)
    log("Part 3: Marker scoring (sc.tl.score_genes on lognorm)")
    log("=" * 60)

    # Take marker_panel; correct Fibroblast_ap
    panel = {k: list(v) for k, v in adata.uns["marker_panel"].items()}
    if "Fibroblast_ap" in panel:
        before = panel["Fibroblast_ap"]
        panel["Fibroblast_ap"] = [g for g in before if g != "HLA-DRA"]
        log(f"Corrected Fibroblast_ap: {before} → {panel['Fibroblast_ap']}")

    # Temporarily set X to lognorm for score_genes
    x_backup = adata.X
    adata.X = adata.layers["lognorm"]

    score_cols = {}
    try:
        for cat, genes in panel.items():
            present = [g for g in genes if g in adata.var_names]
            if not present:
                log(f"{cat}: no matching genes; skip")
                continue
            score_key = f"_score_{cat}"
            sc.tl.score_genes(
                adata, gene_list=present, score_name=score_key,
                use_raw=False, random_state=0,
            )
            score_cols[cat] = adata.obs[score_key].values.astype(np.float32)
            del adata.obs[score_key]
        log(f"Finished scoring {len(score_cols)} categories")
    finally:
        adata.X = x_backup

    # Aggregate to 8 major classes: top-2 mean of member categories per class
    score_df = pd.DataFrame(score_cols, index=adata.obs_names)
    coarse_scores = {}
    for coarse in COARSE_CATEGORIES:
        members = [c for c, tgt in CATEGORY_TO_COARSE.items()
                   if tgt == coarse and c in score_df.columns]
        if not members:
            continue
        sub = score_df[members].values
        if sub.shape[1] >= 2:
            top2 = np.sort(sub, axis=1)[:, -2:]
            coarse_scores[coarse] = top2.mean(axis=1)
        else:
            coarse_scores[coarse] = sub[:, 0]
    coarse_df = pd.DataFrame(coarse_scores, index=adata.obs_names).astype(np.float32)
    adata.obsm["marker_scores_coarse"] = coarse_df
    log(f"  obsm['marker_scores_coarse'] shape: {coarse_df.shape}")

    # Print top-scoring coarse label per cluster (for review)
    log("Top-scoring coarse label per cluster (review):")
    clusters = adata.obs[CLUSTER_KEY].astype(str)
    by_cluster = coarse_df.groupby(clusters).mean()
    argmax_ct = by_cluster.idxmax(axis=1)
    max_score = by_cluster.max(axis=1)
    for cl in sorted(argmax_ct.index, key=lambda s: int(s) if s.isdigit() else s):
        n = (clusters == cl).sum()
        log(f"    cluster {cl} (n={n:,}): marker→{argmax_ct[cl]} (score={max_score[cl]:.3f})")


def part4_decide(adata):
    log("=" * 60)
    log("Part 4: Cluster-level decision")
    log("=" * 60)
    clusters = adata.obs[CLUSTER_KEY].astype(str).values
    ct_coarse = adata.obs["celltype_ct_coarse"].astype(str).values
    orig_mapped = adata.obs["celltype_original_mapped"].astype(str).values

    final_map = {}
    conf_map = {}
    records = []

    for cl in pd.unique(clusters):
        mask = clusters == cl
        n = int(mask.sum())

        # CellTypist mode (ignore Uncertain)
        ct_vals = ct_coarse[mask]
        ct_valid = ct_vals[ct_vals != "Uncertain"]
        ct_label = Counter(ct_valid).most_common(1)[0][0] if len(ct_valid) > 0 else "Uncertain"

        # Author mode (ignore "" and Uncertain)
        orig_vals = orig_mapped[mask]
        orig_valid = orig_vals[(orig_vals != "") & (orig_vals != "Uncertain")]
        if len(orig_valid) > 0:
            orig_counter = Counter(orig_valid)
            orig_label, orig_count = orig_counter.most_common(1)[0]
            orig_frac = orig_count / n   # relative to the whole cluster
        else:
            orig_label = None
            orig_frac = 0.0

        # Decision
        if orig_label is not None and ct_label == orig_label:
            final = ct_label
            conf = "high"
        elif orig_label is not None and ct_label != orig_label:
            if orig_frac > 0.6:
                final = orig_label
                conf = "medium"
            else:
                final = ct_label
                conf = "medium"
        else:
            final = ct_label
            conf = "medium"

        final_map[cl] = final
        conf_map[cl] = conf

        records.append({
            "cluster": cl, "n_cells": n,
            "ct_label": ct_label,
            "orig_label": orig_label if orig_label is not None else "",
            "orig_frac": round(orig_frac, 3),
            "final": final, "conf": conf,
        })

    cluster_df = pd.DataFrame(records)
    try:
        cluster_df["cluster"] = cluster_df["cluster"].astype(int)
        cluster_df = cluster_df.sort_values("cluster").reset_index(drop=True)
        cluster_df["cluster"] = cluster_df["cluster"].astype(str)
    except ValueError:
        cluster_df = cluster_df.sort_values("cluster").reset_index(drop=True)

    log("Cluster decision table (by cluster order):")
    log("  " + cluster_df.to_string(index=False).replace("\n", "\n  "))

    final_arr = np.array([final_map[c] for c in clusters])
    conf_arr = np.array([conf_map[c] for c in clusters])
    adata.obs["celltype_coarse"] = pd.Categorical(
        final_arr, categories=COARSE_CATEGORIES + ["Uncertain"]
    )
    adata.obs["celltype_confidence"] = pd.Categorical(
        conf_arr, categories=["high", "medium", "low"]
    )

    adata.uns["cluster_decision_table"] = cluster_df

    log("Final celltype_coarse distribution:")
    for ct, n in adata.obs["celltype_coarse"].value_counts().items():
        log(f"    {ct}: {n:,}  ({100*n/adata.n_obs:.2f}%)")
    log("celltype_confidence distribution:")
    for c, n in adata.obs["celltype_confidence"].value_counts().items():
        log(f"    {c}: {n:,}")

    return cluster_df


EXPECTED_RANGES = {
    "Epithelial":  (25, 40),
    "T_NK":        (25, 40),
    "Myeloid":     (10, 25),
    "Fibroblast":  (5, 8),
    "B":           (3, 6),
    "Plasma":      (1, 5),
    "Mast":        (0.5, 3),
    "Endothelial": (3, 8),
    "Uncertain":   (0, 5),
}


def part5_report(adata, cluster_df):
    log("=" * 60)
    log("Part 5: Write report")
    log("=" * 60)

    L = ["# Step 3 coarse annotation report (cluster-level)", ""]
    L.append(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- Total cells: {adata.n_obs:,}")
    L.append(f"- cluster key: `{CLUSTER_KEY}` ({adata.obs[CLUSTER_KEY].nunique()} clusters)")
    L.append("")

    # Final distribution vs expected ranges
    L.append("## celltype_coarse distribution vs manuscript expectations")
    L.append("")
    L.append("| CellType | Count | % | Expected range | Status |")
    L.append("|---|---|---|---|---|")
    vc = adata.obs["celltype_coarse"].value_counts()
    for ct in COARSE_CATEGORIES + ["Uncertain"]:
        n = int(vc.get(ct, 0))
        pct = 100 * n / adata.n_obs
        lo, hi = EXPECTED_RANGES.get(ct, (None, None))
        if lo is None:
            status = "-"
        elif lo <= pct <= hi:
            status = ""
        elif pct < lo:
            status = f"↓ (below {lo}%)"
        else:
            status = f"↑ (above {hi}%)"
        rng = f"{lo}-{hi}%" if lo is not None else "-"
        L.append(f"| {ct} | {n:,} | {pct:.2f}% | {rng} | {status} |")
    L.append("")

    L.append("## celltype_confidence distribution")
    for c, n in adata.obs["celltype_confidence"].value_counts().items():
        L.append(f"- {c}: {n:,} ({100*n/adata.n_obs:.2f}%)")
    L.append("")

    # CellTypist vs author conflicts
    L.append("## CellTypist vs author annotation conflicts (cluster decision table only)")
    cd = cluster_df.copy()
    has_orig = cd["orig_label"] != ""
    agree = has_orig & (cd["ct_label"] == cd["orig_label"])
    disagree = has_orig & (cd["ct_label"] != cd["orig_label"])
    L.append(f"- Clusters with author labels: {has_orig.sum()} / {len(cd)}")
    L.append(f"- Agree: {agree.sum()}")
    L.append(f"- Conflict: {disagree.sum()}")
    if disagree.any():
        L.append("")
        L.append("### Conflicting clusters")
        L.append("```")
        L.append(cd[disagree].to_string(index=False))
        L.append("```")
    L.append("")

    # Per-cluster decision table
    L.append("## Per-cluster decision details")
    L.append("```")
    L.append(cluster_df.to_string(index=False))
    L.append("```")
    L.append("")

    # dataset × coarse crosstab
    L.append("## celltype_coarse × dataset crosstab")
    ct = pd.crosstab(adata.obs["dataset"], adata.obs["celltype_coarse"])
    L.append("```")
    L.append(ct.to_string())
    L.append("```")
    L.append("")

    # tissue × coarse
    L.append("## celltype_coarse × tissue_type crosstab")
    ct = pd.crosstab(adata.obs["tissue_type"], adata.obs["celltype_coarse"])
    L.append("```")
    L.append(ct.to_string())
    L.append("```")

    report = "\n".join(L)
    REPORT_PATH.write_text(report, encoding="utf-8")
    log(f"Report: {REPORT_PATH}")
    # Also print to stdout
    print("\n" + "=" * 60)
    print(report)
    print("=" * 60)


def main():
    t0 = datetime.now()
    log(f"Step 3c Annotate started: {t0}")

    log(f"Reading {IN_PATH}")
    adata = sc.read_h5ad(IN_PATH)
    log(f"  shape: {adata.shape}")
    log(f"  layers: {list(adata.layers.keys())}")
    log(f"  obs cols: {list(adata.obs.columns)}")
    assert "lognorm" in adata.layers, "Requires layers['lognorm'] (from 03b)"
    assert CLUSTER_KEY in adata.obs.columns, f"Requires obs['{CLUSTER_KEY}']"

    part1_celltypist(adata)
    part2_original(adata)
    part3_markers(adata)
    cluster_df = part4_decide(adata)

    log(f"Writing {IN_PATH} (overwrite)")
    adata.write_h5ad(IN_PATH, compression="gzip")
    size_gb = IN_PATH.stat().st_size / 1e9
    log(f"Written ({size_gb:.2f} GB)")

    part5_report(adata, cluster_df)

    elapsed = (datetime.now() - t0).total_seconds() / 60
    log(f"\nElapsed: {elapsed:.1f} min")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
