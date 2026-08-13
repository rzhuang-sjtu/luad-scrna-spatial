#!/usr/bin/env python3
"""
Step 2: Merge 7 LUAD scRNA-seq datasets into one h5ad

Usage:
  python 02_merge_datasets.py --dry-run   # check field mapping only; do not load expression
  python 02_merge_datasets.py             # full merge
"""
import argparse
import gc
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

DATA_DIR = Path("${WORK_ROOT}/数据清洗")
OUT_DIR = Path("${PROJECT_ROOT}/data/processed")
RESULTS_DIR = Path("${PROJECT_ROOT}/results")

DATASETS = {
    "GSE123902": "GSE123902_clean.h5ad",
    "GSE131907": "GSE131907_clean.h5ad",
    "GSE143423": "GSE143423_LUAD_clean.h5ad",
    "GSE148071": "GSE148071_LUAD_clean.h5ad",
    "GSE164789": "GSE164789_LUAD_clean.h5ad",
    "GSE189357": "GSE189357_LUAD_clean.h5ad",
    "GSE253013": "GSE253013_clean.h5ad",
}

STANDARD_CATEGORIES = {
    "Normal_Lung", "Adjacent_Normal", "Normal_LN", "Precancerous",
    "Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
    "Distant_Metastasis", "Pleural_Effusion", "Other"
}

CORE_OBS_COLS = [
    "dataset", "patient_id", "sample_id",
    "tissue_type_original", "tissue_type", "tissue_stage",
    "celltype_original",
    "n_genes_by_counts", "total_counts",
    "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",
    "doublet_score", "predicted_doublet",
    "stage", "chemotherapy",
]
NUMERIC_COLS = {
    "n_genes_by_counts", "total_counts",
    "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",
    "doublet_score",
}

def map_tissue_type(dataset, obs):
    """Return an obs copy with tissue_type_new / tissue_stage columns"""
    obs = obs.copy()
    tt = obs["tissue_type"].astype(str) if "tissue_type" in obs.columns else pd.Series([""] * len(obs), index=obs.index)

    if dataset == "GSE123902":
        # Paired _N samples with a T sample → Adjacent_Normal; isolated _N → Normal_Lung
        paired_n = {"LX675_N", "LX682_N", "LX684_N"}
        def _row(sid, t):
            if t == "Tumor": return ("Primary_Tumor", "")
            if t == "Normal":
                return ("Adjacent_Normal" if sid in paired_n else "Normal_Lung", "")
            if t == "Brain_met": return ("Brain_Metastasis", "")
            if t in ("Adrenal_met", "Bone_met"):
                return ("Distant_Metastasis", t.replace("_met", ""))
            return ("Other", "")
        res = [_row(str(s), str(x)) for s, x in zip(obs["sample_id"], tt)]
        obs["tissue_type_new"] = [r[0] for r in res]
        obs["tissue_stage"] = [r[1] for r in res]

    elif dataset == "GSE131907":
        m = {
            "Primary_Tumor": "Primary_Tumor",
            "Normal_Lung": "Normal_Lung",
            "LN_Normal": "Normal_LN",
            "Brain_Met": "Brain_Metastasis",
            "Pleural_Effusion": "Pleural_Effusion",
            "LN_Met": "LN_Metastasis",
            "Advanced_Tumor": "Primary_Tumor",
        }
        obs["tissue_type_new"] = tt.map(m).fillna("Other")
        obs["tissue_stage"] = np.where(tt == "Advanced_Tumor", "Advanced", "")

    elif dataset == "GSE143423":
        obs["tissue_type_new"] = "Brain_Metastasis"
        obs["tissue_stage"] = ""

    elif dataset == "GSE148071":
        obs["tissue_type_new"] = "Primary_Tumor"
        obs["tissue_stage"] = "Advanced"   # stage III/IV advanced LUAD

    elif dataset == "GSE164789":
        def _cls(sid):
            sid = str(sid)
            if sid.endswith(("-T", "-T1", "-T2", "-KT")):
                stage = "K" if sid.endswith("-KT") else ""
                return ("Primary_Tumor", stage)
            if sid.endswith(("-A", "-A1", "-A2", "-KA")):
                stage = "K" if sid.endswith("-KA") else ""
                return ("Precancerous", stage)
            return ("Other", "")
        res = obs["sample_id"].astype(str).apply(_cls)
        obs["tissue_type_new"] = res.apply(lambda x: x[0])
        obs["tissue_stage"] = res.apply(lambda x: x[1])

    elif dataset == "GSE189357":
        pat = re.compile(r"TD(\d+)")
        def _cls(sid):
            m = pat.search(str(sid))
            if not m: return ("Other", "")
            n = int(m.group(1))
            if 1 <= n <= 3: return ("Precancerous", "AIS")
            if 4 <= n <= 6: return ("Precancerous", "MIA")
            if 7 <= n <= 9: return ("Primary_Tumor", "IAC")
            return ("Other", "")
        res = obs["sample_id"].astype(str).apply(_cls)
        obs["tissue_type_new"] = res.apply(lambda x: x[0])
        obs["tissue_stage"] = res.apply(lambda x: x[1])

    elif dataset == "GSE253013":
        m = {"Tumor": "Primary_Tumor", "Normal": "Adjacent_Normal"}
        obs["tissue_type_new"] = tt.map(m).fillna("Other")
        obs["tissue_stage"] = ""

    return obs


def _is_integer_matrix(x, n=200):
    if sparse.issparse(x):
        s = x[:min(n, x.shape[0])].toarray()
    else:
        s = np.asarray(x[:min(n, x.shape[0])])
    if s.size == 0:
        return False
    return np.allclose(s, s.astype(int))


def _build_obs(dataset, raw_obs):
    obs = map_tissue_type(dataset, raw_obs)
    obs["dataset"] = dataset
    obs["tissue_type_original"] = raw_obs["tissue_type"].astype(str) if "tissue_type" in raw_obs.columns else ""
    obs["tissue_type"] = obs["tissue_type_new"].astype(str)
    obs = obs.drop(columns=["tissue_type_new"])

    # Original celltype annotation
    if "Cell_type" in raw_obs.columns:
        obs["celltype_original"] = raw_obs["Cell_type"].astype(str)
    elif "cell_type" in raw_obs.columns:
        obs["celltype_original"] = raw_obs["cell_type"].astype(str)
    else:
        obs["celltype_original"] = ""

    # Fill core columns
    for col in CORE_OBS_COLS:
        if col not in obs.columns:
            obs[col] = np.nan if col in NUMERIC_COLS else ""
    # Force-cast numeric columns
    for col in NUMERIC_COLS:
        obs[col] = pd.to_numeric(obs[col], errors="coerce")
    if "predicted_doublet" in obs.columns:
        obs["predicted_doublet"] = obs["predicted_doublet"].astype(str)

    return obs[CORE_OBS_COLS].copy()


def normalize_one(dataset, fname, dry_run=False):
    path = DATA_DIR / fname
    print(f"\n=== {dataset} ({fname}) ===", flush=True)
    if not path.exists():
        raise FileNotFoundError(path)

    adata = sc.read_h5ad(path)
    print(f"Raw: {adata.shape[0]:,} × {adata.shape[1]:,}", flush=True)

    # 1) Locate raw counts source
    if _is_integer_matrix(adata.X):
        print(f"X is integer counts; use directly")
        counts = adata.X
    elif "counts" in adata.layers and _is_integer_matrix(adata.layers["counts"]):
        print(f"X is not counts; use layers['counts']")
        counts = adata.layers["counts"]
    else:
        raise ValueError(f"{dataset}: integer raw counts not found")

    # 2) Deduplicate var
    if adata.var_names.duplicated().any():
        dup_n = adata.var_names.duplicated().sum()
        print(f"Drop duplicate var_names ({dup_n})")
        keep = ~adata.var_names.duplicated()
        counts = counts[:, keep]
        var_names = adata.var_names[keep]
    else:
        var_names = adata.var_names

    # 3) Standardize obs
    new_obs = _build_obs(dataset, adata.obs)
    new_obs.index = pd.Index(
        [f"{dataset}_{i}" for i in adata.obs.index.astype(str)],
        name="cell_id"
    )

    # 4) Print distributions
    print(f"tissue_type distribution: {dict(new_obs['tissue_type'].value_counts())}")
    stg = new_obs['tissue_stage'].replace("", "(empty)").value_counts().to_dict()
    print(f"tissue_stage distribution: {stg}")

    if dry_run:
        del adata; gc.collect()
        return {
            "dataset": dataset,
            "n_cells": len(new_obs),
            "var_names": list(var_names),
            "obs": new_obs.head(3),
            "tt_counts": new_obs["tissue_type"].value_counts().to_dict(),
        }

    # 5) Build new sparse AnnData
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    counts = counts.astype(np.float32)

    new_ad = ad.AnnData(
        X=counts,
        obs=new_obs,
        var=pd.DataFrame(index=var_names.copy()),
    )
    print(f"After standardization: {new_ad.shape[0]:,} × {new_ad.shape[1]:,}  (dtype={new_ad.X.dtype})")

    del adata; gc.collect()
    return new_ad


def _make_report(adata):
    L = ["# Step 2 merge report", ""]
    L.append(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- Total cells: {adata.shape[0]:,}")
    L.append(f"- Shared genes (var intersection): {adata.shape[1]:,}")
    L.append(f"- X dtype: {adata.X.dtype}, sparse: {sparse.issparse(adata.X)}")
    L.append("")
    L.append("## Cells per dataset")
    for ds, n in adata.obs["dataset"].value_counts().items():
        L.append(f"- {ds}: {n:,}")
    L.append("")
    L.append("## tissue_type × dataset crosstab")
    ct = pd.crosstab(adata.obs["dataset"], adata.obs["tissue_type"])
    L.append("```")
    L.append(ct.to_string())
    L.append("```")
    L.append("")
    L.append("## tissue_stage distribution")
    stg = adata.obs["tissue_stage"].replace("", "(empty)").value_counts()
    for s, n in stg.items():
        L.append(f"- {s}: {n:,}")
    L.append("")
    L.append("## Non-standard tissue_type (should be Other only)")
    bad = adata.obs[~adata.obs["tissue_type"].isin(STANDARD_CATEGORIES)]
    if len(bad) == 0:
        L.append("- All standardized")
    else:
        L.append(f"-  {len(bad):,} cells not in standard categories")
        L.append("```")
        L.append(bad["tissue_type_original"].value_counts().to_string())
        L.append("```")
    L.append("")
    L.append("## Unique patient / sample counts")
    L.append(f"- n_patients: {adata.obs['patient_id'].nunique()}")
    L.append(f"- n_samples: {adata.obs['sample_id'].nunique()}")
    return "\n".join(L)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print("=" * 60); print("DRY RUN mode"); print("=" * 60)
        summaries = [normalize_one(ds, f, dry_run=True) for ds, f in DATASETS.items()]
        common = set(summaries[0]["var_names"])
        for s in summaries[1:]:
            common &= set(s["var_names"])
        print(f"\n>>> Shared genes (intersection of 7 datasets): {len(common):,}")
        print(f">>> Total cells: {sum(s['n_cells'] for s in summaries):,}")
        return

    print("=" * 60)
    print(f"Full merge started: {datetime.now()}")
    print("=" * 60)

    adatas = []
    for ds, f in DATASETS.items():
        adatas.append(normalize_one(ds, f, dry_run=False))

    print(f"\nMerging (ad.concat, join='inner') ...", flush=True)
    merged = ad.concat(
        adatas, axis=0, join="inner", merge="same",
        label=None, index_unique=None,
    )
    # Free intermediate adata
    del adatas; gc.collect()

    # Normalize obs dtypes (categorical columns to string after merge)
    for col in ["dataset", "patient_id", "sample_id",
                "tissue_type", "tissue_type_original", "tissue_stage",
                "celltype_original", "stage", "chemotherapy",
                "predicted_doublet"]:
        if col in merged.obs.columns:
            merged.obs[col] = merged.obs[col].astype(str)

    print(f"\nMerged result: {merged.shape[0]:,} × {merged.shape[1]:,}")
    print(f"X: dtype={merged.X.dtype}, sparse={sparse.issparse(merged.X)}")

    # Write h5ad first (preserve merge result)
    out_path = OUT_DIR / "luad_merged_raw.h5ad"
    print(f"\nWriting {out_path} (gzip) ...", flush=True)
    merged.write_h5ad(out_path, compression="gzip")
    size_gb = out_path.stat().st_size / 1e9
    print(f"h5ad done: {out_path} ({size_gb:.2f} GB)")

    # Then write report (failure does not affect h5ad)
    try:
        report = _make_report(merged)
        report_path = RESULTS_DIR / "step2_merge_summary.md"
        report_path.write_text(report, encoding="utf-8")
        print(f"Report: {report_path}")
    except Exception as e:
        print(f"Report generation failed but h5ad saved: {e}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
