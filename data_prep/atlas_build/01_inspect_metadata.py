#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
01_inspect_metadata.py

Inspect obs metadata of 7 cleaned h5ad files to prepare unified meta-programme integration.
Read-only with backed='r'; does not load X into memory.
"""

from __future__ import annotations

import gc
import re
import traceback
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd


INPUT_FILES: "OrderedDict[str, str]" = OrderedDict([
    ("GSE123902",        "${WORK_ROOT}/数据清洗/GSE123902_clean.h5ad"),
    ("GSE131907",        "${WORK_ROOT}/数据清洗/GSE131907_clean.h5ad"),
    ("GSE143423_LUAD",   "${WORK_ROOT}/数据清洗/GSE143423_LUAD_clean.h5ad"),
    ("GSE148071_LUAD",   "${WORK_ROOT}/数据清洗/GSE148071_LUAD_clean.h5ad"),
    ("GSE164789_LUAD",   "${WORK_ROOT}/数据清洗/GSE164789_LUAD_clean.h5ad"),
    ("GSE189357_LUAD",   "${WORK_ROOT}/数据清洗/GSE189357_LUAD_clean.h5ad"),
    ("GSE253013",        "${WORK_ROOT}/数据清洗/GSE253013_clean.h5ad"),
])

OUT_DIR = Path("${PROJECT_ROOT}/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

REPORT_PATH  = OUT_DIR / "metadata_inspection.md"
MAPPING_PATH = OUT_DIR / "tissue_type_mapping.csv"
MATRIX_PATH  = OUT_DIR / "cell_count_summary.csv"

# Candidate field keywords (case-insensitive match)
TISSUE_KEYWORDS    = ["tissue_type", "sample_type", "origin", "source", "sample_origin"]
PATIENT_KEYWORDS   = ["patient_id", "sample_id", "donor", "samples", "sample", "patientid"]
CELLTYPE_KEYWORDS  = ["celltype", "cell_type"]
MALIGNANT_KEYWORDS = ["malignant", "copykat", "cnv", "infercnv", "tumor_cell", "tumour_cell"]
QC_FIELDS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

STANDARD_CATEGORIES = [
    "Normal_Lung", "Adjacent_Normal", "Normal_LN", "Precancerous",
    "Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
    "Distant_Metastasis", "Pleural_Effusion", "Other",
]


def find_columns(cols: Iterable[str], keywords: Iterable[str]) -> list[str]:
    cols_list = list(cols)
    low = {c.lower(): c for c in cols_list}
    hits = []
    for kw in keywords:
        for lc, c in low.items():
            if kw in lc and c not in hits:
                hits.append(c)
    return hits


def fmt_counts(vc: pd.Series, limit: int = 60) -> str:
    items = [f"  - `{k}`: {int(v)}" for k, v in vc.items()][:limit]
    tail = "" if len(vc) <= limit else f"\n  - ... ({len(vc)} values total; showing first {limit})"
    return "\n".join(items) + tail


def describe_per_group(ser: pd.Series, group: pd.Series) -> pd.Series:
    counts = group.value_counts()
    q = counts.describe(percentiles=[0.25, 0.5, 0.75])
    return q


def qc_summary(obs: pd.DataFrame) -> dict:
    res = {}
    for f in QC_FIELDS:
        if f in obs.columns:
            s = pd.to_numeric(obs[f], errors="coerce").dropna()
            if len(s) > 0:
                res[f] = {
                    "min": float(s.min()),
                    "median": float(s.median()),
                    "max": float(s.max()),
                    "mean": float(s.mean()),
                }
            else:
                res[f] = None
        else:
            res[f] = None
    return res


def prefix_group_131907(s: str) -> str:
    prefixes = ["LUNG_T", "LUNG_N", "LN_", "mBrain", "EBUS_",
                "BRONCHO_", "EFFUSION_", "NS_"]
    for p in prefixes:
        if s.startswith(p):
            return p
    return "Other"


def suffix_group_164789(s: str) -> str:
    # Match longer suffixes first
    for suf in ["-KA", "-KT", "-A1", "-A2", "-T1", "-T2", "-A", "-T"]:
        if s.endswith(suf):
            return suf
    return "Other"


def td_group_189357(s: str) -> str:
    m = re.match(r"TD(\d+)", s)
    if not m:
        return "Other"
    n = int(m.group(1))
    if 1 <= n <= 3:
        return "AIS (TD1-3)"
    if 4 <= n <= 6:
        return "MIA (TD4-6)"
    if 7 <= n <= 9:
        return "IAC (TD7-9)"
    return "Other"


def ant_adeno_group_253013(s: str) -> str:
    s_up = s.upper()
    if "_NAT" in s_up:
        return "NAT (Adjacent_Normal)"
    if "_ADENO" in s_up:
        return "ADENO (Primary_Tumor)"
    if "_NSCLC" in s_up:
        return "NSCLC (Primary_Tumor)"
    if "_ANT" in s_up or s_up.startswith("ANT"):
        return "ANT (Adjacent_Normal)"
    return "Other"


def build_mapping_table(per_dataset_tissue: dict) -> pd.DataFrame:
    """
    Build a mapping table from each dataset's raw tissue_type values.
    Returns four columns: dataset / original_value / standard_category / notes
    """
    rows = []

    # GSE123902: Tumor / Normal / Brain_met / Adrenal_met / Bone_met
    for v in per_dataset_tissue.get("GSE123902", {}):
        if v == "Normal":
            rows.append(("GSE123902", v, "Normal_Lung",
                         "Normal lung from patients sampled at the primary site → Normal_Lung; if from metastatic patients, use Adjacent_Normal (check sample_id)"))
        elif v == "Tumor":
            rows.append(("GSE123902", v, "Primary_Tumor", "Primary lung adenocarcinoma"))
        elif v == "Brain_met":
            rows.append(("GSE123902", v, "Brain_Metastasis", "Brain metastasis"))
        elif v == "Bone_met":
            rows.append(("GSE123902", v, "Distant_Metastasis", "Bone metastasis"))
        elif v == "Adrenal_met":
            rows.append(("GSE123902", v, "Distant_Metastasis", "Adrenal metastasis"))
        else:
            rows.append(("GSE123902", v, "Other", "Unmapped; confirm"))

    # GSE131907: already has standardized tissue_type
    for v in per_dataset_tissue.get("GSE131907", {}):
        if v == "Primary_Tumor":
            rows.append(("GSE131907", v, "Primary_Tumor", ""))
        elif v == "Normal_Lung":
            rows.append(("GSE131907", v, "Normal_Lung", ""))
        elif v == "LN_Normal":
            rows.append(("GSE131907", v, "Normal_LN", ""))
        elif v == "Brain_Met":
            rows.append(("GSE131907", v, "Brain_Metastasis", ""))
        elif v == "Pleural_Effusion":
            rows.append(("GSE131907", v, "Pleural_Effusion", ""))
        elif v == "LN_Met":
            rows.append(("GSE131907", v, "LN_Metastasis", ""))
        elif v == "Advanced_Tumor":
            rows.append(("GSE131907", v, "Primary_Tumor",
                         "tL/B = tLung/Bronchus advanced tumour; clinically primary/locally advanced; may keep a separate Advanced label for sub-analyses"))
        else:
            rows.append(("GSE131907", v, "Other", "Unmapped"))

    # GSE143423: Brain_met only
    for v in per_dataset_tissue.get("GSE143423_LUAD", {}):
        if v == "Brain_met":
            rows.append(("GSE143423_LUAD", v, "Brain_Metastasis", "3 LUAD brain-metastasis samples"))
        else:
            rows.append(("GSE143423_LUAD", v, "Other", "Unexpected category; check"))

    # GSE148071: Tumor (advanced NSCLC / LUAD-filtered)
    for v in per_dataset_tissue.get("GSE148071_LUAD", {}):
        if v == "Tumor":
            rows.append(("GSE148071_LUAD", v, "Primary_Tumor",
                         "Advanced LUAD, stage III/IV; pleural involvement possible at sampling, not further subdivided in obs"))
        else:
            rows.append(("GSE148071_LUAD", v, "Other", "Unmapped"))

    # GSE164789: Normal / Tumor; distinguish by sample_id suffix
    # Task note: "-A suffix is precancerous"; map Normal→Precancerous as a special case
    for v in per_dataset_tissue.get("GSE164789_LUAD", {}):
        if v == "Normal":
            rows.append(("GSE164789_LUAD", v, "Precancerous",
                         "sample_id ending -A/-A1/-A2/-KA treated as precancerous (AAH/AIS) per task note, not true normal lung; if adjacent normal, change to Adjacent_Normal"))
        elif v == "Tumor":
            rows.append(("GSE164789_LUAD", v, "Primary_Tumor",
                         "Primary tumour with sample_id ending -T/-T1/-T2/-KT"))
        else:
            rows.append(("GSE164789_LUAD", v, "Other", "Unmapped"))

    # GSE189357: Tumor, by TD1-3 AIS / TD4-6 MIA / TD7-9 IAC
    for v in per_dataset_tissue.get("GSE189357_LUAD", {}):
        if v == "Tumor":
            rows.append(("GSE189357_LUAD", v, "Precancerous / Primary_Tumor",
                         "TD1-3 (AIS) and TD4-6 (MIA) → Precancerous; TD7-9 (IAC) → Primary_Tumor; split by sample_id when building the matrix"))
        else:
            rows.append(("GSE189357_LUAD", v, "Other", "Unmapped"))

    # GSE253013: Tumor / Normal
    for v in per_dataset_tissue.get("GSE253013", {}):
        if v == "Normal":
            rows.append(("GSE253013", v, "Adjacent_Normal",
                         "sample_id containing _NAT is ANT / adjacent normal tissue"))
        elif v == "Tumor":
            rows.append(("GSE253013", v, "Primary_Tumor",
                         "_NSCLC / _ADENO / US-* primary sites; Site is Primary for all; no metastases"))
        else:
            rows.append(("GSE253013", v, "Other", "Unmapped"))

    df = pd.DataFrame(rows, columns=["dataset", "original_value", "standard_category", "notes"])
    return df


def compute_standard_counts(dataset: str, obs: pd.DataFrame) -> dict:
    """
    Map each obs row to a standard category from tissue_type + sample_id rules; return {standard_category: n_cells}.
    """
    counts = {c: 0 for c in STANDARD_CATEGORIES}

    tissue = obs["tissue_type"].astype(str) if "tissue_type" in obs.columns else pd.Series([""] * len(obs), index=obs.index)
    sid    = obs["sample_id"].astype(str) if "sample_id" in obs.columns else pd.Series([""] * len(obs), index=obs.index)

    if dataset == "GSE123902":
        mp = {"Normal": "Normal_Lung", "Tumor": "Primary_Tumor",
              "Brain_met": "Brain_Metastasis", "Bone_met": "Distant_Metastasis",
              "Adrenal_met": "Distant_Metastasis"}
        mapped = tissue.map(mp).fillna("Other")

    elif dataset == "GSE131907":
        mp = {"Primary_Tumor": "Primary_Tumor", "Normal_Lung": "Normal_Lung",
              "LN_Normal": "Normal_LN", "Brain_Met": "Brain_Metastasis",
              "Pleural_Effusion": "Pleural_Effusion", "LN_Met": "LN_Metastasis",
              "Advanced_Tumor": "Primary_Tumor"}
        mapped = tissue.map(mp).fillna("Other")

    elif dataset == "GSE143423_LUAD":
        mp = {"Brain_met": "Brain_Metastasis"}
        mapped = tissue.map(mp).fillna("Other")

    elif dataset == "GSE148071_LUAD":
        mp = {"Tumor": "Primary_Tumor"}
        mapped = tissue.map(mp).fillna("Other")

    elif dataset == "GSE164789_LUAD":
        mapped = pd.Series(["Other"] * len(obs), index=obs.index)
        mapped[tissue == "Tumor"] = "Primary_Tumor"
        mapped[tissue == "Normal"] = "Precancerous"

    elif dataset == "GSE189357_LUAD":
        mapped = pd.Series(["Other"] * len(obs), index=obs.index)
        td_grp = sid.apply(td_group_189357)
        mapped[td_grp.isin(["AIS (TD1-3)", "MIA (TD4-6)"])] = "Precancerous"
        mapped[td_grp == "IAC (TD7-9)"] = "Primary_Tumor"

    elif dataset == "GSE253013":
        mp = {"Normal": "Adjacent_Normal", "Tumor": "Primary_Tumor"}
        mapped = tissue.map(mp).fillna("Other")

    else:
        mapped = pd.Series(["Other"] * len(obs), index=obs.index)

    vc = mapped.value_counts()
    for k, v in vc.items():
        counts[k] = counts.get(k, 0) + int(v)
    return counts


def analyze_one(dataset: str, path: str) -> dict:
    report_lines: list[str] = []
    info = {"dataset": dataset, "path": path, "ok": False,
            "report": "", "warnings": [], "std_counts": {},
            "tissue_values": [], "n_cells": None, "n_genes": None,
            "n_patients": None, "has_celltype": False, "has_malignant": False}

    def log(s=""):
        report_lines.append(s)

    log(f"### {dataset}")
    log(f"File: `{path}`")

    if not Path(path).exists():
        info["warnings"].append(f"File not found: {path}")
        log(f"\n File not found; skip.")
        info["report"] = "\n".join(report_lines)
        return info

    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as e:
        info["warnings"].append(f"Read failed: {e}")
        log(f"\n Read failed: {e}")
        log("```")
        log(traceback.format_exc())
        log("```")
        info["report"] = "\n".join(report_lines)
        return info

    try:
        obs = adata.obs.copy()
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        info["n_cells"] = int(n_cells)
        info["n_genes"] = int(n_genes)

        log(f"\n- Shape: {n_cells} cells × {n_genes} genes")

        # ---- (a) obs column names ----
        log(f"\n**a. obs columns ({len(obs.columns)}):**")
        log("```")
        log(", ".join(map(str, obs.columns)))
        log("```")

        # ---- (b) tissue / origin fields ----
        tissue_cols = find_columns(obs.columns, TISSUE_KEYWORDS)
        log(f"\n**b. Tissue/origin fields: {tissue_cols if tissue_cols else 'none detected'}**")
        for c in tissue_cols:
            vc = obs[c].astype(str).value_counts()
            log(f"- `{c}` unique values ({len(vc)}):")
            log(fmt_counts(vc))
            if c == "tissue_type":
                info["tissue_values"] = vc.index.tolist()

        # ---- (c) patient / sample fields ----
        pat_cols = find_columns(obs.columns, PATIENT_KEYWORDS)
        log(f"\n**c. Patient/sample fields: {pat_cols if pat_cols else 'none detected'}**")

        n_patients = None
        if "patient_id" in obs.columns:
            n_patients = int(obs["patient_id"].nunique())
            info["n_patients"] = n_patients
            log(f"- Unique patient_id: {n_patients}")
        elif "PatientID" in obs.columns:
            n_patients = int(obs["PatientID"].nunique())
            info["n_patients"] = n_patients
            log(f"- Unique PatientID: {n_patients}")

        if "sample_id" in obs.columns:
            n_samples = int(obs["sample_id"].nunique())
            log(f"- Unique sample_id: {n_samples}")
            counts_per_sample = obs["sample_id"].value_counts()
            q = counts_per_sample.describe(percentiles=[0.25, 0.5, 0.75])
            log(f"- Cells per sample: min={int(q['min'])},"
                f"25%={int(q['25%'])}, median={int(q['50%'])}, "
                f"75%={int(q['75%'])}, max={int(q['max'])}")

        if "patient_id" in obs.columns:
            counts_per_pat = obs["patient_id"].value_counts()
            q = counts_per_pat.describe(percentiles=[0.25, 0.5, 0.75])
            log(f"- Cells per patient: min={int(q['min'])},"
                f"25%={int(q['25%'])}, median={int(q['50%'])}, "
                f"75%={int(q['75%'])}, max={int(q['max'])}")

        # ---- (d) celltype fields ----
        ct_cols = find_columns(obs.columns, CELLTYPE_KEYWORDS)
        if ct_cols:
            info["has_celltype"] = True
        log(f"\n**d. Cell-type fields: {ct_cols if ct_cols else 'none'}**")
        for c in ct_cols:
            vc = obs[c].astype(str).value_counts()
            log(f"- `{c}` unique values ({len(vc)}):")
            log(fmt_counts(vc))

        # ---- (e) malignant / CNV fields ----
        mal_cols = find_columns(obs.columns, MALIGNANT_KEYWORDS)
        info["has_malignant"] = bool(mal_cols)
        log(f"\n**e. Tumour-cell label fields (malignant/copykat/CNV/inferCNV, etc.):**"
            f"{mal_cols if mal_cols else ' none detected'}")

        # ---- (f) QC field distributions ----
        log(f"\n**f. QC field distributions:**")
        qc = qc_summary(obs)
        for k, v in qc.items():
            if v is None:
                log(f"- `{k}`: absent")
            else:
                log(f"- `{k}`: min={v['min']:.2f}, median={v['median']:.2f}, "
                    f"max={v['max']:.2f}, mean={v['mean']:.2f}")

        # ---- Special handling ----
        log(f"\n**Special group stats:**")
        if dataset == "GSE131907" and "sample_id" in obs.columns:
            grp = obs["sample_id"].astype(str).apply(prefix_group_131907)
            grp_counts = grp.value_counts()
            for g, n in grp_counts.items():
                tag = "NS_ sample origin needs further verification" if g == "NS_" else ""
                log(f"- Prefix `{g}`: {int(n)} cells{tag}")
            if "NS_" in grp_counts.index:
                info["warnings"].append("GSE131907 has NS_ prefix samples; origin unclear")

        elif dataset == "GSE164789_LUAD" and "sample_id" in obs.columns:
            grp = obs["sample_id"].astype(str).apply(suffix_group_164789)
            grp_counts = grp.value_counts()
            notes_map = {
                "-A": "Precancerous (per task note)",
                "-A1": "Precancerous (per task note)",
                "-A2": "Precancerous (per task note)",
                "-T": "Primary tumour",
                "-T1": "Primary tumour",
                "-T2": "Primary tumour",
                "-KA": "Kidney-transplant-related A suffix; meaning to verify",
                "-KT": "Kidney-transplant-related T suffix; meaning to verify",
            }
            for g, n in grp_counts.items():
                log(f"- Suffix `{g}`: {int(n)} cells — {notes_map.get(g, '')}")

        elif dataset == "GSE189357_LUAD" and "sample_id" in obs.columns:
            grp = obs["sample_id"].astype(str).apply(td_group_189357)
            grp_counts = grp.value_counts()
            for g, n in grp_counts.items():
                log(f"- {g}: {int(n)} cells")

        elif dataset == "GSE253013" and "sample_id" in obs.columns:
            grp = obs["sample_id"].astype(str).apply(ant_adeno_group_253013)
            grp_counts = grp.value_counts()
            for g, n in grp_counts.items():
                tag = "— adjacent normal tissue" if "NAT" in g or "ANT" in g else ""
                log(f"- `{g}`: {int(n)} cells{tag}")

            if "Site" in obs.columns:
                site_vc = obs["Site"].astype(str).value_counts()
                log(f"- Site values: {site_vc.to_dict()}")
            if "Histology" in obs.columns:
                hist_vc = obs["Histology"].astype(str).value_counts()
                log(f"- Histology values: {hist_vc.to_dict()}")

        elif dataset == "GSE123902" and "sample_id" in obs.columns:
            xt = pd.crosstab(obs["sample_id"].astype(str), obs["tissue_type"].astype(str))
            log("- sample_id × tissue_type crosstab:")
            log("```")
            log(xt.to_string())
            log("```")
            if "Brain_met" in xt.columns:
                bm = xt.index[xt["Brain_met"] > 0].tolist()
                log(f"- Brain-met samples: {bm}")
            if "Bone_met" in xt.columns:
                bn = xt.index[xt["Bone_met"] > 0].tolist()
                log(f"- Bone-met samples: {bn}")
            if "Adrenal_met" in xt.columns:
                ad_ = xt.index[xt["Adrenal_met"] > 0].tolist()
                log(f"- Adrenal-met samples: {ad_}")

        elif dataset == "GSE143423_LUAD":
            n_s = obs["sample_id"].nunique() if "sample_id" in obs.columns else 0
            log(f"- Samples: {n_s}")
            if n_s != 3:
                info["warnings"].append(
                    f"GSE143423 expected 3 LUAD brain-met samples, got {n_s}; check for residual TNBC"
                )
                log(f"Expected 3 samples, got {n_s}; check for residual TNBC")
            if "tissue_type" in obs.columns:
                unique_tt = set(obs["tissue_type"].astype(str).unique())
                if unique_tt - {"Brain_met"}:
                    info["warnings"].append(
                        f"GSE143423 tissue_type has non-Brain_met values: {unique_tt - {'Brain_met'}}")

        elif dataset == "GSE148071_LUAD":
            n_s = obs["sample_id"].nunique() if "sample_id" in obs.columns else 0
            log(f"- LUAD samples: {n_s} (original cohort 42 advanced NSCLC; LUAD only after cleaning)")
            if "cancer_type" in obs.columns:
                vc = obs["cancer_type"].astype(str).value_counts()
                log(f"- cancer_type：{vc.to_dict()}")
                if set(vc.index) - {"LUAD"}:
                    info["warnings"].append(
                        f"GSE148071 cancer_type has non-LUAD: {set(vc.index) - {'LUAD'}}"
                    )
            if n_s != 42:
                log(f"Currently {n_s} LUAD samples (LUSC removed from original 42 NSCLC)")

        # ---- Standard category counts ----
        std_counts = compute_standard_counts(dataset, obs)
        info["std_counts"] = std_counts
        log(f"\n**Cells per standard category:**")
        for k in STANDARD_CATEGORIES:
            if std_counts.get(k, 0) > 0:
                log(f"- {k}: {std_counts[k]}")

        info["ok"] = True

    except Exception as e:
        info["warnings"].append(f"Analysis failed: {e}")
        log(f"\n Analysis error: {e}")
        log("```")
        log(traceback.format_exc())
        log("```")

    finally:
        try:
            adata.file.close()
        except Exception:
            pass
        del adata
        gc.collect()

    info["report"] = "\n".join(report_lines)
    return info


def main():
    all_info = OrderedDict()
    per_dataset_tissue = OrderedDict()

    for ds, path in INPUT_FILES.items():
        print(f"\n>>> Analysing {ds} ...")
        info = analyze_one(ds, path)
        all_info[ds] = info
        per_dataset_tissue[ds] = info.get("tissue_values", [])

    # ---- Summary table (section 1) ----
    overview_rows = []
    for ds, info in all_info.items():
        overview_rows.append({
            "dataset": ds,
            "n_cells": info["n_cells"],
            "n_genes": info["n_genes"],
            "n_patients": info["n_patients"],
            "has_celltype": "yes" if info["has_celltype"] else "否",
            "has_malignant": "yes" if info["has_malignant"] else "否",
        })
    overview_df = pd.DataFrame(overview_rows)

    # ---- Mapping table (section 3) ----
    mapping_df = build_mapping_table(per_dataset_tissue)
    mapping_df.to_csv(MAPPING_PATH, index=False, encoding="utf-8-sig")

    # ---- Cell-count matrix ----
    matrix_rows = []
    for ds, info in all_info.items():
        row = {"dataset": ds}
        for c in STANDARD_CATEGORIES:
            row[c] = info["std_counts"].get(c, 0)
        matrix_rows.append(row)
    matrix_df = pd.DataFrame(matrix_rows)
    matrix_df.to_csv(MATRIX_PATH, index=False, encoding="utf-8-sig")

    # ---- Write markdown report ----
    lines: list[str] = []
    lines.append("# Metadata audit of 7 LUAD / NSCLC single-cell datasets")
    lines.append("")
    lines.append("> Auto-generated by `data_prep/atlas_build/01_inspect_metadata.py`.")
    lines.append(f"> Outputs: `{REPORT_PATH}`, `{MAPPING_PATH}`, `{MATRIX_PATH}`")
    lines.append("")

    # Section 1
    lines.append("## Section 1: Overview of the 7 datasets")
    lines.append("")
    lines.append("| dataset | n_cells | n_genes | n_patients | has celltype | has malignant |")
    lines.append("|---|---:|---:|---:|:---:|:---:|")
    for r in overview_rows:
        lines.append(f"| {r['dataset']} | {r['n_cells']} | {r['n_genes']} | "
                     f"{r['n_patients']} | {r['has_celltype']} | {r['has_malignant']} |")
    lines.append("")

    # Section 2: per-dataset detail
    lines.append("## Section 2: Per-dataset detail")
    lines.append("")
    for ds, info in all_info.items():
        lines.append(info["report"])
        lines.append("")
        lines.append("---")
        lines.append("")

    # Section 3: unified mapping table
    lines.append("## Section 3: Suggested unified mapping")
    lines.append("")
    lines.append("Map each dataset's raw `tissue_type` values to 10 standard categories:")
    lines.append(f"**{', '.join(STANDARD_CATEGORIES)}**")
    lines.append("")
    lines.append("| dataset | original_value | standard_category | notes |")
    lines.append("|---|---|---|---|")
    for _, r in mapping_df.iterrows():
        lines.append(f"| {r['dataset']} | {r['original_value']} | "
                     f"{r['standard_category']} | {r['notes']} |")
    lines.append("")
    lines.append("**Cell-count matrix after mapping (dataset × standard category):**")
    lines.append("")
    header = "| dataset | " + " | ".join(STANDARD_CATEGORIES) + " |"
    sep = "|---|" + "---:|" * len(STANDARD_CATEGORIES)
    lines.append(header)
    lines.append(sep)
    for r in matrix_rows:
        cells = [str(r[c]) for c in STANDARD_CATEGORIES]
        lines.append(f"| {r['dataset']} | " + " | ".join(cells) + " |")
    lines.append("")

    # Section 4: open issues
    warnings_all = []
    for ds, info in all_info.items():
        for w in info["warnings"]:
            warnings_all.append((ds, w))

    lines.append("## Section 4: Issues to resolve before integration")
    lines.append("")
    lines.append("### P0 Required")
    lines.append("")
    p0 = [
        "**Missing tumour-cell labels**: none of the 7 datasets has `malignant` / `copykat` / `inferCNV` in obs. Meta-programme analysis needs tumour vs non-tumour epithelium; run inferCNV or copyKAT first.",
        "**Inconsistent cell-type labels**: only GSE131907 and GSE253013 have obs-level `Cell_type` / `cell_type`; the other 5 do not. Re-annotate before integration (CellTypist immune_Low / Human_Lung_Atlas or scANVI transfer).",
        "**Confirm GSE164789 -A suffix meaning**: task note says `-A` is precancerous, but `tissue_type` is named `Normal`. Check Xing 2021 *Genome Med* whether AAH/AIS/MIA or adjacent normal.",
        "**Unknown GSE164789 `-KA` / `-KT` suffixes**: patient LRH has -KA / -KT; meaning unknown (repeat sampling or special site); cannot map without clarification.",
        "**GSE131907 `Advanced_Tumor` (tL/B) classification**: whether tLung/Bronchus advanced tumours join Primary_Tumor or stay separate must be fixed in the manuscript.",
    ]
    for i, x in enumerate(p0, 1):
        lines.append(f"{i}. {x}")
    lines.append("")

    lines.append("### P1 Recommended")
    lines.append("")
    p1 = [
        "**Origin of GSE123902 `Normal`**: check whether sample_id (`LX675_N` / `LX682_N` / `LX684_N` / `LX685`) patients also have Tumor samples → Adjacent_Normal if yes, else Normal_Lung.",
        "**GSE131907 `NS_` prefix samples**: if present, check the original paper sample table.",
        "**GSE253013 `Site` is all Primary**: align with Normal/Tumor in `tissue_type`; confirm Normal is NAT not distant healthy lung.",
        "**Epithelial fraction in GSE189357 AIS/MIA**: TD1-6 are precancerous; very low epithelial fractions may bias meta-programme learning if included as Precancerous.",
        "**QC threshold differences**: pct_counts_mt / n_genes_by_counts ranges differ; harmonize or use dataset as batch key in Harmony/scVI.",
    ]
    for i, x in enumerate(p1, 1):
        lines.append(f"{i}. {x}")
    lines.append("")

    lines.append("### P2 Optional")
    lines.append("")
    p2 = [
        "Whether GSE131907 `Advanced_Tumor (tL/B)` is a separate class or merged into Primary_Tumor.",
        "Whether Bone_met and Adrenal_met merge as Distant_Metastasis or stay separate.",
        "Whether stage III/IV LUAD in GSE148071 is listed as Advanced_Primary.",
    ]
    for i, x in enumerate(p2, 1):
        lines.append(f"{i}. {x}")
    lines.append("")

    if warnings_all:
        lines.append("### Runtime warnings")
        lines.append("")
        for ds, w in warnings_all:
            lines.append(f"- **{ds}**: {w}")
        lines.append("")

    # Section 5 action list
    lines.append("## Section 5: Suggested next actions")
    lines.append("")
    actions = [
        "**Unified cell-type annotation**: auto-annotate GSE123902 / GSE143423 / GSE148071 / GSE164789 / GSE189357 with CellTypist Human_Lung_Atlas (or scANVI); write `celltype_coarse`.",
        "**inferCNV or copyKAT for malignant labels**: use T/B/NK as reference on all epithelial cells; write `malignant` to obs.",
        "**Apply tissue_type mapping**: add `standard_category` per the section 3 table; split GSE164789 / GSE189357 by sample_id.",
        "**Verify GSE164789 -A/-KA meaning**: read original Methods or GEO GSE164789 `characteristics_ch1`.",
        "**Unify gene symbols**: intersect var.index across 7 datasets, then standardize aliases via HGNC/MyGene.",
        "**Build integrated AnnData**: concat with `dataset` as batch key; `sc.pp.normalize_total(1e4) + log1p`, then batch correction with scVI or Harmony; then meta-programme analysis (NMF / cNMF / Gavish 2023 MP pipeline).",
        "**Three-way contrast Precancerous / Adjacent_Normal / Primary_Tumor**: use GSE164789 + GSE189357 AAH→MIA→IAC samples to test meta-programme dynamics along progression.",
    ]
    for i, a in enumerate(actions, 1):
        lines.append(f"{i}. {a}")
    lines.append("")

    report_md = "\n".join(lines)
    # UTF-8 BOM prefix
    with open(REPORT_PATH, "w", encoding="utf-8-sig") as f:
        f.write(report_md)

    # ---- Print key sections to terminal ----
    print("\n" + "=" * 80)
    print("Report written:", REPORT_PATH)
    print("Mapping CSV :", MAPPING_PATH)
    print("Cell-count matrix CSV :", MATRIX_PATH)
    print("=" * 80)

    # Print sections 3 and 4
    sec3_start = report_md.find("## Section 3: Suggested unified mapping")
    sec5_start = report_md.find("## Section 5")
    print()
    print(report_md[sec3_start:sec5_start].rstrip())


if __name__ == "__main__":
    main()
