#!/usr/bin/env python3
"""
Step 4e: Merge per-patient CopyKAT predictions back into h5ad and assign final malignant labels.

Input:  ${PROJECT_ROOT}/data/processed/luad_integrated.h5ad
        ${PROJECT_ROOT}/data/copykat_input/{dp_key}/copykat_prediction.csv
Output:  ${PROJECT_ROOT}/data/processed/luad_copykat.h5ad
        ${PROJECT_ROOT}/results/step4_copykat_summary.md

Final obs['malignant'] rules (priority top to bottom):
  1. Normal_Lung / Adjacent_Normal / Normal_LN samples → Non-malignant
  2. celltype_coarse != Epithelial              → Non-malignant
     (Rationale: CopyKAT is valid for epithelial tumours; aneuploid calls in immune/stroma are false positives)
  3. Epithelial cells:
     - aneuploid                → Malignant
     - diploid                  → Non-malignant
     - not.defined / missing    → Uncertain
"""
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

IN_H5AD = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
OUT_H5AD = Path("${PROJECT_ROOT}/data/processed/luad_copykat.h5ad")
COPYKAT_DIR = Path("${PROJECT_ROOT}/data/copykat_input")
REPORT = Path("${PROJECT_ROOT}/results/step4_copykat_summary.md")

NORMAL_TT = {"Normal_Lung", "Adjacent_Normal", "Normal_LN"}


def normalize_pred(p: str) -> str:
    """CopyKAT output labels like 'c1:diploid:low.conf' / 'aneuploid' / 'not.defined'
    Map uniformly to aneuploid / diploid / not.defined"""
    p = str(p).strip().lower()
    if "aneuploid" in p:
        return "aneuploid"
    if "diploid" in p:
        return "diploid"
    return "not.defined"


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def load_predictions():
    """Walk copykat_input/{dp_key}/copykat_prediction.csv; return barcode→pred dict"""
    bc2pred = {}
    patient_status = []   # dp_key, status, n_cells_pred
    for pdir in sorted(COPYKAT_DIR.iterdir()):
        if not pdir.is_dir():
            continue
        dp = pdir.name
        pred_file = pdir / "copykat_prediction.csv"
        meta_file = pdir / "metadata.json"
        if not pred_file.exists() or pred_file.stat().st_size == 0:
            # Check for error
            err = pdir / "copykat_error.txt"
            reason = "missing"
            if err.exists():
                reason = "error: " + err.read_text()[:100].strip()
            patient_status.append({"dp_key": dp, "status": reason, "n_cells_pred": 0})
            continue

        try:
            df = pd.read_csv(pred_file)
        except Exception as e:
            patient_status.append({"dp_key": dp, "status": f"csv_error: {e}", "n_cells_pred": 0})
            continue

        # Column names may be cell.names / copykat.pred
        cell_col = "cell.names" if "cell.names" in df.columns else df.columns[0]
        pred_col = "copykat.pred" if "copykat.pred" in df.columns else df.columns[1]

        n = 0
        for bc, p in zip(df[cell_col].astype(str), df[pred_col].astype(str)):
            norm_p = normalize_pred(p)
            # If the same barcode appears more than once (dataset_pool refs from other patients), keep the first
            if bc not in bc2pred:
                bc2pred[bc] = norm_p
                n += 1
        patient_status.append({"dp_key": dp, "status": "ok", "n_cells_pred": n})

    return bc2pred, pd.DataFrame(patient_status)


def main():
    t0 = datetime.now()
    log(f"Step 4e Collect started: {t0}")

    log(f"Reading {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    log(f"  shape: {adata.shape}")

    # Load predictions
    log("Collecting CopyKAT predictions ...")
    bc2pred, status_df = load_predictions()
    n_ok = int((status_df["status"] == "ok").sum())
    n_fail = int(len(status_df) - n_ok)
    log(f"  patient: ok={n_ok}  fail={n_fail}")
    log(f"Total predicted barcodes: {len(bc2pred):,}")

    # Map onto adata
    log("Map to adata.obs['copykat_pred']")
    pred_arr = np.array([bc2pred.get(bc, "missing") for bc in adata.obs_names], dtype=object)
    adata.obs["copykat_pred"] = pd.Categorical(
        pred_arr,
        categories=["aneuploid", "diploid", "not.defined", "missing"],
    )
    log("copykat_pred distribution:")
    for lab, n in adata.obs["copykat_pred"].value_counts().items():
        log(f"    {lab}: {n:,} ({100*n/adata.n_obs:.2f}%)")

    # Cross-check
    log("\nCross-check: celltype_coarse × copykat_pred")
    ct = pd.crosstab(adata.obs["celltype_coarse"], adata.obs["copykat_pred"])
    log(ct.to_string())

    # Final malignant labels (priority top to bottom)
    log("\nAssign final obs['malignant']")
    ct_coarse = adata.obs["celltype_coarse"].astype(str).values
    copy_pred = adata.obs["copykat_pred"].astype(str).values
    tissue = adata.obs["tissue_type"].astype(str).values

    malignant = np.full(adata.n_obs, "Uncertain", dtype=object)

    # Rule 1: normal tissue → Non-malignant
    malignant[np.isin(tissue, list(NORMAL_TT))] = "Non-malignant"

    # Rule 2: non-Epithelial → Non-malignant (CopyKAT valid for epithelium only)
    non_epi_mask = (ct_coarse != "Epithelial") & (malignant == "Uncertain")
    malignant[non_epi_mask] = "Non-malignant"

    # Rule 3: Epithelial cells by CopyKAT call
    epi_aneu = (ct_coarse == "Epithelial") & (copy_pred == "aneuploid") & (malignant == "Uncertain")
    malignant[epi_aneu] = "Malignant"
    epi_dip = (ct_coarse == "Epithelial") & (copy_pred == "diploid") & (malignant == "Uncertain")
    malignant[epi_dip] = "Non-malignant"
    # Remaining Epithelial (not.defined / missing) stay Uncertain

    adata.obs["malignant"] = pd.Categorical(
        malignant,
        categories=["Malignant", "Non-malignant", "Uncertain"],
    )

    log("malignant distribution:")
    for lab, n in adata.obs["malignant"].value_counts().items():
        log(f"  {lab}: {n:,} ({100*n/adata.n_obs:.2f}%)")

    # Write output
    log(f"\nWriting {OUT_H5AD}")
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log(f" {OUT_H5AD.stat().st_size/1e9:.2f} GB")

    L = ["# Step 4 CopyKAT malignant-cell calling report", ""]
    L.append(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- Total cells: {adata.n_obs:,}")
    L.append(f"- Patients succeeded: {n_ok} / {n_ok + n_fail}")
    L.append("")

    L.append("## Final obs['malignant'] distribution")
    for lab, n in adata.obs["malignant"].value_counts().items():
        L.append(f"- {lab}: {n:,} ({100*n/adata.n_obs:.2f}%)")
    L.append("")

    L.append("## Raw obs['copykat_pred'] distribution")
    for lab, n in adata.obs["copykat_pred"].value_counts().items():
        L.append(f"- {lab}: {n:,} ({100*n/adata.n_obs:.2f}%)")
    L.append("")

    L.append("## celltype_coarse × copykat_pred")
    L.append("```")
    L.append(pd.crosstab(adata.obs["celltype_coarse"], adata.obs["copykat_pred"]).to_string())
    L.append("```")
    L.append("")

    L.append("## celltype_coarse × malignant")
    L.append("```")
    L.append(pd.crosstab(adata.obs["celltype_coarse"], adata.obs["malignant"]).to_string())
    L.append("```")
    L.append("")

    L.append("## tissue_type × malignant")
    tm = pd.crosstab(adata.obs["tissue_type"], adata.obs["malignant"])
    tm_pct = tm.div(tm.sum(axis=1), axis=0) * 100
    L.append("### Counts")
    L.append("```")
    L.append(tm.to_string())
    L.append("```")
    L.append("### Percent")
    L.append("```")
    L.append(tm_pct.round(1).to_string())
    L.append("```")
    L.append("")

    L.append("## dataset × malignant")
    dm = pd.crosstab(adata.obs["dataset"], adata.obs["malignant"])
    dm_pct = dm.div(dm.sum(axis=1), axis=0) * 100
    L.append("### Counts")
    L.append("```")
    L.append(dm.to_string())
    L.append("```")
    L.append("### Percent")
    L.append("```")
    L.append(dm_pct.round(1).to_string())
    L.append("```")
    L.append("")

    L.append("## Sanity checks vs expected ranges")
    tumor_tt = ["Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
                "Distant_Metastasis", "Pleural_Effusion"]
    tumor_epi = adata.obs[
        adata.obs["tissue_type"].isin(tumor_tt) &
        (adata.obs["celltype_coarse"] == "Epithelial")
    ]
    if len(tumor_epi) > 0:
        frac_mal = (tumor_epi["malignant"] == "Malignant").mean()
        L.append(f"- Fraction Malignant among Epithelial in tumour samples: {100*frac_mal:.1f}%  (expect >50%)")
    immune_cells = adata.obs[
        adata.obs["celltype_coarse"].isin(["T_NK", "B", "Myeloid", "Mast", "Plasma"])
    ]
    frac_nm = (immune_cells["malignant"] == "Non-malignant").mean()
    L.append(f"- Fraction Non-malignant among immune cells: {100*frac_nm:.1f}%  (expect >99%; forced by rule 2)")
    total_unc = (adata.obs["malignant"] == "Uncertain").mean()
    L.append(f"- Overall Uncertain fraction: {100*total_unc:.2f}%  (expect <10%)")
    L.append("")

    L.append("## Malignant fraction of Epithelial by tissue (key metric)")
    L.append("```")
    epi_by_tt = []
    for tt in sorted(adata.obs["tissue_type"].unique()):
        sub = adata.obs[(adata.obs["tissue_type"] == tt) & (adata.obs["celltype_coarse"] == "Epithelial")]
        if len(sub) == 0:
            continue
        vc = sub["malignant"].value_counts()
        n_mal = int(vc.get("Malignant", 0))
        n_nm  = int(vc.get("Non-malignant", 0))
        n_unc = int(vc.get("Uncertain", 0))
        tot = len(sub)
        epi_by_tt.append([tt, tot, n_mal, n_nm, n_unc,
                          f"{100*n_mal/tot:.1f}%", f"{100*n_unc/tot:.1f}%"])
    header = ["tissue_type", "epi_total", "Malignant", "Non-malignant", "Uncertain", "Mal%", "Unc%"]
    L.append(pd.DataFrame(epi_by_tt, columns=header).to_string(index=False))
    L.append("```")
    L.append("")

    L.append("## Malignant cell counts by dataset")
    L.append("```")
    ds_mal = adata.obs.groupby("dataset", observed=True)["malignant"].apply(
        lambda s: (s == "Malignant").sum()
    )
    L.append(ds_mal.to_string())
    L.append("```")
    L.append("")

    if n_fail > 0:
        L.append("## Failed / missing patients")
        L.append("```")
        L.append(status_df[status_df["status"] != "ok"].to_string(index=False))
        L.append("```")

    REPORT.write_text("\n".join(L), encoding="utf-8")
    log(f"Report: {REPORT}")

    elapsed = (datetime.now() - t0).total_seconds() / 60
    log(f"\nElapsed: {elapsed:.1f} min")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
