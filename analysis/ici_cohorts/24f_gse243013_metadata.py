"""Step 24f: GSE243013 (n=243 NSCLC neoadj chemo-IO scRNA) — metadata-only.

Counts mtx (6.6 GB) NOT downloaded; only cell-level metadata + UMAP info available.
Metadata is immune-cell scRNA only, so MP1-4 (malignant programs) NOT scoreable.

Outputs (fig_treatment/):
  gse243013_sample_metadata.csv     243 sample × clinical
  gse243013_celltype_composition.csv  243 sample × major/sub cell type fraction
  gse243013_response_celltype_compare.csv  per cell type, MPR vs non-MPR Wilcoxon
  gse243013_skip_note.md
"""
from __future__ import annotations
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(os.path.expanduser("~/luad/scripts")))
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from importlib import import_module
C = import_module("24_common")

GSE = "gse243013"
META = Path("${DATA_ROOT}/GSE243013/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")


def main():
    C.log(f"=== {GSE.upper()} metadata-only ===")
    md = pd.read_csv(META, low_memory=False)
    C.log(f"cell-level rows: {len(md)}")

    # Sample-level collapse
    samp = md.groupby("sampleID").agg(
        cancer_type=("cancer_type", "first"),
        gender=("gender", "first"),
        age=("age", "first"),
        smoking=("smoking_history", "first"),
        stage=("pre_treatment_staging", "first"),
        anti_PD1=("anti-PD1_therapy", "first"),
        chemotherapy=("chemotherapy", "first"),
        targeted=("targeted_therapy", "first"),
        cycles=("cycles", "first"),
        pathological_response=("pathological_response", "first"),
        path_response_rate=("pathological_response_rate", "first"),
        radiological_response=("radiological_response", "first"),
        n_cells=("cellID", "count"),
    ).reset_index()
    # 3-way response
    samp["response_3way"] = samp["pathological_response"]
    samp["response_binary"] = samp["pathological_response"].map(
        lambda x: "R" if x in ("pCR", "MPR") else
                  ("NR" if x == "non-MPR" else "unknown"))
    C.log(f"sample n: {len(samp)}; cancer: {samp['cancer_type'].value_counts().to_dict()}")
    C.log(f"3-way response: {samp['response_3way'].value_counts().to_dict()}")
    C.log(f"binary R/NR: {samp['response_binary'].value_counts().to_dict()}")

    # Cell type composition (major + sub)
    counts = (md.groupby(["sampleID", "major_cell_type"])
                .size().unstack(fill_value=0))
    frac = counts.div(counts.sum(axis=1), axis=0)
    frac.columns = ["frac_" + str(c).replace(" ", "_").replace("/", "_")
                     for c in frac.columns]
    sub_counts = (md.groupby(["sampleID", "sub_cell_type"])
                    .size().unstack(fill_value=0))
    sub_frac = sub_counts.div(sub_counts.sum(axis=1), axis=0)
    sub_frac.columns = ["sub_frac_" + str(c).replace(" ", "_").replace("/", "_")
                         for c in sub_frac.columns]
    comp = frac.join(sub_frac, how="outer").fillna(0).reset_index()
    comp_full = samp.merge(comp, on="sampleID", how="left")

    OUT = C.OUT; OUT.mkdir(parents=True, exist_ok=True)
    samp.to_csv(OUT / f"{GSE}_sample_metadata.csv", index=False)
    comp_full.to_csv(OUT / f"{GSE}_celltype_composition.csv", index=False)

    # Response comparison on each cell type fraction
    cmp_full = comp_full[comp_full["response_binary"].isin(["R","NR"])].copy()
    rrows = []
    for col in [c for c in cmp_full.columns if c.startswith("frac_")
                  or c.startswith("sub_frac_")]:
        v_R = cmp_full.loc[cmp_full["response_binary"]=="R", col].dropna().values
        v_NR = cmp_full.loc[cmp_full["response_binary"]=="NR", col].dropna().values
        if len(v_R) < 5 or len(v_NR) < 5:
            continue
        try:
            U, p = mannwhitneyu(v_R, v_NR, alternative="two-sided")
        except Exception:
            U = p = np.nan
        rrows.append({
            "feature": col,
            "n_R": int(len(v_R)), "n_NR": int(len(v_NR)),
            "median_R": float(np.median(v_R)),
            "median_NR": float(np.median(v_NR)),
            "delta_RN": float(np.median(v_R) - np.median(v_NR)),
            "U": float(U), "p": float(p),
        })
    cmp_df = pd.DataFrame(rrows).sort_values("p")
    cmp_df.to_csv(OUT / f"{GSE}_response_celltype_compare.csv", index=False)
    C.log(f"\nTop 10 cell-type features by p (R vs NR):")
    print(cmp_df.head(10).round(4).to_string(index=False))

    note = (f"# GSE243013 — metadata-only analysis\n\n"
            f"- Counts mtx (6.6 GB) NOT downloaded; expression-based MP scoring NOT done.\n"
            f"- This cohort is immune-cell scRNA-seq (T/NK, B, myeloid). No malignant cells.\n"
            f"- MP1-4 are malignant-cell programs → cannot score on immune cells.\n\n"
            f"## Cohort summary\n"
            f"- Samples: {len(samp)} (LUAD: {sum(samp['cancer_type']=='LUAD')}; "
            f"LUSC: {sum(samp['cancer_type']=='LUSC')})\n"
            f"- Total cells: {samp['n_cells'].sum():,}\n"
            f"- 3-way pathological response: "
            f"{samp['response_3way'].value_counts().to_dict()}\n"
            f"- Binary R/NR: {samp['response_binary'].value_counts().to_dict()}\n\n"
            f"## What's available\n"
            f"- Sample × cell type composition vs response (R vs NR Wilcoxon)\n"
            f"- See `{GSE}_response_celltype_compare.csv` for per-feature stats\n")
    (OUT / f"{GSE}_skip_note.md").write_text(note, encoding="utf-8")
    C.log("written outputs")


if __name__ == "__main__":
    main()
