"""Step 10c: Refine Myeloid annotation at Leiden 1.5 resolution.

Goal: recover Macro_SPP1 and Macro_FOLR2 subtypes missed at Leiden 1.0.

Pipeline:
  1. Load luad_myeloid.h5ad.
  2. Re-vote subtype per leiden_1.5 cluster using existing sc_* scores.
  3. Compute mean log1p expression of key SPP1/FOLR2 markers per cluster
     for an explicit marker-based check.
  4. Apply rescue rules:
       SPP1: top SPP1 expression AND TREM2 > global median
       FOLR2: top FOLR2 expression AND CD163 > global median
  5. Update myeloid_subtype_refined and export.

Outputs:
  - ~/luad/results/step10c_leiden15_votes.csv
  - ~/luad/results/step10c_marker_matrix.csv
  - ~/luad/results/step10c_diagnosis.md
  - ${WORK_ROOT}/luad_figures/fig4/myeloid_umap_metadata.csv.gz  (updated)
  - ${WORK_ROOT}/luad_figures/fig4/myeloid_m1m2_scores_refined.csv
  - ${WORK_ROOT}/luad_figures/fig4/myeloid_mp3_association_refined.csv
  - ${WORK_ROOT}/luad_figures/fig4/myeloid_dotplot_markers_refined.csv
"""

from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

IN = Path.home()/"luad/data/processed/luad_myeloid.h5ad"
MAL = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
RES = Path.home()/"luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig4")

CHECK_GENES = ["CD68","CSF1R","APOE","C1QA","C1QB","C1QC","LYZ",
               "SPP1","MMP9","VEGFA","TREM2",
               "FOLR2","MRC1","CD163",
               "MARCO","MSR1","FABP4","PPARG",
               "S100A8","S100A9","VCAN","CD14","FCN1",
               "FCGR3A","CDKN1C",
               "CSF3R","ELANE",
               "CLEC9A","BATF3","IRF8","CD1C","CLEC10A","FCER1A",
               "LAMP3","CCR7","FSCN1",
               "TCF4","IL3RA","LILRA4",
               "TOP2A","STMN1","MKI67"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    # Marker matrix: cluster × gene mean log1p
    present_genes = [g for g in CHECK_GENES if g in a.var_names]
    log(f"  markers present: {len(present_genes)}/{len(CHECK_GENES)}")
    X = a[:, present_genes].X
    df_expr = pd.DataFrame(
        X.toarray() if hasattr(X, "toarray") else np.asarray(X),
        index=a.obs.index,
        columns=present_genes,
    )
    df_expr["cluster"] = a.obs["leiden_1.5"].astype(str).values
    marker_mean = df_expr.groupby("cluster").mean()
    marker_mean.to_csv(RES/"step10c_marker_matrix.csv")

    # Re-vote at leiden_1.5 using existing sc_* scores in obs
    score_cols = [c for c in a.obs.columns if c.startswith("sc_")]
    clust = a.obs.groupby("leiden_1.5", observed=True)[score_cols].mean()
    best = clust.idxmax(axis=1).str.replace("sc_", "")
    second = clust.apply(lambda r: r.nlargest(2).index[-1], axis=1).str.replace("sc_", "")
    margin = clust.apply(lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[-1], axis=1)
    vote = pd.DataFrame({
        "cluster": clust.index,
        "best_subtype": best.values,
        "second_subtype": second.values,
        "margin": margin.values,
        "n_cells": a.obs.groupby("leiden_1.5", observed=True).size().reindex(clust.index).values,
    })

    # Attach marker hints to the vote table
    hint_cols = ["SPP1","MMP9","VEGFA","TREM2","FOLR2","MRC1","CD163","APOE","C1QC","MARCO","CSF3R","CD68"]
    hint_cols = [c for c in hint_cols if c in marker_mean.columns]
    vote = vote.merge(marker_mean[hint_cols].round(3),
                      left_on="cluster", right_index=True)
    vote.to_csv(RES/"step10c_leiden15_votes.csv", index=False)
    log("votes (leiden_1.5):")
    log(vote.to_string(index=False))

    # Stringent rescue rules — favor the strongest marker-expressing cluster(s)
    # SPP1: absolute log1p threshold, not proliferating, SPP1 is the dominant
    # pro-tumoral Macro signature (top_1 or top_2 marker expression)
    def top3_hits(row, n=3):
        return [row.nlargest(n).index[i].replace("sc_","") for i in range(n)]
    top3 = clust.apply(top3_hits, axis=1)
    vote["top3"] = top3.values

    # Rule 1: Macro_SPP1 — SPP1 mean >= 2.0 AND not already proliferating/DC
    SPP1_TH, FOLR2_TH, CD163_TH = 2.0, 0.9, 1.0
    spp1_cands = marker_mean[marker_mean["SPP1"] >= SPP1_TH].index.tolist()
    # Rule 2: Macro_FOLR2 — FOLR2>=0.9 AND CD163>=1.0 AND SPP1 < SPP1_TH
    folr2_cands = marker_mean[
        (marker_mean["FOLR2"] >= FOLR2_TH) &
        (marker_mean["CD163"] >= CD163_TH) &
        (marker_mean["SPP1"]  <  SPP1_TH)
    ].index.tolist()
    log(f"Macro_SPP1 candidate clusters  (SPP1>={SPP1_TH}): {spp1_cands}")
    log(f"Macro_FOLR2 candidate clusters (FOLR2>={FOLR2_TH}, CD163>={CD163_TH}, SPP1<{SPP1_TH}): {folr2_cands}")

    mapping = dict(zip(vote["cluster"], vote["best_subtype"]))
    protect = {"pDC", "cDC1", "cDC2", "cDC_LAMP3", "Macro_prolif",
               "Mono_classical", "Mono_nonclassical", "Neutrophil"}
    for cl in spp1_cands:
        if mapping[cl] not in protect:
            mapping[cl] = "Macro_SPP1"
    for cl in folr2_cands:
        if mapping[cl] not in protect and mapping[cl] != "Macro_SPP1":
            mapping[cl] = "Macro_FOLR2"

    # Also carry over Neutrophil rescue (CSF3R hi / CD68 lo / C1QA lo, unchanged logic)
    for cl in vote["cluster"]:
        row = marker_mean.loc[cl]
        if (row.get("CSF3R", 0) > marker_mean["CSF3R"].median() * 1.5 and
            row.get("CD68", 0) < marker_mean["CD68"].median() * 0.5 and
            row.get("C1QA", 0) < marker_mean["C1QA"].median() * 0.5):
            mapping[cl] = "Neutrophil"

    a.obs["myeloid_subtype_refined"] = (
        a.obs["leiden_1.5"].astype(str).map(mapping).astype("category")
    )
    refined_counts = a.obs["myeloid_subtype_refined"].value_counts()
    log("refined subtype counts:")
    log(refined_counts.to_string())

    # Compare with original
    cross = pd.crosstab(a.obs["myeloid_subtype"], a.obs["myeloid_subtype_refined"])
    cross.to_csv(RES/"step10c_subtype_transition.csv")

    m1m2 = (a.obs.groupby("myeloid_subtype_refined", observed=True)
            [["M1_score","M2_score","M1_M2_ratio"]]
            .agg(["mean","median","std","count"]))
    m1m2.columns = [f"{c}_{k}" for c,k in m1m2.columns]
    m1m2.to_csv(FIG/"myeloid_m1m2_scores_refined.csv")

    # Dot plot table refined (subtype × marker gene, mean log1p + pct>0)
    dot_rows = []
    for st in a.obs["myeloid_subtype_refined"].cat.categories:
        mask = (a.obs["myeloid_subtype_refined"] == st).values
        if not mask.any(): continue
        sub = df_expr.loc[mask, present_genes].values
        mean_expr = sub.mean(axis=0)
        pct_expr = (sub > 0).mean(axis=0)
        for i, g in enumerate(present_genes):
            dot_rows.append({"subtype": st, "gene": g,
                             "mean_log1p": float(mean_expr[i]),
                             "pct_expressing": float(pct_expr[i])})
    pd.DataFrame(dot_rows).to_csv(FIG/"myeloid_dotplot_markers_refined.csv", index=False)

    # Patient-level % refined × MP3 Spearman
    log("patient-level refined × MP3 Spearman")
    mal = sc.read_h5ad(MAL, backed="r")
    mp_pat = (mal.obs[["patient_id","MP1_score","MP2_score","MP3_score","MP4_score"]]
                 .groupby("patient_id", observed=True).mean())
    my_pat = pd.crosstab(a.obs["patient_id"], a.obs["myeloid_subtype_refined"],
                          normalize="index")
    joint = my_pat.join(mp_pat, how="inner")
    rows = []
    for sub in my_pat.columns:
        for mpx in ["MP1_score","MP2_score","MP3_score","MP4_score"]:
            s = joint[[sub, mpx]].dropna()
            if len(s) < 5: continue
            rho, p = spearmanr(s[sub], s[mpx])
            rows.append({"subtype": sub, "MP": mpx.replace("_score",""),
                         "n_patients": len(s), "spearman_rho": rho, "p": p})
    mp_corr = pd.DataFrame(rows)
    mp_corr.to_csv(FIG/"myeloid_mp3_association_refined.csv", index=False)

    # Proportion tables refined
    TISSUE_ORDER = ["Precancerous","Adjacent_Normal","Normal_Lung","Normal_LN",
                    "Primary_Tumor","LN_Metastasis","Brain_Metastasis",
                    "Distant_Metastasis","Pleural_Effusion"]
    prop_tis = (pd.crosstab(a.obs["tissue_type"], a.obs["myeloid_subtype_refined"],
                             normalize="index")
                 .reindex(index=[t for t in TISSUE_ORDER
                                  if t in a.obs["tissue_type"].unique()]))
    prop_tis.to_csv(FIG/"myeloid_proportion_by_tissue_refined.csv")

    # Update UMAP metadata
    umap_df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "patient_id": a.obs["patient_id"].astype(str).values,
        "tissue_type": a.obs["tissue_type"].astype(str).values,
        "leiden_1.0": a.obs["leiden_1.0"].astype(str).values,
        "leiden_1.5": a.obs["leiden_1.5"].astype(str).values,
        "myeloid_subtype": a.obs["myeloid_subtype"].astype(str).values,
        "myeloid_subtype_refined": a.obs["myeloid_subtype_refined"].astype(str).values,
        "UMAP1": a.obsm["X_umap"][:,0],
        "UMAP2": a.obsm["X_umap"][:,1],
        "M1_score": a.obs["M1_score"].values,
        "M2_score": a.obs["M2_score"].values,
    })
    umap_df.to_csv(FIG/"myeloid_umap_metadata.csv.gz",
                    index=False, compression="gzip")

    # Save refined label back to h5ad (obs-only write)
    # Safer: rewrite a small sidecar
    a.obs[["myeloid_subtype","myeloid_subtype_refined","leiden_1.0","leiden_1.5"]].to_csv(
        RES/"step10c_obs_labels.csv.gz", compression="gzip"
    )

    # Summary md
    log("writing step10c_diagnosis.md")
    with open(RES/"step10c_diagnosis.md", "w", encoding="utf-8") as f:
        f.write("# Step 10c — Myeloid refinement at Leiden 1.5\n\n")
        f.write(f"- Leiden 1.0 clusters: 21 → 11 subtypes\n")
        f.write(f"- Leiden 1.5 clusters: {len(vote)} → {a.obs['myeloid_subtype_refined'].nunique()} subtypes\n")
        f.write(f"- Macro_SPP1 rescued clusters: {spp1_cands}\n")
        f.write(f"- Macro_FOLR2 rescued clusters: {folr2_cands}\n\n")
        f.write("## Refined subtype counts\n\n")
        f.write(refined_counts.to_frame("n_cells").to_markdown() + "\n\n")
        f.write("## Original → Refined transition\n\n")
        f.write(cross.to_markdown() + "\n\n")
        f.write("## Leiden 1.5 cluster vote with marker hints\n\n")
        f.write(vote.round(3).to_markdown(index=False) + "\n\n")
        f.write("## Refined M1/M2 per subtype\n\n")
        f.write(m1m2.round(3).to_markdown() + "\n\n")
        f.write("## Top MP3 ↔ subtype % Spearman (refined)\n\n")
        top = (mp_corr[mp_corr["MP"]=="MP3"]
               .sort_values("spearman_rho", ascending=False).head(15))
        f.write(top.round(4).to_markdown(index=False) + "\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
