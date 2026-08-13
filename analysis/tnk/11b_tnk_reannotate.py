"""Step 11b: Re-annotate T_NK subtypes using explicit per-cluster expression rules.

Why: initial score_genes panel voting failed because core lineage markers
(CD8A, CD8B, FOXP3, PDCD1) are only in obsm['marker_counts'], not in var.
The T_general panel (CD3D/E/G — all in var) swamped CD4/CD8 panels, so no
cluster was labeled CD8_generic. gdT was over-assigned (threshold too lax).

This script:
  1. Loads luad_tnk.h5ad (keeps existing Leiden 1.5, UMAP, Harmony).
  2. Builds per-cluster mean expression combining var + marker_counts.
  3. Assigns tnk_subtype via explicit rule cascade (not panel scoring).
  4. Re-exports all downstream Fig_tnk artifacts.
"""

from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

IN = Path.home() / "luad/data/processed/luad_tnk.h5ad"
MAL = Path.home() / "luad/data/processed/luad_malignant_scored.h5ad"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig_tnk")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    # 1) Build combined marker-mean matrix per cluster (var ∪ marker_counts)
    cluster_col = "leiden_1.5"
    mc_df = a.obsm.get("marker_counts", None)

    var_markers = ["CD4", "CD3D", "CD3E", "CD3G", "IL7R", "SELL", "TCF7",
                   "NKG7", "GNLY", "GZMB", "GZMA", "GZMK", "PRF1", "CCL5",
                   "KLRD1", "KLRB1", "FGFBP2", "FCGR3A",
                   "ITGAE", "CD69", "RUNX3",
                   "CTLA4", "TIGIT", "HAVCR2", "LAG3", "IL2RA",
                   "TOP2A", "STMN1", "HMGB2", "TRGC2", "ZBTB16"]
    mc_markers = ["CD8A", "CD8B", "FOXP3", "PDCD1"]  # only in marker_counts

    var_present = [g for g in var_markers if g in a.var_names]
    mc_present = [g for g in mc_markers if mc_df is not None and g in mc_df.columns]
    log(f"  var_present: {len(var_present)}  mc_present: {len(mc_present)}")

    # Build dense expression DF for all cells
    X_var = a[:, var_present].X
    X_var = X_var.toarray() if hasattr(X_var, "toarray") else np.asarray(X_var)
    df = pd.DataFrame(X_var, index=a.obs.index, columns=var_present)

    if mc_present:
        # marker_counts is raw counts → apply log1p for comparability
        mc_raw = mc_df[mc_present].astype("float32")
        mc_log = np.log1p(mc_raw)
        for g in mc_present:
            df[g] = mc_log[g].values

    df["cluster"] = a.obs[cluster_col].astype(str).values
    cluster_mean = df.groupby("cluster").mean()
    log(f"  cluster_mean shape: {cluster_mean.shape}")

    cluster_mean.to_csv(RES / "step11b_cluster_marker_mean.csv")

    # 2) Derive subtype per cluster via explicit rules
    def cd8_score(row):
        return max(row.get("CD8A", 0), row.get("CD8B", 0))

    def exhaust_score(row):
        return (row.get("PDCD1", 0) + row.get("HAVCR2", 0) + row.get("LAG3", 0)) / 3

    def cd8_dom(row, ratio=1.5, floor=0.1):
        cd8 = cd8_score(row); cd4 = row.get("CD4", 0)
        return cd8 >= floor and cd8 >= ratio * max(cd4, 0.05)

    def cd4_dom(row, ratio=1.5, floor=0.2):
        cd8 = cd8_score(row); cd4 = row.get("CD4", 0)
        return cd4 >= floor and cd4 >= ratio * max(cd8, 0.05)

    decisions = {}
    decision_reasons = {}
    for cl, row in cluster_mean.iterrows():
        # 1. Proliferating (cell-cycle signature dominates, not NK)
        if (row.get("TOP2A", 0) >= 0.5 or row.get("STMN1", 0) >= 1.0
            or row.get("HMGB2", 0) >= 1.5) \
           and row.get("NKG7", 0) < 3.0 and row.get("CD3D", 0) >= 0.5:
            decisions[cl] = "T_Proliferating"
            decision_reasons[cl] = f"TOP2A={row.get('TOP2A',0):.2f}"
            continue

        # 2. NK — NKG7 high, CD3D low, CD8 tiny
        if (row.get("NKG7", 0) >= 2.5 and row.get("CD3D", 0) < 0.8
            and cd8_score(row) < 0.2):
            decisions[cl] = "NK"
            decision_reasons[cl] = f"NKG7={row.get('NKG7',0):.2f} CD3D={row.get('CD3D',0):.2f}"
            continue

        # 3. Treg — strict FOXP3 + IL2RA co-expression
        if (row.get("FOXP3", 0) >= 0.3 and row.get("IL2RA", 0) >= 0.2
            and row.get("CD4", 0) >= 0.2 and not cd8_dom(row)):
            decisions[cl] = "Treg"
            decision_reasons[cl] = (f"FOXP3={row.get('FOXP3',0):.2f} "
                                     f"IL2RA={row.get('IL2RA',0):.2f}")
            continue

        # 4. gdT — TRGC2 very high + CD4/CD8 low
        if (row.get("TRGC2", 0) >= 2.0 and cd8_score(row) < 0.3
            and row.get("CD4", 0) < 0.3):
            decisions[cl] = "gdT"
            decision_reasons[cl] = f"TRGC2={row.get('TRGC2',0):.2f}"
            continue

        # 5. NKT
        if (row.get("ZBTB16", 0) >= 1.0 and row.get("NKG7", 0) >= 2.0
            and row.get("CD3D", 0) >= 1.0):
            decisions[cl] = "NKT"
            decision_reasons[cl] = f"ZBTB16={row.get('ZBTB16',0):.2f}"
            continue

        # 6. CD8 lineage branches (exhausted BEFORE TRM; LUAD CD8 is
        # often TRM+exhausted overlap, and LAG3 is the most reliable
        # exhaustion marker here since PDCD1 capture is sparse)
        if cd8_dom(row):
            # 6a. Exhausted — combined checkpoint score OR LAG3 high + CD3+
            if exhaust_score(row) >= 0.4 or \
               (row.get("LAG3", 0) >= 0.8 and row.get("HAVCR2", 0) >= 0.2):
                decisions[cl] = "CD8_Exhausted"
                decision_reasons[cl] = (f"ExScore={exhaust_score(row):.2f} "
                                         f"LAG3={row.get('LAG3',0):.2f}")
                continue
            # 6b. TRM
            if row.get("ITGAE", 0) >= 0.2 or \
               (row.get("CD69", 0) >= 1.0 and row.get("RUNX3", 0) >= 1.5):
                decisions[cl] = "CD8_TRM"
                decision_reasons[cl] = (f"ITGAE={row.get('ITGAE',0):.2f} "
                                         f"CD69={row.get('CD69',0):.2f}")
                continue
            # 6c. Effector (strong GZMB/PRF1)
            if (row.get("GZMB", 0) >= 1.5 or row.get("PRF1", 0) >= 1.5
                or row.get("GNLY", 0) >= 1.5):
                decisions[cl] = "CD8_Effector"
                decision_reasons[cl] = f"GZMB={row.get('GZMB',0):.2f}"
                continue
            # 6d. Naive/CM
            if (row.get("SELL", 0) >= 0.8 or row.get("TCF7", 0) >= 0.6) \
               and row.get("GZMB", 0) < 1.0:
                decisions[cl] = "CD8_Naive_CM"
                decision_reasons[cl] = f"SELL={row.get('SELL',0):.2f}"
                continue
            decisions[cl] = "CD8_other"
            decision_reasons[cl] = f"CD8+={cd8_score(row):.2f} catchall"
            continue

        # 7. CD4 lineage
        if cd4_dom(row):
            if (row.get("SELL", 0) >= 0.5 or row.get("TCF7", 0) >= 0.6) \
               and row.get("GZMB", 0) < 1.0:
                decisions[cl] = "CD4_Naive_CM"
                decision_reasons[cl] = f"CD4={row.get('CD4',0):.2f} naive"
                continue
            decisions[cl] = "CD4_generic"
            decision_reasons[cl] = f"CD4={row.get('CD4',0):.2f}"
            continue

        # 8. Truly ambiguous — catchall
        # Try to recover via NKG7+CD3D pattern (cytotoxic T of unclear lineage)
        if row.get("NKG7", 0) >= 1.5 and row.get("CD3D", 0) >= 1.0:
            decisions[cl] = "T_cytotoxic_ambiguous"
            decision_reasons[cl] = "NKG7+/CD3+ unclear lineage"
            continue

        decisions[cl] = "T_unresolved"
        decision_reasons[cl] = (f"CD4={row.get('CD4',0):.2f} "
                                 f"CD8={cd8_score(row):.2f}")

    # Save cluster decision table
    dec_rows = []
    for cl in sorted(decisions, key=lambda x: int(x)):
        row = cluster_mean.loc[cl]
        dec_rows.append({
            "cluster": cl,
            "n_cells": int((a.obs[cluster_col].astype(str) == cl).sum()),
            "subtype": decisions[cl],
            "reason": decision_reasons[cl],
            "CD8A": row.get("CD8A", np.nan),
            "CD8B": row.get("CD8B", np.nan),
            "CD4":  row.get("CD4",  np.nan),
            "CD3D": row.get("CD3D", np.nan),
            "NKG7": row.get("NKG7", np.nan),
            "FOXP3": row.get("FOXP3", np.nan),
            "IL2RA": row.get("IL2RA", np.nan),
            "PDCD1": row.get("PDCD1", np.nan),
            "HAVCR2": row.get("HAVCR2", np.nan),
            "LAG3": row.get("LAG3", np.nan),
            "GZMB": row.get("GZMB", np.nan),
            "GZMK": row.get("GZMK", np.nan),
            "SELL": row.get("SELL", np.nan),
            "ITGAE": row.get("ITGAE", np.nan),
            "CD69": row.get("CD69", np.nan),
            "TRGC2": row.get("TRGC2", np.nan),
            "ZBTB16": row.get("ZBTB16", np.nan),
            "TOP2A": row.get("TOP2A", np.nan),
        })
    dec_df = pd.DataFrame(dec_rows)
    dec_df.to_csv(RES / "step11b_cluster_decisions.csv", index=False)
    log("cluster decisions (abridged):")
    log(dec_df[["cluster", "n_cells", "subtype", "reason"]].to_string(index=False))

    # 3) Map to obs
    a.obs["tnk_subtype"] = (
        a.obs[cluster_col].astype(str).map(decisions).astype("category")
    )
    counts = a.obs["tnk_subtype"].value_counts()
    log("tnk_subtype counts:")
    log(counts.to_string())

    # 4) Re-export proportions
    TISSUE_ORDER = ["Precancerous", "Adjacent_Normal", "Normal_Lung", "Normal_LN",
                    "Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
                    "Distant_Metastasis", "Pleural_Effusion"]
    prop_tis = (pd.crosstab(a.obs["tissue_type"], a.obs["tnk_subtype"],
                             normalize="index")
                 .reindex(index=[t for t in TISSUE_ORDER if t in
                                  a.obs["tissue_type"].unique()]))
    prop_tis.to_csv(FIG / "tnk_proportion_by_tissue.csv")

    prop_ds = pd.crosstab(a.obs["dataset"], a.obs["tnk_subtype"], normalize="index")
    prop_ds.to_csv(FIG / "tnk_proportion_by_dataset.csv")

    # 5) Dotplot (all markers, subtype × gene)
    all_markers = sorted(set(var_present + mc_present))
    dot_rows = []
    for st in a.obs["tnk_subtype"].cat.categories:
        mask = (a.obs["tnk_subtype"] == st).values
        if not mask.any(): continue
        subdf = df.loc[mask, all_markers]
        mean_log = subdf.mean(axis=0)
        pct = (subdf > 0).mean(axis=0)
        for g in all_markers:
            dot_rows.append({"subtype": st, "gene": g,
                             "mean_log1p": float(mean_log[g]),
                             "pct_expressing": float(pct[g])})
    pd.DataFrame(dot_rows).to_csv(FIG / "tnk_dotplot_markers.csv", index=False)

    # 6) Exhausted vs Treg profile
    chk_genes = ["CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4",
                 "PDCD1", "HAVCR2", "LAG3", "TIGIT"]
    chk_present = [g for g in chk_genes if g in df.columns]
    profile_rows = []
    for st in a.obs["tnk_subtype"].cat.categories:
        mask = (a.obs["tnk_subtype"] == st).values
        if not mask.any(): continue
        sub = df.loc[mask, chk_present]
        m = sub.mean(axis=0)
        pct = (sub > 0).mean(axis=0)
        for g in chk_present:
            profile_rows.append({"subtype": st, "gene": g,
                                 "mean_log1p": float(m[g]),
                                 "pct_expressing": float(pct[g])})
    pd.DataFrame(profile_rows).to_csv(FIG / "tnk_exhausted_treg_profile.csv",
                                        index=False)

    # 7) Patient-level × MP Spearman
    log("patient-level MP correlation")
    mal = sc.read_h5ad(MAL, backed="r")
    mp_pat = (mal.obs[["patient_id", "MP1_score", "MP2_score",
                         "MP3_score", "MP4_score"]]
                  .groupby("patient_id", observed=True).mean())
    tnk_pat = pd.crosstab(a.obs["patient_id"], a.obs["tnk_subtype"],
                            normalize="index")
    joint = tnk_pat.join(mp_pat, how="inner")
    rows = []
    for st in tnk_pat.columns:
        for mpx in ["MP1_score", "MP2_score", "MP3_score", "MP4_score"]:
            s = joint[[st, mpx]].dropna()
            if len(s) < 5: continue
            rho, p = spearmanr(s[st], s[mpx])
            rows.append({"subtype": st, "MP": mpx.replace("_score", ""),
                         "n_patients": len(s),
                         "spearman_rho": rho, "p": p})
    mp_corr = pd.DataFrame(rows)
    mp_corr.to_csv(RES / "step11_tnk_mp3_correlation.csv", index=False)
    mp_corr.to_csv(FIG / "tnk_mp_association.csv", index=False)

    # 8) UMAP metadata
    umap_df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "patient_id": a.obs["patient_id"].astype(str).values,
        "tissue_type": a.obs["tissue_type"].astype(str).values,
        "leiden_1.0": a.obs["leiden_1.0"].astype(str).values,
        "leiden_1.5": a.obs["leiden_1.5"].astype(str).values,
        "tnk_subtype": a.obs["tnk_subtype"].astype(str).values,
        "UMAP1": a.obsm["X_umap"][:, 0],
        "UMAP2": a.obsm["X_umap"][:, 1],
    })
    umap_df.to_csv(FIG / "tnk_umap_metadata.csv.gz", index=False,
                    compression="gzip")

    # 9) Save updated h5ad (in place)
    log("writing back h5ad with new tnk_subtype")
    a.obs[["leiden_1.0", "leiden_1.5", "tnk_subtype"]].to_csv(
        RES / "step11b_obs_labels.csv.gz", compression="gzip"
    )
    a.write(IN)

    # 10) Summary md
    with open(RES / "step11_tnk_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 11b — T_NK explicit-rule annotation\n\n")
        f.write(f"- Total T_NK cells: {a.n_obs}\n")
        f.write(f"- Leiden 1.5 clusters: {a.obs[cluster_col].nunique()}\n")
        f.write(f"- Subtypes: {a.obs['tnk_subtype'].nunique()}\n\n")
        f.write("## Why step11b was needed\n\n")
        f.write("In Step 11 the score_genes panel vote failed because CD8A/CD8B/"
                "FOXP3/PDCD1 are only in obsm['marker_counts'] (not in var). The "
                "T_general (CD3D/E/G) panel out-scored CD4/CD8 lineage panels, so "
                "no CD8 cluster was identified and gdT was over-assigned (threshold "
                "too lax). Step 11b replaces panel voting with explicit per-cluster "
                "expression rules combining var + marker_counts.\n\n")
        f.write("## Subtype cell counts\n\n")
        f.write(counts.to_frame("n_cells").to_markdown() + "\n\n")
        f.write("## Cluster decision table (Leiden 1.5)\n\n")
        f.write(dec_df.drop(columns=["GZMK"]).round(3).to_markdown(index=False) + "\n\n")
        f.write("## Top MP ↔ subtype % Spearman (patient-level, n=88)\n\n")
        top = (mp_corr.sort_values("spearman_rho", ascending=False,
                                    key=lambda s: s.abs()).head(20))
        f.write(top.round(4).to_markdown(index=False) + "\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
