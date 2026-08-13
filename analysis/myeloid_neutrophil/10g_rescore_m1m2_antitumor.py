"""Step 10g: Rescore M1/M2 polarization + Anti-tumor panel for Fig 4F/4M.

Fixes a bug in 10a/10f where missing genes (CXCL9/10/11) from marker_counts
were combined via raw-z-score (with NaNs) producing M1_score range up to 66.

Strategy (clean):
  1. Recover CXCL9/10/11 from obsm['marker_counts'] via log1p(raw_count) imputed
     with 0 where dataset lacks the gene.
  2. Concatenate to lognorm X as 3 extra synthetic var slots.
  3. Use sc.tl.score_genes on the augmented matrix uniformly.

Outputs (overwrites):
  ${WORK_ROOT}/luad_figures/fig4/myeloid_m1m2_scores.csv
  ${WORK_ROOT}/luad_figures/fig4/myeloid_m1m2_scores_refined.csv
  ${WORK_ROOT}/luad_figures/fig4/panel_F_antitumor_genes.csv
  ${WORK_ROOT}/luad_figures/fig4/panel_major_type_metadata.csv.gz
    (refresh M1/M2 cols only)
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

H5 = Path.home() / "luad/data/processed/luad_myeloid.h5ad"
LBL = Path.home() / "luad/results/step10c_obs_labels.csv.gz"
FIG = Path("${WORK_ROOT}/luad_figures/fig4")

# Major-type mapping (per Figure_4_plot.R spec)
MAJ_MAP = {
    "Macro_C1QC": "Macrophage", "Macro_FCN1": "Macrophage",
    "Macro_FOLR2": "Macrophage", "Macro_MARCO": "Macrophage",
    "Macro_SPP1": "Macrophage", "Macro_general": "Macrophage",
    "Macro_prolif": "Macrophage",
    "Mono_nonclassical": "Mono_nonclassical", "Mono_NC": "Mono_nonclassical",
    "Neutrophil": "Neutrophil",
    "cDC1": "cDC1", "cDC2": "cDC2",
    "cDC_LAMP3": "cDC_LAMP3", "pDC": "pDC",
}

M1_MARKERS = ["TNF", "IL1B", "CXCL10", "CXCL9", "CXCL11"]
M2_MARKERS = ["CD163", "MRC1", "MSR1", "TGFB1"]

# ANTI_TUMOR full panel; missing-from-var-AND-marker_counts genes auto-dropped.
ANTI_TUMOR = [
    # MHC-II antigen presentation
    "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", "HLA-DPB1", "HLA-DRA", "HLA-DMA",
    "HLA-DMB", "HLA-DPA1", "HLA-DRB1", "CD74",
    # Pro-inflammatory cytokines
    "IL1A", "IL1B", "IL6", "TNF",
    # Th1 polarization
    "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
    # IRF / TF
    "IRF1", "IRF8", "STAT1",
    # M1 chemokines (CXCL9/10/11 to be recovered from marker_counts)
    "CSF2", "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def add_marker_counts_genes(a: sc.AnnData, genes: list[str]) -> sc.AnnData:
    """Append columns from obsm['marker_counts'] to a as new var entries.
    log1p-transformed, NaN imputed as 0 (treat "not measured" as below detection).
    Returns a NEW AnnData (in-memory) with X augmented.
    """
    mc = a.obsm.get("marker_counts", None)
    if mc is None:
        log("  no marker_counts — skip augment")
        return a
    add = [g for g in genes if g not in a.var_names and g in mc.columns]
    if not add:
        log("  no genes to add from marker_counts")
        return a
    log(f"  augmenting from marker_counts: {add}")

    # log1p(raw count), NaN -> 0
    extra = mc[add].copy()
    extra = extra.fillna(0.0).astype(np.float32)
    extra_log = np.log1p(extra.values)

    # Build augmented X: hstack lognorm X (need dense or csr) + extra_log
    X = a.X
    if sp.issparse(X):
        X = X.tocsr()
        # Convert extra_log dense (float32) to sparse (most are 0)
        extra_sp = sp.csr_matrix(extra_log)
        X_aug = sp.hstack([X, extra_sp], format="csr")
    else:
        X_aug = np.hstack([np.asarray(X), extra_log])

    var_aug = pd.concat([
        a.var,
        pd.DataFrame(index=add)
    ])
    a_aug = sc.AnnData(X=X_aug, obs=a.obs.copy(), var=var_aug,
                       obsm=dict(a.obsm), uns=dict(a.uns))
    log(f"  augmented shape: {a_aug.shape} (was {a.shape})")
    return a_aug


def main() -> None:
    log(f"loading {H5}")
    a = sc.read_h5ad(H5)
    log(f"  shape: {a.shape}")

    # Merge refined labels from 10c
    log(f"loading refined labels: {LBL}")
    lbl = pd.read_csv(LBL, index_col=0)
    a.obs["myeloid_subtype_refined"] = lbl.reindex(a.obs.index)["myeloid_subtype_refined"].astype(str).values
    a.obs["myeloid_major_type"] = a.obs["myeloid_subtype_refined"].map(MAJ_MAP).fillna("Other").astype(str).values
    n_lbl = (a.obs["myeloid_subtype_refined"] != "nan").sum()
    log(f"  refined labels mapped: {n_lbl}/{len(a.obs)}")

    # X is lognorm (per 10a output); confirm
    log(f"  X dtype={a.X.dtype}, sparse={sp.issparse(a.X)}, "
        f"max value (first 1k)={float(a.X[:1000].max() if not sp.issparse(a.X) else a.X[:1000].toarray().max()):.2f}")

    needed = sorted(set(M1_MARKERS + M2_MARKERS + ANTI_TUMOR))
    a_aug = add_marker_counts_genes(a, needed)

    log("M1/M2 polarization scoring (augmented)")
    for name, genes in [("M1", M1_MARKERS), ("M2", M2_MARKERS)]:
        present = [g for g in genes if g in a_aug.var_names]
        log(f"  {name}: {len(present)}/{len(genes)} present: {present}")
        if present:
            sc.tl.score_genes(a_aug, gene_list=present,
                              score_name=f"{name}_score_new",
                              random_state=0, use_raw=False)
        else:
            a_aug.obs[f"{name}_score_new"] = 0.0
    a_aug.obs["M1_M2_ratio_new"] = a_aug.obs["M1_score_new"] - a_aug.obs["M2_score_new"]

    # Diagnostic
    for col in ["M1_score_new", "M2_score_new"]:
        s = a_aug.obs[col]
        log(f"  {col}: mean={s.mean():.3f}, std={s.std():.3f}, "
            f"range=[{s.min():.3f}, {s.max():.3f}]")

    # myeloid_m1m2_scores.csv: groupby myeloid_subtype (raw subtype col)
    # Drop the OLD M1_score / M2_score / M1_M2_ratio columns from obs first to avoid duplicate-column rename collision
    drop_cols = [c for c in ["M1_score", "M2_score", "M1_M2_ratio"] if c in a_aug.obs.columns]
    obs = a_aug.obs.drop(columns=drop_cols).rename(columns={
        "M1_score_new": "M1_score",
        "M2_score_new": "M2_score",
        "M1_M2_ratio_new": "M1_M2_ratio",
    })
    sub_cols = ["M1_score", "M2_score", "M1_M2_ratio"]

    if "myeloid_subtype" in obs.columns:
        m1m2 = (obs.groupby("myeloid_subtype", observed=True)[sub_cols]
                .agg(["mean", "median", "std", "count"]))
        m1m2.columns = [f"{a_}_{b_}" for a_, b_ in m1m2.columns]
        m1m2.to_csv(FIG / "myeloid_m1m2_scores.csv")
        log(f"  myeloid_m1m2_scores.csv: {m1m2.shape}")

    if "myeloid_subtype_refined" in obs.columns:
        m1m2_r = (obs.groupby("myeloid_subtype_refined", observed=True)[sub_cols]
                  .agg(["mean", "median", "std", "count"]))
        m1m2_r.columns = [f"{a_}_{b_}" for a_, b_ in m1m2_r.columns]
        m1m2_r.to_csv(FIG / "myeloid_m1m2_scores_refined.csv")
        log(f"  myeloid_m1m2_scores_refined.csv: {m1m2_r.shape}")

    log("refresh panel_major_type_metadata.csv.gz")
    umap = a_aug.obsm["X_umap"]
    # obs may have duplicate "M1_score" columns after rename; pick first scalar series
    m1s = obs["M1_score"]
    if isinstance(m1s, pd.DataFrame): m1s = m1s.iloc[:, 0]
    m2s = obs["M2_score"]
    if isinstance(m2s, pd.DataFrame): m2s = m2s.iloc[:, 0]
    md = pd.DataFrame({
        "barcode": a_aug.obs.index,
        "dataset": a_aug.obs["dataset"].astype(str).values,
        "patient_id": a_aug.obs["patient_id"].astype(str).values,
        "tissue_type": a_aug.obs["tissue_type"].astype(str).values,
        "myeloid_subtype_refined": a_aug.obs["myeloid_subtype_refined"].astype(str).values,
        "myeloid_major_type": a_aug.obs["myeloid_major_type"].astype(str).values,
        "UMAP1": umap[:, 0], "UMAP2": umap[:, 1],
        "M1_score": np.asarray(m1s.values),
        "M2_score": np.asarray(m2s.values),
    })
    md.to_csv(FIG / "panel_major_type_metadata.csv.gz",
              index=False, compression="gzip")
    log(f"  panel_major_type_metadata.csv.gz: {len(md)} rows")

    log("Panel F: anti-tumor genes mean+pct (augmented)")
    present = [g for g in ANTI_TUMOR if g in a_aug.var_names]
    missing = [g for g in ANTI_TUMOR if g not in a_aug.var_names]
    log(f"  {len(present)}/{len(ANTI_TUMOR)} present: {present}")
    log(f"  missing: {missing}")

    # Pull matrix
    idx = [a_aug.var_names.get_loc(g) for g in present]
    X_sub = a_aug.X[:, idx]
    if sp.issparse(X_sub):
        X_sub = X_sub.toarray()
    df_expr = pd.DataFrame(X_sub, index=a_aug.obs.index, columns=present)
    df_expr["sub"] = a_aug.obs["myeloid_subtype_refined"].astype(str).values

    rows = []
    for st, sub_df in df_expr.groupby("sub", sort=False):
        if st in ("nan", ""):
            continue
        mat = sub_df[present].values
        for i, g in enumerate(present):
            rows.append({
                "subtype": st, "gene": g,
                "mean_log1p": float(mat[:, i].mean()),
                "pct_expressing": float((mat[:, i] > 0).mean()),
            })
    out = pd.DataFrame(rows)
    out.to_csv(FIG / "panel_F_antitumor_genes.csv", index=False)
    log(f"  panel_F_antitumor_genes.csv: {len(rows)} rows, "
        f"{len(present)} genes × {df_expr['sub'].nunique()} subtypes")

    log("DONE")


if __name__ == "__main__":
    main()
