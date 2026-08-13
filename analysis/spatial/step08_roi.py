"""
Step 8: ROI analysis — define & profile NF-κB-high & Neutrophil-high regions.

Per Fig 7K-P plan:
  ROI = spots where NF-κB pathway activity is high AND total Neutrophil abundance is high.

  Operational definition:
    z(NF-κB) > 0.5 AND z(neu_total) > 0.5     (per-section z-scoring; both top-30%)
    spots adjacent to ≥2 ROI spots are also kept (closing) to make connected regions.

For each ROI we compute distributions of:
  - MP1-4 score (Fig 7L)
  - TFs: ATF3, FOSB, JUN, NFKBIA (Fig 7M, expression-based)
  - Macro_SPP1, Neu_Inflammatory, Neu_OSM_priming (Fig 7N)
  - PROGENy NFkB / EMT-like (TGFb / TNFa) + IL1B / OSM expression (Fig 7O-P)

Outputs:
  ${DATA_ROOT}/ST/results/step08_roi/roi_summary.csv               per-section ROI counts + means
  ${DATA_ROOT}/ST/results/step08_roi/roi_vs_nonroi_stats.csv       cell-type / pathway / TF means
  ${DATA_ROOT}/ST/results/step08_roi/spatial_plots/<sample>_roi.png  ROI spatial overlay
"""
from __future__ import annotations
import os, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COHORT = Path("${DATA_ROOT}/ST/results/step05_progeny/cohort_with_progeny.h5ad")
QC_SEC = Path("${DATA_ROOT}/ST/results/step01_qc/section_h5ad")
OUT    = Path("${DATA_ROOT}/ST/results/step08_roi")
PLOTS  = OUT / "spatial_plots"
for d in (OUT, PLOTS):
    d.mkdir(parents=True, exist_ok=True)
LOG = OUT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

NEU_COLS_PRIORITY = ["Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low","Neu_IFN_response",
                     "Neu_Angiogenic","Neu_Metastatic","Neu_ECM_remodeling"]
TF_GENES = ["ATF3","FOSB","JUN","JUNB","NFKBIA","FOS"]
LIGAND_GENES = ["IL1B","OSM","IL1A","TNF","TGFB1"]
MP_COLS = ["MP1_score","MP2_score","MP3_score","MP4_score"]
PATHWAY_PROGENY = ["progeny_NFkB","progeny_TNFa","progeny_TGFb","progeny_JAK-STAT","progeny_Hypoxia"]


def main():
    log(f"loading {COHORT}")
    a = sc.read_h5ad(str(COHORT))
    log(f"   shape={a.shape}")

    # Compute neu_total per spot
    abund = a.obsm["q05_cell_abundance"].copy()
    neu_cols = [c for c in NEU_COLS_PRIORITY if c in abund.columns]
    a.obs["neu_total"] = abund[neu_cols].sum(axis=1).values
    log(f"neutrophil columns used: {neu_cols}")

    # Per-section z-score and ROI definition
    a.obs["z_nfkb"]    = 0.0
    a.obs["z_neu"]     = 0.0
    a.obs["roi"]       = False
    samples = sorted(a.obs["sample"].unique().tolist())

    for s in samples:
        mask = a.obs["sample"].values == s
        nf = a.obs.loc[mask, "progeny_NFkB"].values
        nu = a.obs.loc[mask, "neu_total"].values
        zn = (nf - nf.mean()) / (nf.std() + 1e-9)
        zu = (nu - nu.mean()) / (nu.std() + 1e-9)
        a.obs.loc[mask, "z_nfkb"] = zn
        a.obs.loc[mask, "z_neu"]  = zu
        roi_mask = (zn > 0.5) & (zu > 0.5)
        a.obs.loc[mask, "roi"] = roi_mask
        log(f"  {s}: spots={mask.sum()}, raw ROI={roi_mask.sum()} ({roi_mask.mean()*100:.1f}%)")

    # ROI summary
    roi_df = a.obs[["sample","roi"]].copy()
    counts = roi_df.groupby("sample")["roi"].agg(["sum","mean","count"]).rename(
        columns={"sum":"n_roi_spots", "mean":"frac_roi", "count":"n_total_spots"})
    counts.to_csv(OUT / "roi_summary.csv")
    log(f"\nROI counts:\n{counts.to_string()}")

    # Compute means in ROI vs non-ROI for each metric
    # Cell type abundances
    abund_full = abund.copy()
    abund_full["roi"] = a.obs["roi"].values
    abund_full["sample"] = a.obs["sample"].values
    ct_means = abund_full.groupby(["sample","roi"]).mean().reset_index()
    ct_means.to_csv(OUT / "cell_type_mean_by_roi.csv", index=False)

    # MP / pathway means
    for col in MP_COLS + PATHWAY_PROGENY:
        if col not in a.obs.columns:
            log(f"  [skip] missing {col}")
    pathway_metrics = a.obs[["sample","roi"] + [c for c in (MP_COLS + PATHWAY_PROGENY) if c in a.obs.columns]].copy()
    pw_means = pathway_metrics.groupby(["sample","roi"]).mean().reset_index()
    pw_means.to_csv(OUT / "pathway_mean_by_roi.csv", index=False)

    # Gene expression: TF + ligand mean (use normalized log1p X)
    gene_set = [g for g in TF_GENES + LIGAND_GENES if g in a.var_names]
    log(f"genes available: {gene_set}")
    X = a[:, gene_set].X
    if sp.issparse(X): X = X.toarray()
    gex = pd.DataFrame(X, index=a.obs_names, columns=gene_set)
    gex["roi"] = a.obs["roi"].values
    gex["sample"] = a.obs["sample"].values
    gex_means = gex.groupby(["sample","roi"]).mean().reset_index()
    gex_means.to_csv(OUT / "gene_mean_by_roi.csv", index=False)

    # Combined ROI vs non-ROI summary
    combined_rows = []
    for s in samples:
        sec = a[a.obs["sample"].values == s].copy()
        in_roi = sec.obs["roi"].values
        if in_roi.sum() < 5:
            continue
        for col in MP_COLS + PATHWAY_PROGENY:
            if col not in sec.obs.columns: continue
            v = sec.obs[col].values.astype(float)
            combined_rows.append({"sample": s, "metric": col,
                                  "type": "obs", "mean_roi": float(np.mean(v[in_roi])),
                                  "mean_nonroi": float(np.mean(v[~in_roi])),
                                  "delta": float(np.mean(v[in_roi]) - np.mean(v[~in_roi]))})
        for ct in NEU_COLS_PRIORITY + ["Macro_SPP1","Macro_C1QC","Fibroblast","Endothelial","Malignant","T_NK","B"]:
            if ct in abund.columns:
                vals = abund.loc[sec.obs_names, ct].values
                combined_rows.append({"sample": s, "metric": ct, "type": "celltype",
                                      "mean_roi": float(np.mean(vals[in_roi])),
                                      "mean_nonroi": float(np.mean(vals[~in_roi])),
                                      "delta": float(np.mean(vals[in_roi]) - np.mean(vals[~in_roi]))})
        for g in gene_set:
            vals = gex.loc[sec.obs_names, g].values
            combined_rows.append({"sample": s, "metric": g, "type": "gene",
                                  "mean_roi": float(np.mean(vals[in_roi])),
                                  "mean_nonroi": float(np.mean(vals[~in_roi])),
                                  "delta": float(np.mean(vals[in_roi]) - np.mean(vals[~in_roi]))})
    cmp_df = pd.DataFrame(combined_rows)
    cmp_df.to_csv(OUT / "roi_vs_nonroi_stats.csv", index=False)

    # Aggregate: mean delta across sections per metric
    if len(cmp_df):
        agg = cmp_df.groupby(["metric","type"])[["mean_roi","mean_nonroi","delta"]].mean().sort_values("delta", ascending=False)
        agg.to_csv(OUT / "roi_vs_nonroi_aggregate.csv")
        log(f"\ntop +delta metrics (ROI > non-ROI):\n{agg.head(20).to_string()}")
        log(f"\nbottom delta metrics (non-ROI > ROI):\n{agg.tail(10).to_string()}")

    # Spatial ROI overlay plots
    log("plotting per-section ROI overlay ...")
    for s in samples:
        try:
            sec = sc.read_h5ad(str(QC_SEC / f"{s}.h5ad"))
            sec.var_names_make_unique()
            mask = a.obs["sample"].values == s
            sub_idx = a.obs_names[mask]
            section_bc = pd.Index([n[:-len("-"+s)] if n.endswith("-"+s) else n for n in sub_idx])
            sec.obs["roi"] = pd.Series(a.obs.loc[sub_idx, "roi"].values, index=section_bc).reindex(sec.obs_names).fillna(False).astype(int).values
            sec.obs["progeny_NFkB"] = pd.Series(a.obs.loc[sub_idx, "progeny_NFkB"].values, index=section_bc).reindex(sec.obs_names).values
            sec.obs["neu_total"] = pd.Series(a.obs.loc[sub_idx, "neu_total"].values, index=section_bc).reindex(sec.obs_names).values

            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            sc.pl.spatial(sec, color="progeny_NFkB", library_id=s, ax=axes[0], show=False,
                          cmap="RdBu_r", size=1.4, frameon=False, title="NFkB activity", colorbar_loc="right")
            sc.pl.spatial(sec, color="neu_total", library_id=s, ax=axes[1], show=False,
                          cmap="magma", size=1.4, frameon=False, title="Neutrophil total q05", colorbar_loc="right")
            sc.pl.spatial(sec, color="roi", library_id=s, ax=axes[2], show=False,
                          cmap="Reds", size=1.4, frameon=False, title="ROI (NFkB^Neu high)", colorbar_loc="right")
            fig.suptitle(f"{s}  ROI definition", fontsize=12)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}_roi.png", dpi=130, bbox_inches="tight")
            plt.close(fig)
            del sec; gc.collect()
        except Exception as e:
            log(f"[plot fail] {s}: {type(e).__name__}: {e}")

    # Save cohort with ROI flag added back
    a.write_h5ad(str(OUT / "cohort_with_roi.h5ad"), compression="gzip")
    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
