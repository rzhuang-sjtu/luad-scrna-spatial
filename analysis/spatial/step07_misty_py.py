"""
Step 7: MISTy-equivalent multi-view random forest in Python (sklearn).

Rationale: R `mistyR` depends on `ridge`, which requires libgsl-dev (sudo not
available). Replicate the MISTy algorithm in pure Python:
  Targets:  PROGENy pathway scores (NF-κB, JAK-STAT, TNFa, EMT-like, Hypoxia, ...)
  Predictor views:
    - intraview: cell-type abundances at the same spot (c2l q05 abundance)
    - juxta:    mean abundances over k=6 nearest spatial neighbors
    - para:     mean abundances at radius r=2 hops (over the 6-NN graph)
  Model: RandomForestRegressor (ranger-equivalent), 100 trees, OOB R^2.
  Output: per-section importance matrices for each (target × view), then
          aggregated across the 12 sections (mean importance).

Inputs:
  ${DATA_ROOT}/ST/results/step05_progeny/cohort_with_progeny.h5ad
    obsm['q05_cell_abundance']   spots × 28 cell types
    obs['progeny_<pathway>']     pathway scores per spot
  ${DATA_ROOT}/ST/results/step01_qc/qc_summary.csv  for sample list

Outputs:
  ${DATA_ROOT}/ST/results/step07_misty/per_section_importance.csv
  ${DATA_ROOT}/ST/results/step07_misty/aggregated_importance.csv
  ${DATA_ROOT}/ST/results/step07_misty/per_section_r2.csv
  ${DATA_ROOT}/ST/results/step07_misty/heatmap_<view>.png
  ${DATA_ROOT}/ST/results/step07_misty/run.log
"""
from __future__ import annotations
import os, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

COHORT = Path("${DATA_ROOT}/ST/results/step05_progeny/cohort_with_progeny.h5ad")
OUT    = Path("${DATA_ROOT}/ST/results/step07_misty")
OUT.mkdir(exist_ok=True, parents=True)
LOG = OUT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

K = 6                # immediate neighbors
PATHWAYS = ["NFkB", "JAK-STAT", "TNFa", "TGFb", "Hypoxia", "MAPK", "EGFR", "p53", "VEGF"]
VIEWS = ["intra", "juxta", "para"]


def main():
    log(f"loading {COHORT}")
    a = sc.read_h5ad(str(COHORT))
    samples = sorted(a.obs["sample"].unique().tolist())
    log(f"   spots={a.n_obs}, samples={samples}")

    # cell-type abundance matrix (intraview features)
    abund = a.obsm["q05_cell_abundance"].copy()
    cell_types = list(abund.columns)
    log(f"   cell types: n={len(cell_types)}")

    # PROGENy targets
    targets = [f"progeny_{p}" for p in PATHWAYS if f"progeny_{p}" in a.obs.columns]
    log(f"   targets: {targets}")

    importance_records = []
    r2_records = []

    for s in samples:
        try:
            mask = (a.obs["sample"].values == s)
            sub = a[mask].copy()
            X_intra = abund.loc[sub.obs_names].to_numpy()  # spots × cell_types
            coords = sub.obsm["spatial"]                    # spots × 2

            # Build k-NN graph
            nn = NearestNeighbors(n_neighbors=K+1).fit(coords)
            dist, idx = nn.kneighbors(coords)               # idx[:,0] is self
            idx_juxta = idx[:, 1:K+1]                       # K nearest excluding self
            # para: 2-hop neighbors via the same K-NN graph
            # build neighbor mask, then take 2-hop union excluding self/juxta
            n = sub.n_obs
            juxta_set = [set(idx_juxta[i].tolist()) for i in range(n)]
            para_idx = []
            for i in range(n):
                hop2 = set()
                for j in juxta_set[i]:
                    hop2.update(juxta_set[j])
                hop2.discard(i)
                hop2 = hop2 - juxta_set[i]
                para_idx.append(np.array(sorted(hop2)))
            # Compute juxta and para feature matrices
            X_juxta = np.array([X_intra[idx_juxta[i]].mean(axis=0) for i in range(n)])
            X_para  = np.array([X_intra[p_i].mean(axis=0) if len(p_i) else np.zeros(len(cell_types))
                                for p_i in para_idx])

            view_X = {"intra": X_intra, "juxta": X_juxta, "para": X_para}

            for tgt in targets:
                y = sub.obs[tgt].to_numpy(dtype=float)
                if not np.isfinite(y).all():
                    y = np.nan_to_num(y, nan=np.nanmean(y))
                for vname, X in view_X.items():
                    if X is None or X.shape[0] != y.shape[0]:
                        continue
                    rf = RandomForestRegressor(n_estimators=100, n_jobs=-1, oob_score=True,
                                               random_state=0)
                    rf.fit(X, y)
                    fi = rf.feature_importances_
                    for ct, w in zip(cell_types, fi):
                        importance_records.append({"sample": s, "target": tgt.replace("progeny_",""),
                                                    "view": vname, "cell_type": ct, "importance": float(w)})
                    r2_records.append({"sample": s, "target": tgt.replace("progeny_",""),
                                       "view": vname, "oob_r2": float(rf.oob_score_)})
            log(f"   {s} done (spots={sub.n_obs})")
            del sub; gc.collect()
        except Exception as e:
            log(f"[ERROR] {s}: {type(e).__name__}: {e}\n{traceback.format_exc()}")

    imp_df = pd.DataFrame(importance_records)
    r2_df  = pd.DataFrame(r2_records)
    imp_df.to_csv(OUT / "per_section_importance.csv", index=False)
    r2_df.to_csv(OUT / "per_section_r2.csv", index=False)

    # Aggregate across sections: mean importance per (view, target, cell_type)
    agg = imp_df.groupby(["view", "target", "cell_type"])["importance"].mean().reset_index()
    agg.to_csv(OUT / "aggregated_importance.csv", index=False)

    # Heatmap per view: rows=cell_type, cols=target
    for v in VIEWS:
        sub = agg[agg["view"] == v].pivot(index="cell_type", columns="target", values="importance")
        if sub.empty:
            continue
        # Order rows by mean importance desc
        sub = sub.reindex(sub.mean(axis=1).sort_values(ascending=False).index)
        plt.figure(figsize=(max(6, 0.7*sub.shape[1]+3), max(6, 0.4*sub.shape[0]+1)))
        sns.heatmap(sub, cmap="rocket_r", annot=False, cbar_kws={"label": "RF importance"})
        plt.title(f"MISTy-py {v} view  (mean across sections)")
        plt.ylabel("cell type")
        plt.xlabel("PROGENy pathway")
        plt.tight_layout()
        plt.savefig(OUT / f"heatmap_{v}.png", dpi=140, bbox_inches="tight")
        plt.close()
    log(f"saved heatmaps for views: {VIEWS}")

    # Quick OOB R^2 summary
    if not r2_df.empty:
        per_target_view = r2_df.groupby(["target","view"])["oob_r2"].mean().unstack().round(3)
        per_target_view.to_csv(OUT / "oob_r2_summary.csv")
        log(f"OOB R^2 (per target × view, mean across sections):\n{per_target_view.to_string()}")
    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
