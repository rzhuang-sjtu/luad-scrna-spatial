"""
Step 4: MP1-4 spot-level scoring using LUAD MP signatures (top-100 genes per MP).

Inputs:
  ${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad   28k spots, 9017 genes (raw counts in X)
  ~/luad/results/step6_mp_signatures_top100.csv                        MP × 100 genes

Outputs:
  ${DATA_ROOT}/ST/results/step04_mp_scoring/cohort_with_mp.h5ad         AnnData with normalized X + MP scores in obs
  ${DATA_ROOT}/ST/results/step04_mp_scoring/spot_mp_scores.csv          per-spot (sample, barcode, MP1-5)
  ${DATA_ROOT}/ST/results/step04_mp_scoring/mp_signatures_used.csv      gene lists effectively used
  ${DATA_ROOT}/ST/results/step04_mp_scoring/spatial_plots/<sample>.png  4-panel MP1-4 spatial maps
"""
from __future__ import annotations
import os, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COHORT = Path("${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad")
SIGS   = Path(os.path.expanduser("~/luad/results/step6_mp_signatures_top100.csv"))
QC_SEC = Path("${DATA_ROOT}/ST/results/step01_qc/section_h5ad")
OUT    = Path("${DATA_ROOT}/ST/results/step04_mp_scoring")
PLOTS  = OUT / "spatial_plots"
OUT.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

LOG = OUT / "run.log"
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def main():
    log(f"loading {COHORT}")
    a = sc.read_h5ad(str(COHORT))
    log(f"   shape={a.shape}; X={type(a.X).__name__} dtype={a.X.dtype}")

    # Save raw counts in layer
    a.layers["counts"] = a.X.copy()
    # Normalize + log
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    log("normalized log1p CPM")

    # Read signatures
    sig_df = pd.read_csv(SIGS)
    log(f"signatures: {sig_df.shape}; MPs: {sig_df['MP'].unique().tolist()}")

    used_records = []
    mp_score_cols = []
    for mp in sorted(sig_df["MP"].unique()):
        genes = sig_df.loc[sig_df["MP"] == mp, "gene"].tolist()
        present = [g for g in genes if g in a.var_names]
        used_records.append({"MP": mp, "n_total": len(genes), "n_present": len(present)})
        col = f"{mp}_score"
        sc.tl.score_genes(a, gene_list=present, score_name=col, ctrl_size=200, random_state=0)
        mp_score_cols.append(col)
        log(f"   {mp}: {len(present)}/{len(genes)} genes present -> obs['{col}']")

    # dominant MP per spot (within MP1-MP4 to match Fig 7 plan)
    main_cols = [f"MP{i}_score" for i in (1,2,3,4)]
    a.obs["dominant_MP_4"] = a.obs[main_cols].idxmax(axis=1).str.replace("_score", "")
    a.obs["dominant_MP_4_score"] = a.obs[main_cols].max(axis=1)

    pd.DataFrame(used_records).to_csv(OUT / "mp_signatures_used.csv", index=False)

    # Save per-spot CSV
    per_spot = a.obs[["sample"] + mp_score_cols + ["dominant_MP_4", "dominant_MP_4_score"]].copy()
    per_spot.index.name = "spot"
    per_spot.to_csv(OUT / "spot_mp_scores.csv")
    log(f"per-spot scores -> spot_mp_scores.csv  ({per_spot.shape})")

    # Save cohort with MP scores
    out_h5 = OUT / "cohort_with_mp.h5ad"
    a.write_h5ad(str(out_h5), compression="gzip")
    log(f"cohort with MP scores -> {out_h5}  ({out_h5.stat().st_size/1e9:.2f} GB)")

    # ---- spatial plots: per-section MP1-4 ----
    log("plotting per-section MP1-4 spatial maps ...")
    samples = sorted(a.obs["sample"].unique().tolist())
    for s in samples:
        try:
            sec = sc.read_h5ad(str(QC_SEC / f"{s}.h5ad"))
            sec.var_names_make_unique()
            mask = a.obs["sample"].values == s
            sub_idx = a.obs_names[mask]
            section_bc = pd.Index([n[:-len("-"+s)] if n.endswith("-"+s) else n for n in sub_idx])
            for col in main_cols:
                sec.obs[col] = pd.Series(a.obs.loc[sub_idx, col].values, index=section_bc).reindex(sec.obs_names).values
            fig, axes = plt.subplots(2, 2, figsize=(10, 9))
            axes = axes.flatten()
            for i, col in enumerate(main_cols):
                sc.pl.spatial(sec, color=col, library_id=s, ax=axes[i], show=False,
                              cmap="RdBu_r", size=1.4, frameon=False, title=col, colorbar_loc="right")
            fig.suptitle(f"{s}  MP1-4 spot scores", fontsize=12)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}.png", dpi=130, bbox_inches="tight")
            plt.close(fig)
            del sec; gc.collect()
        except Exception as e:
            log(f"[plot fail] {s}: {type(e).__name__}: {e}")

    # Per-sample mean MP scores
    means = a.obs.groupby("sample")[mp_score_cols].mean().round(3)
    means.to_csv(OUT / "per_sample_mean_mp.csv")
    log(f"per-sample mean MP -> per_sample_mean_mp.csv\n{means.to_string()}")

    # Per-sample dominant_MP_4 distribution
    dist = a.obs.groupby("sample")["dominant_MP_4"].value_counts().unstack(fill_value=0)
    dist.to_csv(OUT / "per_sample_dominant_mp4.csv")
    log(f"per-sample dominant_MP_4 dist:\n{dist.to_string()}")

    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
