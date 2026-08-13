"""
Step 5: PROGENy 14-pathway scoring on spots via decoupler-py.

Inputs:
  ${DATA_ROOT}/ST/results/step04_mp_scoring/cohort_with_mp.h5ad

Outputs:
  ${DATA_ROOT}/ST/results/step05_progeny/cohort_with_progeny.h5ad
  ${DATA_ROOT}/ST/results/step05_progeny/spot_progeny_scores.csv
  ${DATA_ROOT}/ST/results/step05_progeny/per_sample_mean_progeny.csv
  ${DATA_ROOT}/ST/results/step05_progeny/spatial_plots/<sample>.png  (NF-κB, JAK-STAT, TNFa, EMT-related)
"""
from __future__ import annotations
import os, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_H5  = Path("${DATA_ROOT}/ST/results/step04_mp_scoring/cohort_with_mp.h5ad")
QC_SEC = Path("${DATA_ROOT}/ST/results/step01_qc/section_h5ad")
OUT    = Path("${DATA_ROOT}/ST/results/step05_progeny")
PLOTS  = OUT / "spatial_plots"
OUT.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)
LOG = OUT / "run.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

# Highlighted pathways for spatial plots
HIGHLIGHT = ["NFkB", "JAK-STAT", "TNFa", "TGFb", "Hypoxia", "MAPK"]


def main():
    log(f"loading {IN_H5}")
    a = sc.read_h5ad(str(IN_H5))
    log(f"   shape={a.shape}")

    import decoupler as dc
    log(f"decoupler version: {dc.__version__}")

    # Load PROGENy collectri-like network: top 500 genes per pathway
    log("loading PROGENy net (top=500) ...")
    progeny = dc.op.progeny(organism="human", top=500)
    log(f"   net shape={progeny.shape}, pathways={sorted(progeny['source'].unique())}")

    # Run MLM (multivariate linear model) — standard PROGENy method
    log("dc.mt.mlm on log-normalized cohort ...")
    res = dc.mt.mlm(data=a, net=progeny, verbose=True)
    # decoupler 2.x: returns (estimate, pvals) when given AnnData OR writes to obsm directly
    if isinstance(res, tuple) and len(res) == 2:
        estimate, pvals = res
    else:
        # if writes to obsm
        estimate = a.obsm["score_mlm"] if "score_mlm" in a.obsm else None
        pvals    = a.obsm["padj_mlm"]  if "padj_mlm"  in a.obsm else None
        if estimate is None:
            raise RuntimeError(f"unexpected dc.mt.mlm return: {type(res)}")
    if not isinstance(estimate, pd.DataFrame):
        estimate = pd.DataFrame(estimate, index=a.obs_names)
    a.obsm["progeny_mlm"] = estimate
    if pvals is not None:
        a.obsm["progeny_mlm_pvals"] = pvals

    # Add main pathway scores to obs for plotting
    for col in estimate.columns:
        a.obs[f"progeny_{col}"] = estimate[col].values

    # Save
    out_h5 = OUT / "cohort_with_progeny.h5ad"
    a.write_h5ad(str(out_h5), compression="gzip")
    log(f"saved cohort -> {out_h5}  ({out_h5.stat().st_size/1e9:.2f} GB)")

    # Save tabular summary
    estimate.to_csv(OUT / "spot_progeny_scores.csv")
    means = estimate.copy()
    means["sample"] = a.obs["sample"].values
    per_sample = means.groupby("sample").mean().round(3)
    per_sample.to_csv(OUT / "per_sample_mean_progeny.csv")
    log(f"per-sample mean PROGENy:\n{per_sample.to_string()}")

    # Per-section spatial plots
    log("plotting per-section PROGENy spatial maps ...")
    samples = sorted(a.obs["sample"].unique().tolist())
    avail = [p for p in HIGHLIGHT if p in estimate.columns]
    for s in samples:
        try:
            sec = sc.read_h5ad(str(QC_SEC / f"{s}.h5ad"))
            sec.var_names_make_unique()
            mask = a.obs["sample"].values == s
            sub_idx = a.obs_names[mask]
            section_bc = pd.Index([n[:-len("-"+s)] if n.endswith("-"+s) else n for n in sub_idx])
            for col in avail:
                key = f"progeny_{col}"
                sec.obs[key] = pd.Series(a.obs.loc[sub_idx, key].values, index=section_bc).reindex(sec.obs_names).values
            ncols = 3
            nrows = (len(avail) + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
            axes = np.array(axes).flatten()
            for i, col in enumerate(avail):
                key = f"progeny_{col}"
                sc.pl.spatial(sec, color=key, library_id=s, ax=axes[i], show=False,
                              cmap="RdBu_r", size=1.4, frameon=False, title=col, colorbar_loc="right")
            for j in range(len(avail), len(axes)):
                axes[j].axis("off")
            fig.suptitle(f"{s}  PROGENy pathway activity (MLM)", fontsize=12)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}.png", dpi=130, bbox_inches="tight")
            plt.close(fig)
            del sec; gc.collect()
        except Exception as e:
            log(f"[plot fail] {s}: {type(e).__name__}: {e}")

    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
