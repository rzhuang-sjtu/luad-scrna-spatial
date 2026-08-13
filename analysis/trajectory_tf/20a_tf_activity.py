"""Step 20a (Fig 3A): TF activity per MP via decoupler + CollecTRI.

Pseudo-bulk per MP (mean log1p expression across cells) → decoupler ULM
against CollecTRI prior network → TF activity score per MP.

Outputs:
  ${WORK_ROOT}/luad_figures/fig3/tf_activity_mp_matrix.csv  TF × MP
  ${WORK_ROOT}/luad_figures/fig3/tf_activity_top_per_mp.csv top TFs per MP
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
COLLECTRI = Path.home()/"luad/data/reference/collectri_symbols.tsv"
FIG3 = Path("${WORK_ROOT}/luad_figures/fig3")
TOPN = 25


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    log("loading CollecTRI")
    net = pd.read_csv(COLLECTRI, sep="\t")
    net = net[~net["source"].astype(str).str.contains("_")]
    net = net[net["source"].notna() & net["target"].notna()]
    net["weight"] = net["weight"].fillna(0).astype(float)
    net.loc[net["weight"] == 0, "weight"] = 1.0
    # Collapse duplicates: take signed-mean weight
    net = net.groupby(["source", "target"], as_index=False).agg(weight=("weight", "mean"))
    log(f"  CollecTRI net (deduped): {net.shape}; unique TFs: {net['source'].nunique()}")

    # Build pseudo-bulk per MP (mean expression across cells)
    log("building pseudo-bulks per dominant_MP")
    a.X = a.X.astype("float32")
    if a.X.min() < 0:
        log("  X has negatives → using layers['counts'] + log1p")
        a.X = np.log1p(a.layers["counts"].astype("float32"))
    pb = []
    mps = sorted([m for m in a.obs["dominant_MP"].astype(str).unique()
                  if m in {"MP1","MP2","MP3","MP4"}])
    for mp in mps:
        mask = (a.obs["dominant_MP"].astype(str) == mp).values
        X = a.X[mask]
        m = np.asarray(X.mean(axis=0)).ravel() if hasattr(X, "mean") \
            else X.mean(axis=0)
        pb.append(m)
    pb_mat = pd.DataFrame(np.array(pb), index=mps, columns=a.var_names)
    log(f"  pseudo-bulk: {pb_mat.shape}")

    log("running decoupler ULM (TF × MP)")
    import decoupler as dc
    log(f"  decoupler version: {dc.__version__}")
    log(f"  pb_mat type: {type(pb_mat)}, shape: {pb_mat.shape}")
    # decoupler 2.x: dc.mt.ulm returns (acts, pvals) tuple of DataFrames
    out = dc.mt.ulm(data=pb_mat, net=net, verbose=False)
    log(f"  ulm return type: {type(out)}")
    if isinstance(out, tuple):
        est, pval = out[0], out[1]
    elif hasattr(out, "obsm"):
        # AnnData-like return
        est = pd.DataFrame(out.obsm.get("score_ulm", out.X),
                            index=out.obs_names, columns=out.var_names)
        pval = None
    else:
        est = out
        pval = None
    log(f"  TF activity shape: {est.shape}")

    # Determine orientation: rows should be MP samples; columns are TFs
    if est.shape[0] == len(mps):
        tf_by_mp = est.T
    else:
        tf_by_mp = est
    tf_by_mp.to_csv(FIG3/"tf_activity_mp_matrix.csv")
    log(f"  tf_activity_mp_matrix.csv saved")

    # Top 25 TFs per MP (by absolute activity)
    rows = []
    for mp in mps:
        col = tf_by_mp[mp]
        top_pos = col.sort_values(ascending=False).head(TOPN)
        top_neg = col.sort_values(ascending=True).head(TOPN)
        for tf, score in top_pos.items():
            rows.append({"MP": mp, "TF": tf, "activity": float(score),
                         "direction": "Up", "rank": top_pos.index.get_loc(tf)+1})
        for tf, score in top_neg.items():
            rows.append({"MP": mp, "TF": tf, "activity": float(score),
                         "direction": "Down", "rank": top_neg.index.get_loc(tf)+1})
    top_df = pd.DataFrame(rows)
    top_df.to_csv(FIG3/"tf_activity_top_per_mp.csv", index=False)
    log(f"  tf_activity_top_per_mp.csv ({len(top_df)} rows)")
    for mp in mps:
        sub = top_df[(top_df["MP"]==mp) & (top_df["direction"]=="Up")].head(8)
        log(f"\n  Top-8 UP TFs in {mp}:")
        log(sub[["TF","activity"]].round(3).to_string(index=False))

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
