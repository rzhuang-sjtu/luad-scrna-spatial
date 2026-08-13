"""Step 15: PAGA + force-directed layout trajectory for Fig 2G.

Pipeline:
  1. Load malignant h5ad (45k cells); rebuild neighbors on X_pca_mal_harmony.
  2. sc.tl.paga(adata, groups='dominant_MP') — coarse PAGA between MP states.
  3. sc.tl.draw_graph(adata, init_pos='X_umap') — ForceAtlas2 / FR force-layout.
     Falls back to igraph 'fr' if fa2 backend missing.
  4. Bin pseudotime (60 bins) → median FA1/FA2 per bin → ordered skeleton path.

Outputs:
  ${WORK_ROOT}/luad_figures/fig2/trajectory_cells.csv.gz
    barcode, FA1, FA2, pseudotime, pseudotime_winsorized, dominant_MP
  ${WORK_ROOT}/luad_figures/fig2/trajectory_graph.csv
    node_id, x, y, pseudotime, n_cells
  ${WORK_ROOT}/luad_figures/fig2/paga_connectivities.csv
    PAGA group×group connectivity matrix
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
PT_CSV = Path.home()/"luad/results/step10b_pseudotime.csv.gz"
FIG = Path("${WORK_ROOT}/luad_figures/fig2")
N_BINS = 60


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    log("rebuilding neighbors on X_pca_mal_harmony")
    sc.pp.neighbors(a, use_rep="X_pca_mal_harmony", n_neighbors=30, random_state=0)

    log("PAGA on dominant_MP")
    sc.tl.paga(a, groups="dominant_MP")
    paga_conn = pd.DataFrame(
        a.uns["paga"]["connectivities"].toarray(),
        index=a.obs["dominant_MP"].cat.categories
              if hasattr(a.obs["dominant_MP"], "cat") else
              sorted(a.obs["dominant_MP"].unique()),
        columns=a.obs["dominant_MP"].cat.categories
              if hasattr(a.obs["dominant_MP"], "cat") else
              sorted(a.obs["dominant_MP"].unique()),
    )
    paga_conn.to_csv(FIG/"paga_connectivities.csv")
    log("PAGA connectivity:")
    log(paga_conn.round(3).to_string())

    log("running sc.tl.draw_graph (force-directed)")
    t1 = time.time()
    try:
        # Default: uses fa2 if available, else fr via igraph
        sc.tl.draw_graph(a, init_pos="X_umap", random_state=0)
        layout_used = "fa"
    except Exception as e:
        log(f"  default draw_graph failed: {e}; falling back to fr (igraph)")
        sc.tl.draw_graph(a, init_pos="X_umap", random_state=0, layout="fr")
        layout_used = "fr"
    log(f"  draw_graph done in {time.time()-t1:.1f}s (layout={layout_used})")

    fa = a.obsm["X_draw_graph_fa"] if "X_draw_graph_fa" in a.obsm else \
         a.obsm[next(k for k in a.obsm if k.startswith("X_draw_graph"))]
    log(f"  FA layout shape: {fa.shape}")

    # Get pseudotime from saved CSV
    log("merging pseudotime")
    pt = pd.read_csv(PT_CSV)[["barcode","pseudotime"]]
    pt_lookup = dict(zip(pt["barcode"], pt["pseudotime"]))
    pseudotime = np.array([pt_lookup.get(bc, np.nan) for bc in a.obs.index])

    q99 = float(np.nanquantile(pseudotime, 0.99))
    pseudotime_w = np.clip(pseudotime, None, q99)
    log(f"  q99={q99:.4f}; cells with pt: {int(np.isfinite(pseudotime).sum())}")

    # Per-cell output
    cells_df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "FA1": fa[:, 0],
        "FA2": fa[:, 1],
        "pseudotime": pseudotime,
        "pseudotime_winsorized": pseudotime_w,
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
    })
    cells_df.to_csv(FIG/"trajectory_cells.csv.gz",
                     index=False, compression="gzip")
    log(f"  trajectory_cells.csv.gz written ({len(cells_df)} rows)")

    # Skeleton: bin pseudotime, take median FA per bin, order by bin
    log(f"building skeleton ({N_BINS} pseudotime bins)")
    valid = cells_df.dropna(subset=["pseudotime"]).copy()
    valid["pt_bin"] = pd.qcut(valid["pseudotime"], q=N_BINS,
                                labels=False, duplicates="drop")
    skel = (valid.groupby("pt_bin")
                  .agg(x=("FA1", "median"),
                       y=("FA2", "median"),
                       pseudotime=("pseudotime", "median"),
                       n_cells=("barcode", "size"))
                  .reset_index())
    skel = skel.sort_values("pt_bin").reset_index(drop=True)
    skel.insert(0, "node_id", [f"node_{i+1}" for i in range(len(skel))])
    skel = skel[["node_id", "x", "y", "pseudotime", "n_cells"]]
    skel.to_csv(FIG/"trajectory_graph.csv", index=False)
    log(f"  trajectory_graph.csv written ({len(skel)} nodes)")
    log("first/last 5 skeleton nodes:")
    log(skel.head(5).round(4).to_string(index=False))
    log("...")
    log(skel.tail(5).round(4).to_string(index=False))

    log(f"\nDONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
