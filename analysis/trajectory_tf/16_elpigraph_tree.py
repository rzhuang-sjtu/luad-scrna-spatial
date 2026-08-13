"""Step 16: elastic principal tree on malignant UMAP (Monocle3-style skeleton).

Replaces unavailable Monocle3/R. Uses elpigraph-python to fit an elastic
principal tree on the malignant cells' UMAP coords, then exports edges
in (x_start, y_start, x_end, y_end) form for direct R plotting via
geom_segment.

Cells: 45,791 malignant LUAD cells with X_umap_mal + DPT pseudotime + dominant_MP
Tree: 60 nodes, fitted on 10,000 stratified-subsampled cells for speed,
      then projected back to all 45k cells for cell-level pseudotime mapping.

Outputs to ${WORK_ROOT}/luad_figures/fig2/:
  - monocle3_cells.csv.gz       barcode, UMAP1, UMAP2, monocle3_pseudotime, dominant_MP
  - monocle3_graph_segments.csv x_start, y_start, x_end, y_end, edge_id
  - monocle3_nodes.csv          node_id, x, y, n_cells_assigned
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
N_NODES = 60
N_SUBSAMPLE = 10000


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN, backed="r")
    log(f"  shape={a.shape}")

    # UMAP coords (use mal-specific UMAP)
    umap_key = "X_umap_mal" if "X_umap_mal" in a.obsm_keys() else "X_umap"
    umap = np.asarray(a.obsm[umap_key])
    log(f"  using {umap_key}, shape={umap.shape}")

    barcodes = a.obs.index.to_numpy()
    dominant_mp = a.obs["dominant_MP"].astype(str).to_numpy()

    # DPT pseudotime
    pt_df = pd.read_csv(PT_CSV)[["barcode","pseudotime"]]
    pt_lookup = dict(zip(pt_df["barcode"], pt_df["pseudotime"]))
    pseudotime_full = np.array([pt_lookup.get(bc, np.nan) for bc in barcodes])
    log(f"  pseudotime cells: {int(np.isfinite(pseudotime_full).sum())}/{len(pseudotime_full)}")

    # Stratified subsample by dominant_MP for tree fitting
    log(f"stratified subsample to {N_SUBSAMPLE} cells")
    rng = np.random.default_rng(0)
    idx_by_mp = {mp: np.where(dominant_mp == mp)[0]
                  for mp in np.unique(dominant_mp)}
    n_total = len(dominant_mp)
    keep_idx = []
    for mp, idx in idx_by_mp.items():
        n_take = max(50, int(N_SUBSAMPLE * len(idx) / n_total))
        if len(idx) <= n_take:
            keep_idx.extend(idx)
        else:
            keep_idx.extend(rng.choice(idx, size=n_take, replace=False))
    keep_idx = np.array(sorted(keep_idx))
    log(f"  subsampled to {len(keep_idx)} cells")

    umap_sub = umap[keep_idx].astype("float64")

    # Fit elastic principal tree
    log(f"fitting elastic principal tree (NumNodes={N_NODES})")
    import elpigraph
    t1 = time.time()
    res = elpigraph.computeElasticPrincipalTree(
        umap_sub, NumNodes=N_NODES, verbose=False,
        Lambda=0.01, Mu=0.1,
    )
    log(f"  tree fit done in {time.time()-t1:.1f}s")

    info = res[0] if isinstance(res, list) else res
    nodes = np.asarray(info["NodePositions"])
    edges = np.asarray(info["Edges"][0]).astype(int)  # (n_edges, 2)
    log(f"  nodes shape: {nodes.shape}; edges shape: {edges.shape}")

    # Build segments (start/end coords per edge)
    seg_rows = []
    for i, (a_i, b_i) in enumerate(edges):
        x0, y0 = nodes[a_i]
        x1, y1 = nodes[b_i]
        seg_rows.append({"edge_id": i, "node_start": int(a_i), "node_end": int(b_i),
                         "x_start": float(x0), "y_start": float(y0),
                         "x_end": float(x1), "y_end": float(y1)})
    seg_df = pd.DataFrame(seg_rows)
    seg_df.to_csv(FIG/"monocle3_graph_segments.csv", index=False)
    log(f"  monocle3_graph_segments.csv written ({len(seg_df)} edges)")

    # Assign each cell to its nearest tree node, propagate node-pseudotime to cells
    log("assigning each (full) cell to nearest tree node")
    from scipy.spatial import cKDTree
    tree = cKDTree(nodes)
    dist, node_assign = tree.query(umap.astype("float64"), k=1)
    log(f"  cell-to-node mean dist: {dist.mean():.3f}")

    # Compute per-node pseudotime as median of assigned cells' DPT pseudotime
    node_pt = np.full(len(nodes), np.nan)
    node_n_cells = np.zeros(len(nodes), dtype=int)
    for i in range(len(nodes)):
        mask = node_assign == i
        node_n_cells[i] = int(mask.sum())
        if mask.any() and np.isfinite(pseudotime_full[mask]).any():
            node_pt[i] = float(np.nanmedian(pseudotime_full[mask]))

    # Cell-level monocle3-style pseudotime = the node pseudotime they were assigned to
    monocle3_pt = node_pt[node_assign]

    nodes_df = pd.DataFrame({
        "node_id": [f"node_{i+1}" for i in range(len(nodes))],
        "x": nodes[:,0], "y": nodes[:,1],
        "n_cells_assigned": node_n_cells,
        "node_pseudotime": node_pt,
    })
    nodes_df.to_csv(FIG/"monocle3_nodes.csv", index=False)
    log(f"  monocle3_nodes.csv written ({len(nodes_df)} nodes)")

    # Per-cell output
    cells_df = pd.DataFrame({
        "barcode": barcodes,
        "UMAP1": umap[:,0],
        "UMAP2": umap[:,1],
        "monocle3_pseudotime": monocle3_pt,
        "dpt_pseudotime": pseudotime_full,
        "dominant_MP": dominant_mp,
    })
    cells_df.to_csv(FIG/"monocle3_cells.csv.gz",
                     index=False, compression="gzip")
    log(f"  monocle3_cells.csv.gz written ({len(cells_df)} cells)")

    # Sanity: correlation between monocle3 (node-aligned) pt and DPT
    from scipy.stats import spearmanr
    valid = np.isfinite(monocle3_pt) & np.isfinite(pseudotime_full)
    if valid.sum() > 100:
        rho, p = spearmanr(monocle3_pt[valid], pseudotime_full[valid])
        log(f"  Spearman(monocle3_pt, dpt_pseudotime) = {rho:.3f}, p={p:.2e}")
    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
