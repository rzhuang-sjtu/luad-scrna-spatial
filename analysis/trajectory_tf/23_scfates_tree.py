"""Step 23: scFates principal-tree trajectory (Monocle-style).

Replaces previous Fig 2G / Fig 3B trajectory data with a proper tree-shaped
output produced by scFates SimplePPT in the diffusion-map-reduced space.

Pipeline:
  1. Load malignant h5ad (45k cells with X_pca_mal_harmony).
  2. Compute neighbors + diffmap (5 components) on Harmony space.
  3. scf.tl.curve(...) OR scf.tl.tree(..., method="ppt") to fit principal tree.
  4. Set root from MP4-dominant cells (highest MP4 score).
  5. scf.tl.pseudotime() → cell-level pseudotime + branch ID.
  6. Export:
     - tree_cells.csv.gz   : barcode, Component_1/2 (diffmap DC1/DC2),
                              UMAP1/UMAP2, pseudotime, branch_id, dominant_MP
     - tree_nodes.csv      : node_id, x, y, parent_id, n_cells_assigned
     - tree_edges.csv      : edge_id, src, tgt, x_src, y_src, x_tgt, y_tgt

All saved to ${WORK_ROOT}/luad_figures/fig2/ (overwrites stale files).
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
FIG2 = Path("${WORK_ROOT}/luad_figures/fig2")
N_NODES = 50


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    # Rebuild neighbors on mal-specific Harmony PCs
    log("rebuilding neighbors on X_pca_mal_harmony")
    sc.pp.neighbors(a, use_rep="X_pca_mal_harmony", n_neighbors=30, random_state=0)

    log("computing diffmap (5 components)")
    sc.tl.diffmap(a, n_comps=5, random_state=0)
    # Skip first diffusion component (constant); use DC1, DC2 as 2D embedding
    a.obsm["X_dm"] = a.obsm["X_diffmap"][:, 1:]
    log(f"  X_dm shape: {a.obsm['X_dm'].shape}")

    log(f"fitting scFates ElPiGraph tree (Nodes={N_NODES})")
    import scFates as scf
    # Use 4D diffmap (DC1-DC4) so tree can branch in higher-dim;
    # later project to DC1/DC2 for 2D visualization
    a.obsm["X_emb"] = a.obsm["X_dm"][:, :4].astype(np.float64)
    t1 = time.time()
    scf.tl.tree(a, Nodes=N_NODES, method="epg", use_rep="X_emb",
                 epg_lambda=0.01, epg_mu=0.1, seed=0)
    log(f"  tree fit done in {time.time()-t1:.1f}s")
    log(f"  uns['graph'] keys: {list(a.uns.get('graph', {}).keys())}")

    # Set root: tree NODE closest to MP4-max cell (in X_emb 4D space)
    mp4 = a.obs["MP4_score"].values
    icell = int(np.argmax(mp4))
    log(f"  MP4-max cell idx={icell} ({a.obs.index[icell]}), MP4={mp4[icell]:.3f}")
    nodes_F = np.asarray(a.uns["graph"]["F"])
    n_nodes_actual = max(nodes_F.shape)
    if nodes_F.shape[0] != n_nodes_actual:
        nodes_F = nodes_F.T
    log(f"  tree nodes: {nodes_F.shape}")
    cell_pos = a.obsm["X_emb"][icell]
    from scipy.spatial.distance import cdist
    dists = cdist(cell_pos.reshape(1, -1), nodes_F).ravel()
    iroot_node = int(np.argmin(dists))
    log(f"  iroot_node={iroot_node} (closest tree node to MP4-max cell)")
    scf.tl.root(a, root=iroot_node)

    log("computing pseudotime")
    try:
        scf.tl.pseudotime(a, n_jobs=8, n_map=1, seed=0)
    except (IndexError, KeyError) as e:
        log(f"  pseudotime auto-color failed ({e}); proceeding manually")
        # Pre-populate milestones_colors to bypass scFates color setup bug
        if "milestones" in a.obs.columns:
            n_ms = a.obs["milestones"].astype(str).nunique()
            a.uns["milestones_colors"] = ["#"+f"{i:06x}" for i in
                                            np.linspace(0, 0xffffff, max(n_ms, 1)).astype(int)]
            scf.tl.pseudotime(a, n_jobs=8, n_map=1, seed=0)

    # Inspect graph structure
    g = a.uns["graph"]
    nodes = g["F"]                       # node positions in X_emb space (n_nodes × 2)
    edges = np.asarray(g["B"]).astype(int)  # adjacency (n_nodes × n_nodes)
    nodes = np.asarray(nodes).T if nodes.shape[0] != N_NODES else nodes
    if nodes.shape[1] != 2 and nodes.shape[0] == 2:
        nodes = nodes.T
    log(f"  tree nodes: {nodes.shape}; edge matrix: {edges.shape}")
    n_edges_total = int(edges.sum() // 2) if edges.ndim == 2 else len(edges)
    log(f"  n_edges (undirected): {n_edges_total}")

    # Build edge list from adjacency matrix; project node coords to first 2 dims
    nodes_2d = nodes[:, :2]
    edge_list = []
    if edges.ndim == 2:
        seen = set()
        for i in range(edges.shape[0]):
            for j in range(edges.shape[1]):
                if edges[i, j] > 0 and i < j and (i, j) not in seen:
                    seen.add((i, j))
                    edge_list.append({
                        "edge_id": len(edge_list),
                        "src": i, "tgt": j,
                        "x_src": float(nodes_2d[i, 0]), "y_src": float(nodes_2d[i, 1]),
                        "x_tgt": float(nodes_2d[j, 0]), "y_tgt": float(nodes_2d[j, 1]),
                    })

    edges_df = pd.DataFrame(edge_list)
    edges_df.to_csv(FIG2/"tree_edges.csv", index=False)
    log(f"  tree_edges.csv: {edges_df.shape}")

    # Per-cell tree segment assignment (which node is closest)
    seg_assign = a.obs.get("seg") if "seg" in a.obs else None
    milestones = a.obs.get("milestones") if "milestones" in a.obs else None

    # Cell metadata
    pt = a.obs["t"].values if "t" in a.obs else a.obs.get("pseudotime", np.full(a.n_obs, np.nan)).values
    log(f"  pseudotime range: [{np.nanmin(pt):.4f}, {np.nanmax(pt):.4f}]")
    pt_rank = pd.Series(pt).rank(pct=True).values

    # For 2D visualization: Component_1/2 = first 2 diffmap components
    # (the X_emb is 4D used for tree fitting only)
    cells_df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "Component_1": a.obsm["X_dm"][:, 0],
        "Component_2": a.obsm["X_dm"][:, 1],
        "UMAP1": a.obsm["X_umap_mal" if "X_umap_mal" in a.obsm else "X_umap"][:, 0],
        "UMAP2": a.obsm["X_umap_mal" if "X_umap_mal" in a.obsm else "X_umap"][:, 1],
        "pseudotime": pt,
        "pseudotime_rank": pt_rank,
        "branch_id": (a.obs["seg"].astype(str).values
                       if "seg" in a.obs else "0"),
        "milestone": (a.obs["milestones"].astype(str).values
                       if "milestones" in a.obs else ""),
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
    })
    cells_df.to_csv(FIG2/"tree_cells.csv.gz", index=False, compression="gzip")
    log(f"  tree_cells.csv.gz: {cells_df.shape}")

    # Nodes table with cell counts (assign each cell to nearest node)
    from scipy.spatial import cKDTree
    kd = cKDTree(nodes)
    _, node_assign = kd.query(a.obsm["X_emb"], k=1)
    n_cells_per_node = np.bincount(node_assign, minlength=len(nodes))
    nodes_df = pd.DataFrame({
        "node_id": [f"node_{i}" for i in range(len(nodes))],
        "x": nodes_2d[:, 0], "y": nodes_2d[:, 1],
        "n_cells_assigned": n_cells_per_node,
    })
    nodes_df.to_csv(FIG2/"tree_nodes.csv", index=False)
    log(f"  tree_nodes.csv: {nodes_df.shape}")

    log(f"\nbranch breakdown (cells per branch):")
    if "seg" in a.obs:
        log(a.obs["seg"].astype(str).value_counts().to_string())

    log(f"\nMP × branch crosstab (top):")
    if "seg" in a.obs:
        ct = pd.crosstab(a.obs["dominant_MP"], a.obs["seg"])
        log(ct.to_string())
        ct.to_csv(FIG2/"tree_mp_branch_crosstab.csv")

    log(f"\nDONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
