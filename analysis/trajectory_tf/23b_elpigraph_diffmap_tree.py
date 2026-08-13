"""Step 23b: ElPiGraph principal tree on diffusion-map space + manual pseudotime.

Replaces scFates (broken pseudotime). Pipeline:
  1. Load malignant h5ad; rebuild diffmap on X_pca_mal_harmony.
  2. ElPiGraph elastic principal tree on first 4 diffusion components.
  3. Project to 2D (DC1, DC2) for Monocle-style Component_1/Component_2.
  4. Identify root node (closest to MP4-max cell).
  5. Compute pseudotime via shortest-path on tree from root.
  6. Assign each cell to nearest tree node; cell pseudotime = node's distance.
  7. Identify branches (each leaf is the end of a branch from root).

Outputs (overwriting Fig 2/3 trajectory data):
  ${WORK_ROOT}/luad_figures/fig2/tree_cells.csv.gz
  ${WORK_ROOT}/luad_figures/fig2/tree_nodes.csv
  ${WORK_ROOT}/luad_figures/fig2/tree_edges.csv
  ${WORK_ROOT}/luad_figures/fig2/tree_branches.csv
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
N_NODES = 60
N_DM_DIMS = 4   # use first 4 diffusion components for tree fit


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    log("rebuilding neighbors + diffmap")
    sc.pp.neighbors(a, use_rep="X_pca_mal_harmony", n_neighbors=30, random_state=0)
    sc.tl.diffmap(a, n_comps=5, random_state=0)
    DM = a.obsm["X_diffmap"][:, 1:N_DM_DIMS+1]   # skip trivial DC0
    log(f"  diffmap (DC1..DC{N_DM_DIMS}): {DM.shape}")
    log(f"  range DC1=[{DM[:,0].min():.4f},{DM[:,0].max():.4f}] "
        f"DC2=[{DM[:,1].min():.4f},{DM[:,1].max():.4f}]")

    # --- ElPiGraph tree fit on 4D diffmap ---
    log(f"fitting ElPiGraph tree (Nodes={N_NODES}) on {N_DM_DIMS}D diffmap")
    import elpigraph
    t1 = time.time()
    res = elpigraph.computeElasticPrincipalTree(
        DM.astype(np.float64), NumNodes=N_NODES,
        verbose=False, Lambda=0.01, Mu=0.1,
    )
    info = res[0] if isinstance(res, list) else res
    nodes_4d = np.asarray(info["NodePositions"])
    edges = np.asarray(info["Edges"][0]).astype(int)
    log(f"  tree fit done in {time.time()-t1:.1f}s")
    log(f"  nodes={nodes_4d.shape}; edges={edges.shape}")

    # --- Find root: tree node closest to MP4-max cell ---
    mp4 = a.obs["MP4_score"].values
    icell = int(np.argmax(mp4))
    cell_pos = DM[icell]
    from scipy.spatial.distance import cdist
    iroot = int(np.argmin(cdist(cell_pos.reshape(1, -1), nodes_4d).ravel()))
    log(f"  root = node {iroot} (closest to MP4-max cell @ {a.obs.index[icell]})")

    # --- Build NetworkX graph for shortest-path pseudotime ---
    log("computing pseudotime via shortest-path on tree")
    import networkx as nx
    G = nx.Graph()
    for u, v in edges:
        # Edge weight = euclidean distance in 4D diffmap space
        w = float(np.linalg.norm(nodes_4d[u] - nodes_4d[v]))
        G.add_edge(int(u), int(v), weight=w)
    log(f"  graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Pseudotime per node = shortest-path distance from root
    pt_per_node = nx.single_source_dijkstra_path_length(G, iroot)
    node_pt = np.array([pt_per_node.get(i, np.nan) for i in range(len(nodes_4d))])
    log(f"  node pseudotime range: [{np.nanmin(node_pt):.4f},{np.nanmax(node_pt):.4f}]")

    # --- Cell → nearest node assignment + cell pseudotime ---
    from scipy.spatial import cKDTree
    kd = cKDTree(nodes_4d)
    cell_dist, cell_node = kd.query(DM, k=1)
    cell_pt = node_pt[cell_node]
    cell_pt_rank = pd.Series(cell_pt).rank(pct=True).values
    log(f"  cell pseudotime range: [{np.nanmin(cell_pt):.4f},{np.nanmax(cell_pt):.4f}]")

    # --- Branch identification: each leaf = end of one branch ---
    log("identifying branches")
    leaves = [n for n in G.nodes if G.degree(n) == 1 and n != iroot]
    forks = [n for n in G.nodes if G.degree(n) >= 3]
    log(f"  leaves (tips): {leaves}")
    log(f"  forks: {forks}")

    # Each cell's branch = the leaf at the end of its path from root
    # Find which leaf each node leads to (closest leaf in tree)
    branch_of_node = {}
    for leaf in leaves:
        path = nx.shortest_path(G, iroot, leaf)
        for n in path:
            # Multiple branches may share early nodes; later assignment wins
            # Better: each node assigned to the leaf with closest leaf-distance
            cur = branch_of_node.get(n)
            if cur is None:
                branch_of_node[n] = leaf
            else:
                # If this node is in multiple paths, prefer the LATER (closer-to-leaf) one
                d_new = nx.shortest_path_length(G, n, leaf, weight="weight")
                d_old = nx.shortest_path_length(G, n, cur, weight="weight")
                if d_new < d_old:
                    branch_of_node[n] = leaf
    branch_per_node = np.array([branch_of_node.get(i, iroot) for i in range(len(nodes_4d))])
    cell_branch = branch_per_node[cell_node]
    log(f"  branch counts: {pd.Series(cell_branch).value_counts().head(10).to_dict()}")

    # --- Save outputs ---
    log("writing outputs")
    # Tree nodes (project to 2D = DC1, DC2)
    nodes_df = pd.DataFrame({
        "node_id": np.arange(len(nodes_4d)),
        "x": nodes_4d[:, 0],
        "y": nodes_4d[:, 1],
        "DC3": nodes_4d[:, 2] if nodes_4d.shape[1] > 2 else np.nan,
        "DC4": nodes_4d[:, 3] if nodes_4d.shape[1] > 3 else np.nan,
        "node_pseudotime": node_pt,
        "n_cells_assigned": np.bincount(cell_node, minlength=len(nodes_4d)),
        "is_root": np.arange(len(nodes_4d)) == iroot,
        "is_leaf": np.isin(np.arange(len(nodes_4d)), leaves),
        "is_fork": np.isin(np.arange(len(nodes_4d)), forks),
        "branch_leaf": branch_per_node,
    })
    nodes_df.to_csv(FIG2/"tree_nodes.csv", index=False)

    # Tree edges
    edges_df = pd.DataFrame({
        "edge_id": np.arange(len(edges)),
        "src": edges[:, 0], "tgt": edges[:, 1],
        "x_src": nodes_4d[edges[:, 0], 0], "y_src": nodes_4d[edges[:, 0], 1],
        "x_tgt": nodes_4d[edges[:, 1], 0], "y_tgt": nodes_4d[edges[:, 1], 1],
        "length": [G.edges[u, v]["weight"] for u, v in edges],
    })
    edges_df.to_csv(FIG2/"tree_edges.csv", index=False)

    # Branch summary table
    branch_rows = []
    for i, leaf in enumerate(leaves):
        n_cells = int((cell_branch == leaf).sum())
        path = nx.shortest_path(G, iroot, leaf)
        path_len = sum(G.edges[path[k], path[k+1]]["weight"] for k in range(len(path)-1))
        # MP composition along this branch
        mp_in_branch = a.obs["dominant_MP"].astype(str).values[cell_branch == leaf]
        if len(mp_in_branch):
            mp_top = pd.Series(mp_in_branch).value_counts().head(3).to_dict()
        else:
            mp_top = {}
        branch_rows.append({
            "branch_id": f"branch_{i}",
            "leaf_node": leaf, "path_length": path_len,
            "n_path_nodes": len(path),
            "n_cells": n_cells,
            "MP_top3": ";".join(f"{k}:{v}" for k, v in mp_top.items()),
        })
    pd.DataFrame(branch_rows).to_csv(FIG2/"tree_branches.csv", index=False)
    log(f"  tree_branches.csv: {len(branch_rows)} branches")

    # Per-cell output
    umap_key = "X_umap_mal" if "X_umap_mal" in a.obsm else "X_umap"
    cells_df = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "Component_1": DM[:, 0],
        "Component_2": DM[:, 1],
        "UMAP1": a.obsm[umap_key][:, 0],
        "UMAP2": a.obsm[umap_key][:, 1],
        "pseudotime": cell_pt,
        "pseudotime_rank": cell_pt_rank,
        "tree_node": cell_node,
        "branch_leaf": cell_branch,
        "dominant_MP": a.obs["dominant_MP"].astype(str).values,
    })
    cells_df.to_csv(FIG2/"tree_cells.csv.gz", index=False, compression="gzip")
    log(f"  tree_cells.csv.gz: {cells_df.shape}")

    # MP × branch crosstab
    ct = pd.crosstab(cells_df["dominant_MP"], cells_df["branch_leaf"])
    ct.to_csv(FIG2/"tree_mp_branch_crosstab.csv")
    log("\nMP × branch_leaf:")
    log(ct.to_string())

    log(f"\nDONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
