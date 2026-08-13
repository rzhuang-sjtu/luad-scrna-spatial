"""Augment monocle3_graph_segments.csv with per-edge pseudotime and direction."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

BASE = Path("${WORK_ROOT}/luad_figures/fig2")
seg = pd.read_csv(BASE / "monocle3_graph_segments.csv")
cells = pd.read_csv(BASE / "pseudotime_umap.csv.gz")

# unify column names
u1, u2 = "UMAP1", "UMAP2"
if u1 not in cells.columns:
    u1 = next(c for c in cells.columns if c.lower() in ("umap1", "umap_1"))
    u2 = next(c for c in cells.columns if c.lower() in ("umap2", "umap_2"))
pt_col = next(c for c in cells.columns
              if c.lower() in ("pt_winsorized", "pseudotime", "pt"))

cells = cells.dropna(subset=[u1, u2, pt_col])
tree = cKDTree(cells[[u1, u2]].values)

# unique node coordinates: union of (node_start, x_start, y_start) and (node_end, ...)
node_rows = (
    pd.concat([
        seg.rename(columns={"node_start":"node","x_start":"x","y_start":"y"})[["node","x","y"]],
        seg.rename(columns={"node_end":"node","x_end":"x","y_end":"y"})[["node","x","y"]],
    ])
    .drop_duplicates("node")
    .sort_values("node")
    .reset_index(drop=True)
)

K = 25
node_pt = []
for _, r in node_rows.iterrows():
    _, idx = tree.query([r["x"], r["y"]], k=K)
    node_pt.append(float(np.median(cells[pt_col].values[idx])))
node_rows["pseudotime"] = node_pt
node_pt_lookup = dict(zip(node_rows["node"], node_rows["pseudotime"]))

# Per-edge pseudotime + direction
seg["pt_start"] = seg["node_start"].map(node_pt_lookup)
seg["pt_end"]   = seg["node_end"].map(node_pt_lookup)
seg["pt_mean"]  = (seg["pt_start"] + seg["pt_end"]) / 2

# orient each edge so it points from low to high pseudotime
flip = seg["pt_start"] > seg["pt_end"]
for low_col, high_col in [("x_start","x_end"), ("y_start","y_end"),
                           ("pt_start","pt_end"), ("node_start","node_end")]:
    a = seg[low_col].copy(); b = seg[high_col].copy()
    seg.loc[flip, low_col]  = b[flip]
    seg.loc[flip, high_col] = a[flip]
seg.to_csv(BASE / "monocle3_graph_segments_oriented.csv", index=False)

# Annotate root (min pt) + tip (max pt) nodes
root = node_rows.loc[node_rows["pseudotime"].idxmin()]
tips = node_rows[node_rows["pseudotime"] >= node_rows["pseudotime"].quantile(0.97)]
markers = pd.DataFrame({
    "label": ["Root"] + [f"Tip {i+1}" for i in range(len(tips))],
    "x":     [root["x"]] + tips["x"].tolist(),
    "y":     [root["y"]] + tips["y"].tolist(),
    "pt":    [root["pseudotime"]] + tips["pseudotime"].tolist(),
    "kind":  ["root"] + ["tip"] * len(tips),
})
markers.to_csv(BASE / "monocle3_graph_root_tip.csv", index=False)

print(f"node count: {len(node_rows)}, edges: {len(seg)}")
print(f"pt range: {node_rows['pseudotime'].min():.2f} - {node_rows['pseudotime'].max():.2f}")
print(f"root node {int(root['node'])} at ({root['x']:.2f}, {root['y']:.2f}), pt={root['pseudotime']:.2f}")
print(f"{len(tips)} tip(s) at top 3% pseudotime")
print("wrote monocle3_graph_segments_oriented.csv + monocle3_graph_root_tip.csv")
