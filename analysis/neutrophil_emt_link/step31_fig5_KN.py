"""step31: Build Fig 5 K-N data — pseudotime trajectory + GO/Hallmark enrichment.

Outputs (added to ~/luad/results/fig5_plot_data/):
  fig5k_pseudotime_umap.csv.gz       (barcode, UMAP1, UMAP2, dpt_pseudotime, neu_subtype)
  fig5k_paga_connectivity.csv        (source_subtype, target_subtype, connectivity)
  fig5k_paga_positions.csv           (subtype, x, y, n_cells)
  fig5l_enrichment_neu_inflammatory.csv
  fig5m_enrichment_neu_metastatic.csv
  fig5n_enrichment_neu_ecm_remodeling.csv

Trajectory: PAGA on scANVI neighbors → diffusion pseudotime rooted at the
            highest-confidence Neu_unclassified cell (proxy for "naive" neutrophil).
Enrichment: gseapy.enrichr on top-200 wilcoxon DE genes per subtype.
            Libraries: GO_Biological_Process_2023, MSigDB_Hallmark_2020.
"""
import os, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

t0 = time.time()
DATA = Path("${PROJECT_ROOT}/data/processed")
RES  = Path("${PROJECT_ROOT}/results")
OUT  = RES / "fig5_plot_data"
OUT.mkdir(exist_ok=True, parents=True)

# 5K — Trajectory (PAGA + DPT)
print("[5K] PAGA + DPT")
a = sc.read_h5ad(DATA / "luad_neutrophil_own_annotated.h5ad")
print(f"  loaded: {a.shape}; obsm: {list(a.obsm.keys())}")

a.obs["neu_subtype"] = a.obs["neu_subtype"].astype(str)
# rebuild kNN graph on scANVI latent (more stable than reusing existing)
sc.pp.neighbors(a, use_rep="X_scANVI", n_neighbors=15, random_state=0)

# PAGA on neu_subtype categories
sc.tl.paga(a, groups="neu_subtype")
# extract connectivity matrix
conn = a.uns["paga"]["connectivities"].toarray()
groups = a.obs["neu_subtype"].astype("category").cat.categories.tolist()
paga_long = []
for i in range(len(groups)):
    for j in range(i + 1, len(groups)):
        if conn[i, j] > 0.01:
            paga_long.append({
                "source_subtype": groups[i],
                "target_subtype": groups[j],
                "connectivity":   float(conn[i, j]),
            })
paga_df = pd.DataFrame(paga_long).sort_values("connectivity", ascending=False)
paga_df.to_csv(OUT / "fig5k_paga_connectivity.csv", index=False)
print(f"  PAGA edges (conn>0.01): {len(paga_df)}")

# Position each subtype as median UMAP coordinate (R will use these for nodes)
um = a.obsm["X_umap"]
pos = []
for g in groups:
    mask = a.obs["neu_subtype"].values == g
    pos.append({
        "subtype": g,
        "x": float(np.median(um[mask, 0])),
        "y": float(np.median(um[mask, 1])),
        "n_cells": int(mask.sum()),
    })
pos_df = pd.DataFrame(pos)
pos_df.to_csv(OUT / "fig5k_paga_positions.csv", index=False)

# Diffusion pseudotime: root at most-typical Neu_unclassified cell (median UMAP dist within group)
sc.tl.diffmap(a, n_comps=15, random_state=0)
unclass_mask = a.obs["neu_subtype"].values == "Neu_unclassified"
if unclass_mask.sum() > 0:
    # pick cell closest to median X_scANVI within Neu_unclassified
    sub_lat = np.asarray(a.obsm["X_scANVI"])[unclass_mask]
    centroid = sub_lat.mean(axis=0)
    dists = np.linalg.norm(sub_lat - centroid, axis=1)
    root_local = int(np.argmin(dists))
    root_idx = int(np.where(unclass_mask)[0][root_local])
else:
    root_idx = 0
a.uns["iroot"] = root_idx
print(f"  iroot = {root_idx} (Neu_unclassified centroid cell)")
sc.tl.dpt(a)

dpt_df = pd.DataFrame({
    "barcode":         a.obs.index.astype(str),
    "UMAP1":           um[:, 0],
    "UMAP2":           um[:, 1],
    "neu_subtype":     a.obs["neu_subtype"].astype(str).values,
    "dpt_pseudotime":  a.obs["dpt_pseudotime"].astype(float).values,
})
dpt_df.to_csv(OUT / "fig5k_pseudotime_umap.csv.gz", index=False, compression="gzip")
print(f"  pseudotime range: {dpt_df['dpt_pseudotime'].min():.3f} - {dpt_df['dpt_pseudotime'].max():.3f}")

# 5L/M/N — GO BP + Hallmark enrichment per anchor subtype
print("\n[5L/M/N] GO + Hallmark enrichment")
de = pd.read_csv(RES / "step25e_markers_scanvi_fullgene.csv")
print(f"  DE rows: {len(de)}; groups: {sorted(de['group'].unique())}")

# Map original Salcher labels → functional names for the 3 anchors
ANCHORS = [
    ("TAN-1", "Neu_Inflammatory", "fig5l_enrichment_neu_inflammatory.csv"),
    ("TAN-4", "Neu_Metastatic",   "fig5m_enrichment_neu_metastatic.csv"),
    ("NAN-1", "Neu_ECM_remodeling", "fig5n_enrichment_neu_ecm_remodeling.csv"),
]

import gseapy as gp
LIB = ["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
N_TOP = 200

for old, new, fname in ANCHORS:
    sub = de[de["group"] == old].copy()
    if len(sub) == 0:
        print(f"  WARN: no DE rows for {old} ({new})")
        continue
    # take top-N by score, padj<0.05 + logfc>0
    sub = sub[(sub["padj"] < 0.05) & (sub["logfc"] > 0)].sort_values("score", ascending=False)
    genes = sub["gene"].head(N_TOP).tolist()
    if len(genes) < 10:
        # relax: top-N by score regardless
        genes = de[de["group"] == old].sort_values("score", ascending=False)["gene"].head(N_TOP).tolist()
    print(f"  {new}: {len(genes)} genes → enrichr")

    enr_rows = []
    for lib in LIB:
        try:
            enr = gp.enrichr(gene_list=genes, gene_sets=lib,
                              organism="human", outdir=None, no_plot=True)
            r = enr.results.copy()
            r["library"] = lib
            r["subtype"] = new
            enr_rows.append(r)
        except Exception as e:
            print(f"    {lib} failed: {e}")
    if not enr_rows:
        continue
    out = pd.concat(enr_rows, ignore_index=True)
    # keep informative columns
    keep_cols = ["subtype", "library", "Term", "Overlap", "P-value",
                  "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"]
    keep_cols = [c for c in keep_cols if c in out.columns]
    out = out[keep_cols].copy()
    out.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in out.columns]
    # filter to padj<0.05 then top 30 per library by combined_score
    out = out[out.get("adjusted_p_value", 1) < 0.05].copy() if "adjusted_p_value" in out.columns else out
    out["minus_log10_padj"] = -np.log10(out["adjusted_p_value"].clip(lower=1e-300)) \
        if "adjusted_p_value" in out.columns else np.nan
    if "combined_score" in out.columns:
        out = out.sort_values(["library", "combined_score"], ascending=[True, False])
    # keep top 30 per library
    out = out.groupby("library").head(30).reset_index(drop=True)
    out.to_csv(OUT / fname, index=False)
    print(f"  {new}: {len(out)} significant terms → {fname}")

print(f"\nelapsed: {(time.time()-t0)/60:.1f} min")
print("DONE.")
