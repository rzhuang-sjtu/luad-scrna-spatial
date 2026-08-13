#!/usr/bin/env python
"""H7 — Write two neutrophil-related tables (neutrophil subtype composition and ambient-RNA controls).

Archive sections 9.5 and 3.8 report these numbers without CSV; response-letter figures must be file-driven.
  Table 1 detection rates of CXCR1/CXCR2/CSF3R/FCGR3B within each dataset
       — response to 'receptors missing from scRNA-seq': not missing in the data; dropped when intersecting to 9,881 shared genes
  Table 2 per-subtype cell counts, median genes detected, and scANVI uncertainty
"""
import os
import numpy as np, pandas as pd, scanpy as sc

OUT = "${PROJECT_ROOT}/results/neutrophil_tables"
os.makedirs(OUT, exist_ok=True)
a = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad")
ann = sc.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad",
                   backed="r")
print(f"neutrophils {a.shape}", flush=True)

GENES = ["CXCR2", "CXCR1", "CSF3R", "FCGR3B", "S100A8", "OSM"]
ds = a.obs["dataset"].astype(str).values
rows = []
for d in sorted(set(ds)):
    m = ds == d
    r = {"dataset": d, "n": int(m.sum())}
    for g in GENES:
        if g in a.var_names:
            v = a[m, g].X
            v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
            r[g] = float((v > 0).mean() * 100)
        else:
            r[g] = np.nan
    rows.append(r)
D = pd.DataFrame(rows).sort_values("n", ascending=False)
tot = {"dataset": "weighted mean", "n": int(D.n.sum())}
for g in GENES:
    tot[g] = float((D[g] * D.n).sum() / D.n.sum())
D = pd.concat([D, pd.DataFrame([tot])], ignore_index=True)
D.to_csv(f"{OUT}/cxcr_detection_by_dataset.csv", index=False)
print(D.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

o = ann.obs
S = o.groupby("neu_subtype", observed=True).agg(
    n=("neu_subtype", "size"),
    median_genes=("scanvi_uncertainty", "size"),
    scanvi_uncertainty=("scanvi_uncertainty", "median")).reset_index()
if "n_genes" in o.columns:
    S["median_genes"] = o.groupby("neu_subtype", observed=True)["n_genes"].median().values
else:
    ng = np.asarray((a.X > 0).sum(1)).ravel()
    tmp = pd.DataFrame({"s": o["neu_subtype"].astype(str).values[:len(ng)], "g": ng})
    S["median_genes"] = tmp.groupby("s").g.median().reindex(S.neu_subtype).values
S = S.sort_values("n", ascending=False)
S.to_csv(f"{OUT}/subtype_summary.csv", index=False)
print()
print(S.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print(f"\nWrote {OUT}/")
