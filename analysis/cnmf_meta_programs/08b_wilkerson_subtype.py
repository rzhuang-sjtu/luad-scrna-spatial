"""Step 8b: Wilkerson TRU/PP/PI subtype scoring of LUAD malignant cells.

Method: for each of 3 centroids, take genes with max centroid in that subtype
AND margin > 0.5 relative to the second-highest. Score each malignant cell
with sc.tl.score_genes against these 3 up-regulated signatures. Assign
dominant subtype per cell. Cross-tab against dominant MP.

Outputs:
  - ${WORK_ROOT}/luad_figures/fig2/wilkerson_subtype_signatures.csv
  - ${WORK_ROOT}/luad_figures/fig2/wilkerson_cell_scores.csv.gz
  - ${WORK_ROOT}/luad_figures/fig2/wilkerson_MP_crosstab_count.csv
  - ${WORK_ROOT}/luad_figures/fig2/wilkerson_MP_crosstab_pct.csv
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

CENTROIDS = Path.home() / "luad/data/reference/wilkerson_LAD_centroids.csv"
MAL_H5AD = Path.home() / "luad/data/processed/luad_malignant_scored.h5ad"
MP_ASSIGN = Path.home() / "luad/results/step7_mp_cell_scores.csv"
OUTDIR = Path("${WORK_ROOT}/luad_figures/fig2")
OUTDIR.mkdir(parents=True, exist_ok=True)

SUBTYPE_MAP = {"bronchioid": "TRU", "magnoid": "PP", "squamoid": "PI"}
MARGIN_THRESH = 0.5  # centroid margin for signature membership

cent = pd.read_csv(CENTROIDS, index_col=0)
cent.index.name = "gene"
print(f"Wilkerson centroids: {cent.shape}")

# Signature per subtype: gene's max column is this subtype AND margin > thresh
sigs = {}
for sub_raw, sub_std in SUBTYPE_MAP.items():
    others = [c for c in cent.columns if c != sub_raw]
    max_other = cent[others].max(axis=1)
    margin = cent[sub_raw] - max_other
    mask = (cent[sub_raw] == cent.max(axis=1)) & (margin > MARGIN_THRESH)
    sigs[sub_std] = cent.index[mask].tolist()
    print(f"  {sub_std:4s} ({sub_raw:10s}): n_genes={len(sigs[sub_std])}")

sig_long = pd.concat(
    [pd.DataFrame({"subtype": k, "gene": v}) for k, v in sigs.items()], ignore_index=True
)
sig_long.to_csv(OUTDIR / "wilkerson_subtype_signatures.csv", index=False)

print(f"\nLoading malignant h5ad: {MAL_H5AD}")
ad = sc.read_h5ad(MAL_H5AD)
print(f"  shape={ad.shape}  obs cols={list(ad.obs.columns)[:12]}...")

# Ensure log-normalized for score_genes
if "log1p" not in ad.uns:
    if ad.X.max() > 50:
        print("  Applying normalize_total + log1p...")
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
    else:
        print("  X appears pre-log; skipping normalize")

for sub_std, genes in sigs.items():
    present = [g for g in genes if g in ad.var_names]
    print(f"  scoring {sub_std}: {len(present)}/{len(genes)} genes present")
    sc.tl.score_genes(ad, gene_list=present, score_name=f"Wilkerson_{sub_std}",
                      random_state=0, use_raw=False)

score_df = ad.obs[[f"Wilkerson_{s}" for s in SUBTYPE_MAP.values()]].copy()
score_df.columns = list(SUBTYPE_MAP.values())
score_df["dominant_subtype"] = score_df.idxmax(axis=1)

# Merge MP assignment
mp = pd.read_csv(MP_ASSIGN, index_col=0)
mp_col = "dominant_MP" if "dominant_MP" in mp.columns else mp.columns[-1]
print(f"  MP column used: {mp_col}")
score_df = score_df.join(mp[[mp_col]].rename(columns={mp_col: "dominant_MP"}), how="inner")
print(f"  joined rows: {len(score_df)}")

score_df.to_csv(OUTDIR / "wilkerson_cell_scores.csv.gz", compression="gzip")

# Cross-tab (exclude MP5 if present, since it's the discarded cluster)
valid_mps = sorted([m for m in score_df["dominant_MP"].unique() if m != "MP5"])
sub = score_df[score_df["dominant_MP"].isin(valid_mps)]

ct = pd.crosstab(sub["dominant_MP"], sub["dominant_subtype"])
ct = ct.reindex(index=valid_mps, columns=["TRU", "PP", "PI"], fill_value=0)
ct.to_csv(OUTDIR / "wilkerson_MP_crosstab_count.csv")

pct = ct.div(ct.sum(axis=1), axis=0).round(4)
pct.to_csv(OUTDIR / "wilkerson_MP_crosstab_pct.csv")

print("\n=== MP × Wilkerson subtype (counts) ===")
print(ct)
print("\n=== MP × Wilkerson subtype (row %) ===")
print((pct * 100).round(1))
print(f"\nOutputs → {OUTDIR}")
