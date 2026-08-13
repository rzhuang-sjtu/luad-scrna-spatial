"""Step 8a: Align LUAD 4 MPs to Gavish 2023 pan-cancer 41 MPs.

Inputs:
  - ~/luad/results/step6_mp_signatures_top100.csv  (LUAD MP top-100)
  - ~/luad/data/reference/gavish2023_MPs.csv       (Gavish top-50)

Outputs:
  - ${WORK_ROOT}/luad_figures/fig2/gavish_overlap.csv         (MP x MP overlap coefficient)
  - ${WORK_ROOT}/luad_figures/fig2/gavish_jaccard.csv         (MP x MP Jaccard)
  - ${WORK_ROOT}/luad_figures/fig2/gavish_top_matches.csv     (each LUAD MP top-5 Gavish hits)
  - ${WORK_ROOT}/luad_figures/fig2/gavish_hypergeom_fdr.csv   (BH-FDR over 4x41 hypergeometric tests)
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

LUAD_SIG = Path.home() / "luad/results/step6_mp_signatures_top100.csv"
GAVISH = Path.home() / "luad/data/reference/gavish2023_MPs.csv"
OUTDIR = Path("${WORK_ROOT}/luad_figures/fig2")
OUTDIR.mkdir(parents=True, exist_ok=True)

TOPN_LUAD = 50  # match Gavish's top-50 so overlap/Jaccard on same footing

luad = pd.read_csv(LUAD_SIG)
luad = luad[luad["rank"] <= TOPN_LUAD].copy()
luad_mps = {mp: set(df.gene) for mp, df in luad.groupby("MP")}

gav = pd.read_csv(GAVISH)
gav_mps = {mp: set(df.gene) for mp, df in gav.groupby("MP", sort=False)}

# gene universe = union of all genes appearing in either side's top lists
universe = set().union(*luad_mps.values(), *gav_mps.values())
N = len(universe)
print(f"LUAD MPs: {list(luad_mps.keys())}")
print(f"Gavish MPs: n={len(gav_mps)}")
print(f"universe N={N}")

rows_overlap, rows_jacc, rows_hyper = [], [], []
for lname, lset in luad_mps.items():
    for gname, gset in gav_mps.items():
        inter = lset & gset
        n_inter = len(inter)
        # overlap coefficient = |A ∩ B| / min(|A|, |B|)
        oc = n_inter / min(len(lset), len(gset)) if min(len(lset), len(gset)) else 0
        # Jaccard = |A ∩ B| / |A ∪ B|
        ja = n_inter / len(lset | gset) if (lset | gset) else 0
        # hypergeometric: pick |gset| from N, how many hit |lset|
        pval = hypergeom.sf(n_inter - 1, N, len(lset), len(gset)) if n_inter > 0 else 1.0
        rows_overlap.append({"LUAD_MP": lname, "Gavish_MP": gname, "overlap_coef": oc, "n_intersect": n_inter})
        rows_jacc.append({"LUAD_MP": lname, "Gavish_MP": gname, "jaccard": ja})
        rows_hyper.append({"LUAD_MP": lname, "Gavish_MP": gname, "pvalue": pval,
                           "n_intersect": n_inter, "shared_genes": ",".join(sorted(inter))})

df_oc = pd.DataFrame(rows_overlap).pivot(index="LUAD_MP", columns="Gavish_MP", values="overlap_coef")
df_oc = df_oc[[c for c in gav_mps.keys()]]  # preserve MP1..MP41 order
df_oc.to_csv(OUTDIR / "gavish_overlap.csv")

df_ja = pd.DataFrame(rows_jacc).pivot(index="LUAD_MP", columns="Gavish_MP", values="jaccard")
df_ja = df_ja[[c for c in gav_mps.keys()]]
df_ja.to_csv(OUTDIR / "gavish_jaccard.csv")

df_h = pd.DataFrame(rows_hyper)
_, fdr, _, _ = multipletests(df_h["pvalue"].values, method="fdr_bh")
df_h["fdr_bh"] = fdr
df_h.to_csv(OUTDIR / "gavish_hypergeom_fdr.csv", index=False)

df_oc_long = pd.DataFrame(rows_overlap)[["LUAD_MP", "Gavish_MP", "overlap_coef"]]
merged = df_h.merge(df_oc_long, on=["LUAD_MP", "Gavish_MP"])
top_hits = (
    merged.sort_values(["LUAD_MP", "overlap_coef", "fdr_bh"], ascending=[True, False, True])
          .groupby("LUAD_MP")
          .head(5)
          .reset_index(drop=True)
)
top_hits.to_csv(OUTDIR / "gavish_top_matches.csv", index=False)

# Console report
print("\n=== Top-3 Gavish matches per LUAD MP ===")
for lm in luad_mps:
    sub = top_hits[top_hits["LUAD_MP"] == lm].head(3)
    print(f"\n{lm}:")
    for _, r in sub.iterrows():
        print(f"  {r['Gavish_MP']:40s} overlap={r['overlap_coef']:.2f} "
              f"n={r['n_intersect']:2d} FDR={r['fdr_bh']:.2e}")

print(f"\n{OUTDIR}/ outputs written.")
