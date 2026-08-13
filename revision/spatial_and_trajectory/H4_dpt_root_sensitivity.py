#!/usr/bin/env python
"""
H4 — Sensitivity of diffusion pseudotime to root-cell choice.

Motivation: the root cell for diffusion pseudotime is the cell with the highest MP4
score. That is biologically motivated but still arbitrary, and the sensitivity of the
trajectory topology to the choice of root had not been assessed.

Also: MP3 score does not increase monotonically, yet the text risks implying a
directional MP4-to-MP3 trajectory.

Design:
  Main-text root = single cell with highest MP4 score.
  Alternative roots:
    A. 10 cells drawn at random from the top 1% MP4 (same biological definition; only the specific cell changes)
    B. "Centre" of high-MP4 cells (nearest to the top-1% MP4 centroid) — more stable than a single extreme
    C. Highest-scoring cell for MP1 / MP2 / MP3 each — different biological starting points
    D. 5 cells drawn at random from all cells — no prior

Report two things:
  1. **Topological stability**: Spearman correlation among pseudotimes from different roots.
     If correlations under the same biological definition (A) are high, conclusions do not depend on which specific cell is chosen.
  2. **Whether the directional claim holds**: correlation of MP3 with pseudotime under each root, plus binned mean curves.
     If MP3 does not rise monotonically with pseudotime, the text must not imply an MP4→MP3 direction.

Usage: python H4_dpt_root_sensitivity.py [--n-comps 15]
"""
import os, argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats as sps

H5 = "${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad"
OUT = "${PROJECT_ROOT}/results/dpt_root_sensitivity"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-comps", type=int, default=15)
    ap.add_argument("--seed", type=int, default=20260805)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    print("Loading malignant-cell object ...", flush=True)
    ad = sc.read_h5ad(H5)
    print(f"  {ad.shape}", flush=True)

    # Neighbourhood graph present; compute diffusion map directly
    if "X_diffmap" not in ad.obsm:
        print("Computing diffusion map ...", flush=True)
        sc.tl.diffmap(ad, n_comps=a.n_comps)
    mp = {m: ad.obs[f"{m}_score"].values.astype(float)
          for m in ["MP1", "MP2", "MP3", "MP4"]}

    roots = {}
    roots["main text: highest MP4"] = int(np.argmax(mp["MP4"]))

    top = np.where(mp["MP4"] >= np.quantile(mp["MP4"], 0.99))[0]
    for i, c in enumerate(rng.choice(top, 10, replace=False), 1):
        roots[f"A{i}: random from top 1% MP4"] = int(c)

    # B: nearest to top-1% MP4 centroid in diffusion space; more stable than a single extreme
    D = ad.obsm["X_diffmap"][:, 1:a.n_comps]
    cen = D[top].mean(0)
    roots["B: centre of high-MP4 cells"] = int(top[np.argmin(((D[top] - cen) ** 2).sum(1))])

    for m in ["MP1", "MP2", "MP3"]:
        roots[f"C: highest {m}"] = int(np.argmax(mp[m]))
    for i, c in enumerate(rng.choice(len(D), 5, replace=False), 1):
        roots[f"D{i}: random from all cells"] = int(c)

    print(f"\n{len(roots)} roots; computing DPT for each ...", flush=True)
    pt = {}
    for name, r in roots.items():
        ad.uns["iroot"] = r
        sc.tl.dpt(ad, n_dcs=a.n_comps)
        v = ad.obs["dpt_pseudotime"].values.astype(float)
        v[~np.isfinite(v)] = np.nan
        pt[name] = v
    P = pd.DataFrame(pt)
    P.to_csv(f"{OUT}/pseudotime_by_root.csv.gz", index=False, compression="gzip")

    ref = P.columns[0]
    rows = []
    for c in P.columns:
        ok = P[ref].notna() & P[c].notna()
        rho = sps.spearmanr(P.loc[ok, ref], P.loc[ok, c])[0]
        rows.append(dict(根=c, 与正文根的Spearman=rho,
                         有效细胞=int(ok.sum())))
    S = pd.DataFrame(rows)
    print("\n【1. Pseudotime from each root vs main-text root】")
    print(S.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    grpA = [c for c in P.columns if c.startswith("A")]
    if grpA:
        sub = P[[ref] + grpA].dropna()
        cc = sub.corr(method="spearman").values
        iu = np.triu_indices_from(cc, 1)
        print(f"\n  Within the same biological definition (main-text root + {len(grpA)} random top-1% MP4)"
              f"pairwise Spearman: median {np.median(cc[iu]):.3f},"
              f"range {cc[iu].min():.3f}–{cc[iu].max():.3f}")

    print("\n【2. MP3 vs pseudotime (testing the MP4→MP3 directional claim)】")
    rows2 = []
    for c in P.columns:
        ok = P[c].notna()
        r3 = sps.spearmanr(P.loc[ok, c], mp["MP3"][ok.values])[0]
        r4 = sps.spearmanr(P.loc[ok, c], mp["MP4"][ok.values])[0]
        # Whether mean MP3 is monotone after decile binning
        q = pd.qcut(P.loc[ok, c], 10, labels=False, duplicates="drop")
        m3 = pd.Series(mp["MP3"][ok.values]).groupby(q.values).mean()
        mono = bool((m3.diff().dropna() > 0).all())
        rows2.append(dict(根=c, MP3_rho=r3, MP4_rho=r4,
                          MP3十分位单调上升=mono,
                          MP3首末分位差=float(m3.iloc[-1] - m3.iloc[0])))
    T = pd.DataFrame(rows2)
    print(T.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    S.merge(T, on="根").to_csv(f"{OUT}/root_sensitivity_summary.csv", index=False)

    print(f"\n  Roots with positive MP3–pseudotime correlation: {int((T.MP3_rho>0).sum())}/{len(T)}")
    print(f"Roots with strictly monotone increasing MP3 deciles: {int(T.MP3十分位单调上升.sum())}/{len(T)}")
    print(f"\nWrote {OUT}/", flush=True)


if __name__ == "__main__":
    main()
