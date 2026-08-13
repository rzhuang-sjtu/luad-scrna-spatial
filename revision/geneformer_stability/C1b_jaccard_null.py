#!/usr/bin/env python
"""
C1b — Null model for split-half top-200 Jaccard.

Why needed: J=0.337 alone cannot be judged high or low.
Two common misreadings both fail:
  «0.337 is only one third, so ranking is unstable» — incorrect.
  «0.337 far exceeds random expectation 200²/9400 ≈ 0.01» — also incorrect; that denominator uses all tested genes,
    but only the few hundred genes that pass significance filtering compete for top-200, and the two halves' significant sets
    already overlap substantially (~250 of ~450 shared in practice), so random ranks still yield non-trivial overlap.

Correct null model: **keep each half's significant gene set fixed and shuffle ranks within each set only**,
then take top-200 Jaccard. This asks whether **ranking itself is informative**,
not whether the two halves share significant genes.

Reads per-replicate caches from C1p (sh_<rep>_A.pkl / sh_<rep>_B.pkl); no need to re-run stats.

Usage: python C1b_jaccard_null.py <cache_dir> [output csv]
"""
import os, sys, pickle
import numpy as np
import pandas as pd

CACHE = sys.argv[1] if len(sys.argv) > 1 else (
    "${PROJECT_ROOT}/results/fig8_geneformer/rank_stability/"
    "_cache_rent_macro1500_macro_spp1_to_c1qc")
DEST = sys.argv[2] if len(sys.argv) > 2 else (
    "${PROJECT_ROOT}/results/fig8_geneformer/rank_stability/"
    "rent_macro1500_splithalf_null.csv")
TOP, NPERM, SEED = 200, 300, 0


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a or b) else np.nan


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for r in range(100):
        fa, fb = f"{CACHE}/sh_{r}_A.pkl", f"{CACHE}/sh_{r}_B.pkl"
        if not (os.path.exists(fa) and os.path.exists(fb)):
            continue
        with open(fa, "rb") as fh:
            ga = pickle.load(fh).Gene_name.values
        with open(fb, "rb") as fh:
            gb = pickle.load(fh).Gene_name.values
        obs = jac(ga[:TOP], gb[:TOP])
        null = np.array([jac(rng.permutation(ga)[:TOP], rng.permutation(gb)[:TOP])
                         for _ in range(NPERM)])
        rows.append(dict(rep=r, n_sig_A=len(ga), n_sig_B=len(gb),
                         n_sig_overlap=len(set(ga) & set(gb)),
                         J_obs=obs, J_null_mean=null.mean(),
                         J_null_p95=float(np.percentile(null, 95)),
                         ratio=obs / null.mean() if null.mean() else np.nan,
                         above_null_p95=bool(obs > np.percentile(null, 95))))
    if not rows:
        sys.exit(f"No sh_*_A.pkl / sh_*_B.pkl in cache dir: {CACHE}")
    R = pd.DataFrame(rows)
    R.to_csv(DEST, index=False)
    print(R.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print()
    print(f"Mean observed J {R.J_obs.mean():.3f}   mean null J {R.J_null_mean.mean():.3f}"
          f"fold {R.ratio.mean():.2f}×")
    print(f"Above null 95th percentile: {int(R.above_null_p95.sum())} / {len(R)}")
    print(f"\nWriting {DEST}")


if __name__ == "__main__":
    main()
