#!/usr/bin/env python
"""
C1c — Multi-metric robustness of ranking reproducibility.

Motivation: top-200 Jaccard is **one metric with an arbitrary cutoff** (why 200 not 50?).
It alone is not a reliable basis for conclusions. Re-test the same split-half results with six statistics,
three of which need no cutoff. If a claim holds under only one metric, do not put it in the revision.

Metrics:
  1. Jaccard          |A∩B| / |A∪B|, scanned over K
  2. Overlap coefficient          |A∩B| / min(|A|,|B|); more robust when gene-set sizes differ
  3. Hypergeometric test        p-value for top-K overlap (population = union of significant genes from both halves)
  4. RBO              rank-biased overlap, **no cutoff**; rank-weighted, higher ranks weigh more
  5. Spearman         rank correlation of ΔS on shared genes, **no cutoff**
  6. Kendall tau      same, more robust to ties, **no cutoff**

Null model: keep each half's significant gene set fixed; shuffle ranks within each set only.
This asks whether ranking itself is informative, not whether the halves share significant genes—
the latter would treat filter-driven set overlap as ranking reproducibility.

Usage: python C1c_stability_robustness.py [cache_dir] [output_prefix]
"""
import os, sys, pickle
import numpy as np
import pandas as pd
from scipy import stats as sps

CACHE = sys.argv[1] if len(sys.argv) > 1 else (
    "${PROJECT_ROOT}/results/fig8_geneformer/rank_stability/"
    "_cache_rent_macro1500_macro_spp1_to_c1qc")
DEST = sys.argv[2] if len(sys.argv) > 2 else (
    "${PROJECT_ROOT}/results/fig8_geneformer/rank_stability/rent_macro1500")
KS = [20, 50, 100, 200, 300, 500]
RBO_PS = [0.90, 0.95, 0.98]
NPERM, SEED = 200, 0


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a or b) else np.nan


def overlap_coef(a, b):
    a, b = set(a), set(b)
    m = min(len(a), len(b))
    return len(a & b) / m if m else np.nan


def rbo(l1, l2, p, depth=None):
    """Rank-biased overlap (Webber 2010), non-extrapolated form.
    Rank-weighted geometric decay; no cutoff: smaller p focuses more on the head."""
    d = depth or min(len(l1), len(l2))
    if d == 0:
        return np.nan
    s1, s2, agree, total = set(), set(), 0.0, 0.0
    for i in range(d):
        s1.add(l1[i]); s2.add(l2[i])
        w = p ** i
        agree += w * len(s1 & s2) / (i + 1)
        total += w
    return agree / total if total else np.nan


def main():
    rng = np.random.default_rng(SEED)
    pairs = []
    for r in range(100):
        fa, fb = f"{CACHE}/sh_{r}_A.pkl", f"{CACHE}/sh_{r}_B.pkl"
        if os.path.exists(fa) and os.path.exists(fb):
            with open(fa, "rb") as fh:
                A = pickle.load(fh)
            with open(fb, "rb") as fh:
                B = pickle.load(fh)
            pairs.append((r, A, B))
    if not pairs:
        sys.exit(f"No sh_*_A.pkl / sh_*_B.pkl in cache dir: {CACHE}")
    print(f"Split-half pairs: {len(pairs)}", flush=True)

    rows = []
    for r, A, B in pairs:
        ga, gb = A.Gene_name.values, B.Gene_name.values
        univ = len(set(ga) | set(gb))
        for K in KS:
            a, b = ga[:K], gb[:K]
            ka, kb = len(a), len(b)
            inter = len(set(a) & set(b))
            # Hypergeometric: population univ, success state kb, draw ka, observed inter
            p_hyp = sps.hypergeom.sf(inter - 1, univ, kb, ka) if ka and kb else np.nan
            nj, no = [], []
            for _ in range(NPERM):
                pa = rng.permutation(ga)[:K]; pb = rng.permutation(gb)[:K]
                nj.append(jaccard(pa, pb)); no.append(overlap_coef(pa, pb))
            rows.append(dict(rep=r, K=K, n_A=ka, n_B=kb, universe=univ, inter=inter,
                             jaccard=jaccard(a, b), jaccard_null=np.mean(nj),
                             overlap=overlap_coef(a, b), overlap_null=np.mean(no),
                             p_hypergeom=p_hyp))
    S = pd.DataFrame(rows)
    S["jaccard_ratio"] = S.jaccard / S.jaccard_null
    S["overlap_ratio"] = S.overlap / S.overlap_null
    S.to_csv(f"{DEST}_robustness_bykt.csv", index=False)

    g = S.groupby("K").agg(
        Jaccard=("jaccard", "mean"), J_null=("jaccard_null", "mean"),
        J倍数=("jaccard_ratio", "mean"),
        重叠系数=("overlap", "mean"), overlap_null=("overlap_null", "mean"),
        重叠倍数=("overlap_ratio", "mean"),
        超几何p中位=("p_hypergeom", "median"))
    print("\n[1-3] Scan by top-K (mean over 10 pairs)")
    print(g.to_string(float_format=lambda x: f"{x:.4g}"))

    rows2 = []
    for r, A, B in pairs:
        ga, gb = list(A.Gene_name.values), list(B.Gene_name.values)
        m = A[["Gene_name", "Shift_to_goal_end"]].merge(
            B[["Gene_name", "Shift_to_goal_end"]], on="Gene_name", suffixes=("_a", "_b"))
        sp = sps.spearmanr(m.Shift_to_goal_end_a, m.Shift_to_goal_end_b)[0] if len(m) > 20 else np.nan
        kt = sps.kendalltau(m.Shift_to_goal_end_a, m.Shift_to_goal_end_b)[0] if len(m) > 20 else np.nan
        # Null for rank correlations: shuffle one column
        nsp, nkt = [], []
        v = m.Shift_to_goal_end_b.values
        for _ in range(NPERM):
            q = rng.permutation(v)
            nsp.append(sps.spearmanr(m.Shift_to_goal_end_a, q)[0])
            nkt.append(sps.kendalltau(m.Shift_to_goal_end_a, q)[0])
        row = dict(rep=r, n_common=len(m), spearman=sp, spearman_null=np.mean(nsp),
                   kendall=kt, kendall_null=np.mean(nkt))
        for p in RBO_PS:
            row[f"rbo_p{p}"] = rbo(ga, gb, p)
            row[f"rbo_p{p}_null"] = np.mean(
                [rbo(list(rng.permutation(ga)), list(rng.permutation(gb)), p)
                 for _ in range(20)])          # RBO is slower; fewer null draws
        rows2.append(row)
    T = pd.DataFrame(rows2)
    T.to_csv(f"{DEST}_robustness_nocut.csv", index=False)

    print("\n[4-6] Metrics that need no cutoff (mean over 10 pairs)")
    print(f"Shared genes            {T.n_common.mean():.0f}")
    print(f"Spearman  observed {T.spearman.mean():.3f}   null {T.spearman_null.mean():+.3f}")
    print(f"Kendall   observed {T.kendall.mean():.3f}   null {T.kendall_null.mean():+.3f}")
    for p in RBO_PS:
        print(f"RBO p={p}  observed {T[f'rbo_p{p}'].mean():.3f}"
              f"null {T[f'rbo_p{p}_null'].mean():.3f}"
              f"fold {T[f'rbo_p{p}'].mean()/T[f'rbo_p{p}_null'].mean():.2f}×")

    # Abstract set overlap is not the same as whether selected genes are trustworthy.
    # Ask directly: in how many of 10 split-half pairs does a gene enter top-K on both halves.
    rec = {}
    for K in [50, 200]:
        cnt = {}
        for r, A, B in pairs:
            both = set(A.Gene_name.values[:K]) & set(B.Gene_name.values[:K])
            for gname in both:
                cnt[gname] = cnt.get(gname, 0) + 1
        rec[K] = pd.Series(cnt).sort_values(ascending=False)
    R = pd.DataFrame({f"两半同时进 top{K}的组数": rec[K] for K in rec}).fillna(0).astype(int)
    R.index.name = "gene"
    R.to_csv(f"{DEST}_gene_recurrence.csv")
    print(f"\n[7] Gene recurrence ({len(pairs)} split-half pairs)")
    print(f"Genes in both-halves top200 in all 10/10 pairs:"
          f"{int((R['两半同时进 top200的组数'] == len(pairs)).sum())}")
    print(f"  ≥8/10 组的: {int((R['两半同时进 top200的组数'] >= 0.8*len(pairs)).sum())}")
    print(f"≥8/10 pairs and also both-halves top50:"
          f"{int(((R.get('两半同时进 top50的组数', 0) >= 0.8*len(pairs))).sum())}")
    print("\n  Top 15 most stable genes:")
    print(R.sort_values(list(R.columns), ascending=False).head(15).to_string())
    print(f"\nWriting {DEST}_robustness_bykt.csv / _robustness_nocut.csv / _gene_recurrence.csv")


if __name__ == "__main__":
    main()
