#!/usr/bin/env python
"""
C1d — Confidence intervals for ranks of specified genes (repeat count configurable).

Motivation: the rank «SEC61G 143/688» from `ref.pkl` is a **single point estimate** with no error range.
Using it alone in the revision treats a point estimate as a conclusion.

Two resampling schemes, both run (each has bias):
  bootstrap (with replacement, n=all)
      duplicates cells (pseudo-replication) and inflates significance — genes passing the filter rose from 688 to 750–1236 in practice,
      so absolute ranks are not comparable; report **percentiles** instead.
  subsample (without replacement, n=frac×all)
      no cell duplication, no that inflation; cost is smaller effective n, which slightly lowers detection.

Report both; if a claim holds under only one scheme, do not put it in the revision.

Usage:
  ~/miniforge3/envs/geneformer/bin/python C1d_gene_rank_ci.py \
    --perturb-dir ${PROJECT_ROOT}/results/fig8_geneformer/perturb_rent_macro \
    --transition macro_spp1_to_c1qc --tag macro1500 \
    --genes SEC61G SEC61B ENO2 ANGPTL4 CTSK ALDOA \
    --n-boot 50 --n-sub 50 --sub-frac 0.8 --jobs 8
"""
import os, sys, time, argparse

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import multiprocessing as mp                                      # noqa: E402
import pickle                                                     # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from C1_rank_stability import cell_files, stats_for, OUT          # noqa: E402

_CF = _CACHE = None


def _init(cf, cache):
    global _CF, _CACHE
    _CF, _CACHE = cf, cache


def _run(job):
    key, cells, dup = job
    p = os.path.join(_CACHE, "_".join(map(str, key)) + ".pkl")
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return key, pickle.load(fh)
    d = stats_for(cells, _CF, dup_tag=dup)
    tmp = p + ".part"
    with open(tmp, "wb") as fh:
        pickle.dump(d, fh)
    os.replace(tmp, p)
    return key, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturb-dir", required=True)
    ap.add_argument("--transition", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--n-boot", type=int, default=50)
    ap.add_argument("--n-sub", type=int, default=50)
    ap.add_argument("--sub-frac", type=float, default=0.8)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260805)
    a = ap.parse_args()
    if not os.path.isabs(a.perturb_dir):
        sys.exit("--perturb-dir must be an absolute path")

    cf = cell_files(os.path.join(a.perturb_dir, a.transition))
    cells = np.array(sorted(cf))
    rng = np.random.default_rng(a.seed)
    nsub = int(round(a.sub_frac * len(cells)))

    jobs = [(("full",), cells, None)]
    for r in range(a.n_boot):
        jobs.append((("boot", r), rng.choice(cells, len(cells), replace=True), r))
    for r in range(a.n_sub):
        jobs.append((("sub", r), rng.choice(cells, nsub, replace=False), None))

    cache = os.path.join(OUT, f"_cache_ci_{a.tag}_{a.transition}")
    os.makedirs(cache, exist_ok=True)
    hit = sum(os.path.exists(os.path.join(cache, "_".join(map(str, j[0])) + ".pkl"))
              for j in jobs)
    print(f"{a.transition}: {len(cells)} cells, {len(jobs)} rounds (cache hits {hit}),"
          f"{a.jobs}-way parallel", flush=True)

    res, t0, done = {}, time.time(), 0
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init,
                             initargs=(cf, cache),
                             mp_context=mp.get_context("spawn")) as ex:
        for f in as_completed([ex.submit(_run, j) for j in jobs]):
            k, d = f.result()
            res[k] = d
            done += 1
            if done % 10 == 0 or done == len(jobs):
                print(f"{done}/{len(jobs)}  elapsed {(time.time()-t0)/60:.1f} min", flush=True)

    rows = []
    for k, d in res.items():
        g = list(d.Gene_name.values)
        kind = k[0]
        for gene in a.genes:
            r = g.index(gene) + 1 if gene in g else np.nan
            rows.append(dict(kind=kind, rep=k[1] if len(k) > 1 else -1, gene=gene,
                             n_pass=len(g), rank=r,
                             pct=r / len(g) * 100 if r == r else np.nan))
    R = pd.DataFrame(rows)
    dest = f"{OUT}/{a.tag}_{a.transition}_gene_rank_ci.csv"
    R.to_csv(dest, index=False)

    print("\nRank percentile per gene (smaller = higher rank) — absolute ranks not comparable because n genes passing the filter differs across rounds")
    for kind in ["full", "boot", "sub"]:
        s = R[R.kind == kind]
        if s.empty:
            continue
        lab = {"full": "Full single run", "boot": f"bootstrap 有放回 ×{a.n_boot}",
               "sub": f"Subsample without replacement {a.sub_frac:.0%} ×{a.n_sub}"}[kind]
        print(f"\n[{lab}] genes passing filter: median {s.n_pass.median():.0f}")
        t = s.groupby("gene").agg(
            通过次数=("rank", lambda v: int(v.notna().sum())),
            总次数=("rank", "size"),
            百分位中位=("pct", "median"),
            百分位2_5=("pct", lambda v: np.nanpercentile(v, 2.5) if v.notna().any() else np.nan),
            百分位97_5=("pct", lambda v: np.nanpercentile(v, 97.5) if v.notna().any() else np.nan),
            排名中位=("rank", "median"))
        t = t.reindex(a.genes)
        print(t.to_string(float_format=lambda x: f"{x:.1f}"))
    print(f"\nWriting {dest}")


if __name__ == "__main__":
    main()
