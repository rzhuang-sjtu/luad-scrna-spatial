#!/usr/bin/env python
"""
C1p — Parallel version of C1_rank_stability. Computation identical to C1; only scheduling changes.

Why: each stats round in C1 (each split-half half, each saturation subset, each bootstrap)
is fully independent and runs InSilicoPerturberStats in its own temp dir. C1 is serial;
~2.4 min per round on 1 core — 31 of 32 threads idle; 1500-cell run took ~4 h.

**Numerics match the serial version bit-for-bit**: random subsets are still drawn in the main process with the same seed and the same
call order as C1 (including details such as not consuming RNG when `2*n > len(cells)`), then dispatched to the pool.
Parallelism only changes who computes, not what is computed.

Usage (two constraints):
  Must use the geneformer env — transformers in scst is incompatible with geneformer
  --perturb-dir must be absolute — stats_for symlinks pickles into /tmp; relative paths do not resolve

  ~/miniforge3/envs/geneformer/bin/python C1p_rank_stability_parallel.py \
    --perturb-dir ${PROJECT_ROOT}/results/fig8_geneformer/perturb_rent_macro \
    --transitions macro_spp1_to_c1qc --tag rent_macro1500 --jobs 12 \
    --sat-ns 100 200 300 500 750 --agree-ns 100 200 300 500 750 1000 1250
"""
import os, sys, time, argparse

# Must pin BLAS thread count before import numpy.
# This workload is mostly pandas and loops; each worker uses ~1 core in practice;
# but if any path hits large matrix ops, OpenBLAS opens one thread per core by default —
# 28 processes × 32 threads ≈ 900 threads contending; slower than serial.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import multiprocessing as mp                                      # noqa: E402
import pickle                                                     # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from C1_rank_stability import (            # noqa: E402  reuse all compute logic; do not rewrite
    cell_files, stats_for, jac, rank_corr, OUT, TOPN,
)

_CF = None          # per-worker {cell: [pickle paths]}
_CACHE = None       # on-disk cache directory for per-round results


def _init(cf, cache):
    global _CF, _CACHE
    _CF, _CACHE = cf, cache


def _key_path(key):
    return os.path.join(_CACHE, "_".join(map(str, key)) + ".pkl")


def _run(job):
    """job = (key, cell-index array, dup_tag); returns (key, filtered rank table)

    Disk cache for results: if 82 rounds fail at round 70, all progress is lost
    (SIGBUS / stack segment observed; suspected fork of a parent that already loaded torch).
    With a cache the run can resume; re-running skips completed rounds.
    """
    key, cells, dup = job
    p = _key_path(key)
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return key, pickle.load(fh)
    d = stats_for(cells, _CF, dup_tag=dup)
    tmp = p + ".part"
    with open(tmp, "wb") as fh:
        pickle.dump(d, fh)
    os.replace(tmp, p)          # atomic replace so partial files are not treated as valid cache
    return key, d


def build_jobs(cells, a):
    """Consume RNG in the same order as C1 so subsets match exactly."""
    rng = np.random.default_rng(20260804)
    jobs = [(("ref",), cells, None)]

    for r in range(a.n_splithalf):
        p = rng.permutation(cells); h = len(cells) // 2
        jobs.append((("sh", r, "A"), p[:h], None))
        jobs.append((("sh", r, "B"), p[h:2 * h], None))

    for n in a.sat_ns:
        if 2 * n > len(cells):        # C1 continues before drawing RNG; keep the same behaviour here
            continue
        for r in range(a.sat_reps):
            p = rng.permutation(cells)
            jobs.append((("sat", n, r, "A"), p[:n], None))
            jobs.append((("sat", n, r, "B"), p[n:2 * n], None))

    for n in a.agree_ns:
        if n >= len(cells):
            continue
        for r in range(a.sat_reps):
            s = rng.choice(cells, n, replace=False)
            jobs.append((("agree", n, r), s, None))

    for r in range(a.n_boot):
        s = rng.choice(cells, len(cells), replace=True)
        jobs.append((("boot", r), s, r))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturb-dir", required=True)
    ap.add_argument("--transitions", nargs="*", required=True)
    ap.add_argument("--tag", default="parallel")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--n-splithalf", type=int, default=10)
    ap.add_argument("--n-boot", type=int, default=10)
    ap.add_argument("--sat-reps", type=int, default=3)
    ap.add_argument("--sat-ns", type=int, nargs="*", default=[100, 200, 300, 500, 750])
    ap.add_argument("--agree-ns", type=int, nargs="*",
                    default=[100, 200, 300, 500, 750, 1000, 1250])
    a = ap.parse_args()
    if not os.path.isabs(a.perturb_dir):
        sys.exit("--perturb-dir must be absolute (stats_for uses symlinks; relative paths do not resolve under /tmp)")
    os.makedirs(OUT, exist_ok=True)

    rows_ref, rows_sh, rows_sat, rows_boot = [], [], [], []
    for t in a.transitions:
        src = os.path.join(a.perturb_dir, t)
        cf = cell_files(src)
        cells = np.array(sorted(cf))
        jobs = build_jobs(cells, a)
        print(f"\n=== {t}: {len(cells)} cells, {len(jobs)} stats rounds, {a.jobs}-way parallel ===",
              flush=True)

        cache = os.path.join(OUT, f"_cache_{a.tag}_{t}")
        os.makedirs(cache, exist_ok=True)
        n_cached = len([1 for j in jobs if os.path.exists(
            os.path.join(cache, "_".join(map(str, j[0])) + ".pkl"))])
        if n_cached:
            print(f"Cache hits {n_cached}/{len(jobs)} rounds; will skip", flush=True)

        res, t0, done = {}, time.time(), 0
        # Must use spawn: geneformer imports torch, which starts threads at import time;
        # forking a multi-threaded parent leaves inconsistent lock state — default fork
        # with 28 workers previously failed (SIGBUS, trap stack segment). spawn restarts each worker cleanly;
        # cost is one re-import per worker (~tens of seconds), paid once.
        with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init,
                                 initargs=(cf, cache),
                                 mp_context=mp.get_context("spawn")) as ex:
            futs = [ex.submit(_run, j) for j in jobs]
            for f in as_completed(futs):
                k, d = f.result()
                res[k] = d
                done += 1
                el = time.time() - t0
                # ETA must account for parallelism. elapsed/done*remaining is the serial formula
                # and overestimates by >10× before the first batch finishes (e.g. reported 153 min, actual 15).
                # Use remaining_batches × time_per_batch: time_per_batch ≈ elapsed / completed_batches.
                waves_done = max(done / a.jobs, 1e-9)
                waves_left = (len(jobs) - done) / a.jobs
                eta = el / waves_done * waves_left
                print(f"  [{done}/{len(jobs)}] {'/'.join(map(str, k))}  "
                      f"Passed filter {len(d)} genes  elapsed {el/60:.1f} min"
                      f"ETA remaining {eta/60:.1f} min", flush=True)

        ref = res[("ref",)]
        ref_top = list(ref.Gene_name.head(TOPN))
        rows_ref.append({"transition": t, "n_cells": len(cells), "n_pass_filter": len(ref),
                         "seconds": round(time.time() - t0, 1)})

        for r in range(a.n_splithalf):
            A, B = res[("sh", r, "A")], res[("sh", r, "B")]
            rho, nov = rank_corr(A, B)
            rows_sh.append({"transition": t, "rep": r, "n_per_half": len(cells) // 2,
                            "jaccard_top200": jac(A.Gene_name.head(TOPN), B.Gene_name.head(TOPN)),
                            "spearman_shift": rho, "n_common_genes": nov,
                            "n_pass_A": len(A), "n_pass_B": len(B)})

        for n in a.sat_ns:
            if 2 * n > len(cells):
                continue
            for r in range(a.sat_reps):
                A, B = res[("sat", n, r, "A")], res[("sat", n, r, "B")]
                rows_sat.append({"transition": t, "n": n, "rep": r,
                                 "kind": "disjoint_consistency",
                                 "jaccard_top200": jac(A.Gene_name.head(TOPN),
                                                       B.Gene_name.head(TOPN)),
                                 "spearman_shift": rank_corr(A, B)[0]})
        for n in a.agree_ns:
            if n >= len(cells):
                continue
            for r in range(a.sat_reps):
                A = res[("agree", n, r)]
                rows_sat.append({"transition": t, "n": n, "rep": r,
                                 "kind": "agreement_with_full",
                                 "jaccard_top200": jac(A.Gene_name.head(TOPN), ref_top),
                                 "spearman_shift": rank_corr(A, ref)[0]})

        for r in range(a.n_boot):
            A = res[("boot", r)]
            rows_boot.append({"transition": t, "rep": r,
                              "jaccard_top200": jac(A.Gene_name.head(TOPN), ref_top),
                              "spearman_shift": rank_corr(A, ref)[0], "n_pass": len(A)})

    for name, rows in [("reference", rows_ref), ("splithalf", rows_sh),
                       ("saturation", rows_sat), ("bootstrap", rows_boot)]:
        pd.DataFrame(rows).to_csv(f"{OUT}/{a.tag}_{name}.csv", index=False)

    sh, sat, bo = pd.DataFrame(rows_sh), pd.DataFrame(rows_sat), pd.DataFrame(rows_boot)
    print("\n" + "=" * 70)
    print("Rank stability summary (top-200 Jaccard)")
    print("=" * 70)
    for t in a.transitions:
        s, b = sh[sh.transition == t], bo[bo.transition == t]
        n = rows_ref[[x["transition"] for x in rows_ref].index(t)]["n_cells"]
        print(f"\n{t} ({n} cells)")
        print(f"  split-half（{n//2} vs {n//2}）: J = {s.jaccard_top200.mean():.3f} "
              f"± {s.jaccard_top200.std():.3f}   Spearman = {s.spearman_shift.mean():.3f}")
        print(f"bootstrap (vs full set)        : J = {b.jaccard_top200.mean():.3f}"
              f"± {b.jaccard_top200.std():.3f}")
        c = sat[(sat.transition == t) & (sat.kind == "disjoint_consistency")]
        if len(c):
            print("Consistency between non-overlapping subsets (extrapolable; whether n cells suffice):")
            for n_, g in c.groupby("n"):
                print(f"    n={n_:>5}: J = {g.jaccard_top200.mean():.3f} "
                      f"± {g.jaccard_top200.std():.3f}   "
                      f"Spearman = {g.spearman_shift.mean():.3f}")
        g2 = sat[(sat.transition == t) & (sat.kind == "agreement_with_full")]
        if len(g2):
            print("Agreement of subset vs full result (subset fraction only; not alone proof of convergence):")
            for n_, g in g2.groupby("n"):
                print(f"    n={n_:>5}: J = {g.jaccard_top200.mean():.3f}")
    print(f"\nWriting {OUT}/{a.tag}_*.csv", flush=True)


if __name__ == "__main__":
    main()
