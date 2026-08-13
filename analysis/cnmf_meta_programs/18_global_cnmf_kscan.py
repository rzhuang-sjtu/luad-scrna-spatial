"""Step 18: global K-selection curve for Fig S3A.

Subsample 15k malignant cells (stratified by dominant_MP), use 2000 HVG
(taken from existing HVG flag in h5ad), run sklearn NMF for K=8..20 with
N_RUNS=8 random seeds per K. Compute:
  - stability: mean pairwise matched-cosine of H components across seeds
  - error: mean Frobenius reconstruction error

Output: ${WORK_ROOT}/luad_figures/fig_s2/global_cnmf_k_selection.csv
  Columns: K, n_runs, stability, stability_std, error_mean, error_std
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
OUT = Path("${WORK_ROOT}/luad_figures/fig_s2")
OUT.mkdir(parents=True, exist_ok=True)

N_SUBSAMPLE = 15000
N_HVG = 2000
K_RANGE = list(range(8, 21))   # 8..20 inclusive
N_RUNS = 8                      # seeds per K
N_JOBS = 8                      # joblib parallel fits


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stability_matched_cosine(Hs):
    """Mean pairwise matched cosine similarity between H matrices.

    Hs: list of (K, n_genes) ndarrays. Uses Hungarian matching.
    """
    from scipy.optimize import linear_sum_assignment
    n = len(Hs)
    if n < 2: return np.nan
    # L2-normalize rows for cosine
    Hns = [H / (np.linalg.norm(H, axis=1, keepdims=True) + 1e-12) for H in Hs]
    sims = []
    for i in range(n):
        for j in range(i+1, n):
            cos = Hns[i] @ Hns[j].T   # (K, K)
            row, col = linear_sum_assignment(-cos)  # maximize → negate
            sims.append(cos[row, col].mean())
    return float(np.mean(sims)), float(np.std(sims))


def fit_one(X_dense, K, seed):
    from sklearn.decomposition import NMF
    nmf = NMF(n_components=K, init="random", solver="cd",
              beta_loss="frobenius", max_iter=300, tol=1e-3,
              random_state=seed)
    W = nmf.fit_transform(X_dense)
    return nmf.components_, float(nmf.reconstruction_err_)


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    # Subsample stratified by dominant_MP
    log(f"stratified subsample to {N_SUBSAMPLE}")
    rng = np.random.default_rng(0)
    mp = a.obs["dominant_MP"].astype(str).values
    keep = []
    for label in np.unique(mp):
        idx = np.where(mp == label)[0]
        n_take = max(50, int(N_SUBSAMPLE * len(idx) / len(mp)))
        if len(idx) <= n_take:
            keep.extend(idx)
        else:
            keep.extend(rng.choice(idx, size=n_take, replace=False))
    keep = np.array(sorted(keep))
    log(f"  retained {len(keep)} cells")

    # HVG selection: use existing flag if available, else compute fresh
    if "highly_variable" in a.var.columns and a.var["highly_variable"].sum() >= N_HVG:
        hv = a.var.sort_values("highly_variable", ascending=False).head(N_HVG).index.tolist()
        log(f"  using existing HVG flag → {len(hv)} genes")
    else:
        log(f"  computing HVG ({N_HVG})")
        sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat",
                                     batch_key="dataset", inplace=True)
        hv = a.var[a.var["highly_variable"]].index.tolist()

    sub = a[keep, hv].copy()
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    # Ensure non-negative for NMF (lognorm shouldn't have negatives, but clip)
    if X.min() < 0:
        log(f"  X min={X.min():.3f} → clipping to 0")
        X = np.clip(X, 0, None)
    X = X.astype(np.float32)
    log(f"  X for NMF: {X.shape}, range [{X.min():.3f}, {X.max():.3f}]")

    # Run K-scan
    log(f"K scan: K in {K_RANGE}, {N_RUNS} runs each, {N_JOBS} parallel jobs")
    from joblib import Parallel, delayed
    rows = []
    for K in K_RANGE:
        tk = time.time()
        results = Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
            delayed(fit_one)(X, K, seed) for seed in range(N_RUNS)
        )
        Hs = [r[0] for r in results]
        errs = [r[1] for r in results]
        stab_mean, stab_std = stability_matched_cosine(Hs)
        rows.append({
            "K": K, "n_runs": N_RUNS,
            "stability": stab_mean, "stability_std": stab_std,
            "error_mean": float(np.mean(errs)),
            "error_std": float(np.std(errs)),
            "elapsed_sec": time.time()-tk,
        })
        log(f"  K={K:2d}  stability={stab_mean:.4f}±{stab_std:.4f}  "
            f"error={np.mean(errs):.2f}  in {time.time()-tk:.1f}s")

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT/"global_cnmf_k_selection.csv", index=False)
    log(f"\nresult:\n{out_df.round(4).to_string(index=False)}")
    # Also pick K_optimal
    out_df["delta_stab"] = out_df["stability"].diff()
    log(f"\nK with peak stability: {int(out_df.loc[out_df['stability'].idxmax(), 'K'])}")
    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
