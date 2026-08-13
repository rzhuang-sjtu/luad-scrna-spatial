#!/usr/bin/env python
"""
H3 — Recompute post-integration batch-mixing metrics (post-integration batch-mixing metrics).

The post-integration kBET rejection rate reported first time round (0.947) was high,
and came from an internal approximation rather than the original kBET definition.
The authors should discuss whether this degree of residual batch effect is typical
for an atlas of this scale and whether it could influence downstream biological conclusions.

**Fact that must be disclosed**: the reported 0.947 is not true kBET.
`scripts/export_fig_s1_data.py` notes «true kBET needs an R package; Shannon entropy /
chi-square approximation used here». This script reimplements kBET from the original definition and recomputes.

kBET definition (Büttner et al. 2019 Nat Methods):
  For each sampled cell, take k nearest neighbours; chi-square goodness-of-fit tests whether neighbour batch composition
  departs from global batch composition; rejection rate = fraction of cells with significant departure.
  Lower rejection = better mixing; 0 is perfect mixing.

**Important: global kBET is necessarily high in multi-dataset atlases and does not alone mean integration failed.**
Datasets differ in tissue origin and cell composition (some tumour-only, some with adjacent normal),
so the null «local batch composition = global» is biologically false.
What matters is mixing **within the same cell type** — this script therefore also reports stratified results.

Also report iLISI / cLISI as complements:
  high iLISI = good batch mixing; low cLISI = cell types not mixed away (no over-integration).

Usage: python H3_batch_integration_metrics.py [--n-sample 5000] [--k 50]
"""
import os, argparse
import numpy as np
import pandas as pd
import h5py
from scipy.stats import chi2 as chi2_dist
from sklearn.neighbors import NearestNeighbors

INT = "${PROJECT_ROOT}/data/processed/luad_integrated.h5ad"
CK = "${PROJECT_ROOT}/data/processed/luad_copykat.h5ad"   # correct cell-type labels
OUT = "${PROJECT_ROOT}/results/batch_metrics"


def load_cat(path, key):
    with h5py.File(path, "r") as h:
        g = h["obs"][key]
        if isinstance(g, h5py.Group):
            cats = np.array([x.decode() for x in g["categories"][:]])
            cd = g["codes"][:]
            return np.where(cd >= 0, cats[np.clip(cd, 0, None)], "NA")
        v = g[:]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


def kbet(emb, batch, k=50, n_sample=5000, seed=42):
    """kBET rejection rate: chi-square GOF of local vs global batch composition."""
    rng = np.random.default_rng(seed)
    n = len(emb)
    idx = rng.choice(n, min(n_sample, n), replace=False)
    nn = NearestNeighbors(n_neighbors=min(k, n - 1), n_jobs=-1).fit(emb)
    nb = nn.kneighbors(emb[idx], return_distance=False)

    ub, inv = np.unique(batch, return_inverse=True)
    gfreq = np.bincount(inv, minlength=len(ub)) / n
    keep = gfreq > 0
    gfreq = gfreq[keep]
    df = len(gfreq) - 1
    if df < 1:
        return np.nan, 0

    codes = inv[nb]                                   # (n_sample, k)
    obs = np.stack([(codes == j).sum(1) for j in np.where(keep)[0]], 1)
    exp = gfreq * codes.shape[1]
    stat = ((obs - exp) ** 2 / exp).sum(1)
    p = 1 - chi2_dist.cdf(stat, df)
    return float((p < 0.05).mean()), len(idx)


def lisi(emb, labels, perplexity=30, n_sample=5000, seed=42):
    """Simplified LISI: effective number of labels (1 = no mixing, L = fully mixed)."""
    rng = np.random.default_rng(seed)
    n = len(emb)
    idx = rng.choice(n, min(n_sample, n), replace=False)
    k = min(3 * perplexity, n - 1)
    nn = NearestNeighbors(n_neighbors=k, n_jobs=-1).fit(emb)
    dist, nb = nn.kneighbors(emb[idx])
    ub, inv = np.unique(labels, return_inverse=True)
    codes = inv[nb]
    # Gaussian-kernel weights (fixed bandwidth approximation; skip perplexity binary search)
    sigma = np.maximum(dist[:, perplexity - 1:perplexity], 1e-12)
    w = np.exp(-(dist / sigma) ** 2)
    w /= w.sum(1, keepdims=True)
    out = np.empty(len(idx))
    for i in range(len(idx)):
        p = np.bincount(codes[i], weights=w[i], minlength=len(ub))
        out[i] = 1.0 / np.sum(p ** 2)
    return float(np.median(out)), len(ub)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sample", type=int, default=5000)
    ap.add_argument("--k", type=int, default=50)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    print("Loading embeddings and labels ...", flush=True)
    with h5py.File(INT, "r") as h:
        pre = h["obsm"]["X_pca"][:]
        post = h["obsm"]["X_pca_harmony"][:]
    ds = load_cat(INT, "dataset")
    # Cell types from luad_copykat.h5ad labels — the luad_merged/integrated copy is incorrect
    # (58% of 240k cells labelled Fibroblast express PTPRC; see archive 11.4b)
    ct = load_cat(CK, "celltype_coarse")
    print(f"{pre.shape[0]:,} cells, {pre.shape[1]} dims;"
          f"{len(set(ds))} datasets, {len(set(ct))} cell types", flush=True)
    assert len(ct) == len(ds), "Cell counts differ between objects; align barcodes first"

    rows = []
    print(f"\n[Global kBET] k={a.k}, sample {a.n_sample}", flush=True)
    for lab, emb in [("Pre-integration X_pca", pre), ("Post-integration X_pca_harmony", post)]:
        r, n = kbet(emb, ds, a.k, a.n_sample)
        print(f"{lab:<22} rejection rate {r:.3f}", flush=True)
        rows.append(dict(层次="Global", 细胞类型="ALL", 嵌入=lab,
                         kbet拒绝率=r, n=n))

    print(f"\n[Per-cell-type kBET (post-integration)] — the informative readout", flush=True)
    print(f"Global kBET is high by construction in multi-dataset atlases: tissue origins differ,"
          f"so the null «local batch composition = global» does not hold.", flush=True)
    for c in sorted(set(ct)):
        m = ct == c
        if m.sum() < 500:
            continue
        r_pre, _ = kbet(pre[m], ds[m], a.k, a.n_sample)
        r_post, n = kbet(post[m], ds[m], a.k, a.n_sample)
        nb = len(set(ds[m]))
        print(f"{c:<14} n={int(m.sum()):>7,}  batches {nb}"
              f"pre {r_pre:.3f} → post {r_post:.3f}", flush=True)
        rows.append(dict(层次="By cell type", 细胞类型=c, 嵌入="X_pca_harmony",
                         kbet拒绝率=r_post, kbet整合前=r_pre,
                         n=int(m.sum()), n_batch=nb))
    pd.DataFrame(rows).to_csv(f"{OUT}/kbet.csv", index=False)

    print(f"\n[LISI] higher iLISI better (cap = n datasets); lower cLISI better (floor 1)",
          flush=True)
    li = []
    for lab, emb in [("Pre-integration X_pca", pre), ("Post-integration X_pca_harmony", post)]:
        i_, nb = lisi(emb, ds, n_sample=a.n_sample)
        c_, nc = lisi(emb, ct, n_sample=a.n_sample)
        print(f"  {lab:<22} iLISI {i_:.2f} / {nb}   cLISI {c_:.2f} / {nc}", flush=True)
        li.append(dict(嵌入=lab, iLISI=i_, iLISI上限=nb, cLISI=c_, cLISI上限=nc))
    pd.DataFrame(li).to_csv(f"{OUT}/lisi.csv", index=False)
    print(f"\nWriting {OUT}/", flush=True)


if __name__ == "__main__":
    main()
