"""Step 22b (Fig 3C v3): GeneSwitches with Nagelkerke R² + relaxed binarization.

Final tuning:
  - Binarization at log1p > 0.5 (vs v2's 1.0; many genes had only 5% ON before)
  - pct_on filter relaxed to [5%, 95%]
  - All genes passing filter are evaluated (no top-N variance cut)
  - Nagelkerke pseudo-R² (bounded [0, 1], comparable to standard R²) instead
    of McFadden (which is conservatively scaled 0-0.4)
  - Keep both R² metrics for transparency

Output overwrites geneswitches_results.csv / geneswitches_top.csv in fig3/.
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

IN = Path.home()/"luad/data/processed/luad_malignant_scored.h5ad"
PT_CSV = Path.home()/"luad/results/step10b_pseudotime.csv.gz"
FIG3 = Path("${WORK_ROOT}/luad_figures/fig3")
PCT_LO, PCT_HI = 0.05, 0.95
LOG1P_THRESHOLD = 0.5
R2_MIN = 0.05  # Nagelkerke


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fit_one(y_int, x):
    """Returns slope, intercept, mcfadden_R2, nagelkerke_R2, switch_pt."""
    from sklearn.linear_model import LogisticRegression
    if y_int.std() == 0: return None
    lr = LogisticRegression(C=1e6, solver="lbfgs", max_iter=500)
    lr.fit(x.reshape(-1, 1), y_int)
    slope = float(lr.coef_[0, 0]); intercept = float(lr.intercept_[0])
    if abs(slope) < 1e-9: return None
    p = 1.0 / (1.0 + np.exp(-(slope * x + intercept)))
    p = np.clip(p, 1e-9, 1-1e-9)
    LL_full = float(np.sum(y_int * np.log(p) + (1 - y_int) * np.log(1 - p)))
    p0 = float(y_int.mean()); p0 = max(min(p0, 1-1e-9), 1e-9)
    LL_null = float(np.sum(y_int * np.log(p0) + (1 - y_int) * np.log(1 - p0)))
    n = len(y_int)
    # McFadden
    R2_mc = 1 - LL_full / LL_null if LL_null != 0 else 0.0
    # Cox-Snell: 1 - (L0/L1)^(2/n) = 1 - exp((LL0 - LL1) * 2/n)
    R2_cs = 1.0 - np.exp((LL_null - LL_full) * 2.0 / n)
    # Nagelkerke: R2_cs / (1 - exp(2*LL_null/n))
    denom = 1.0 - np.exp(2.0 * LL_null / n)
    R2_nag = R2_cs / denom if denom > 0 else R2_cs
    switch_pt = -intercept / slope
    return slope, intercept, R2_mc, R2_nag, switch_pt


def main():
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  shape={a.shape}")

    pt = pd.read_csv(PT_CSV)[["barcode","pseudotime"]]
    pt_lookup = dict(zip(pt["barcode"], pt["pseudotime"]))
    pseudotime = np.array([pt_lookup.get(bc, np.nan) for bc in a.obs.index])
    keep = np.isfinite(pseudotime)
    a = a[keep].copy(); pseudotime = pseudotime[keep]
    log(f"  cells with pt: {len(a)}")

    pt_rank = pd.Series(pseudotime).rank(pct=True).values.astype(np.float64)

    if hasattr(a.X, "toarray"):
        X = a.X.toarray()
    else:
        X = np.asarray(a.X)
    if X.max() > 30:
        X = np.log1p(X.astype(np.float32))
    log(f"  X shape: {X.shape}, range [{X.min():.2f}, {X.max():.2f}]")

    Y = (X > LOG1P_THRESHOLD).astype(np.int8)
    pct_on = Y.mean(axis=0)
    log(f"  binarization at log1p > {LOG1P_THRESHOLD}; pct_on percentiles: "
        f"q25={np.percentile(pct_on,25):.3f} "
        f"q50={np.percentile(pct_on,50):.3f} "
        f"q75={np.percentile(pct_on,75):.3f}")

    keep_pct = (pct_on >= PCT_LO) & (pct_on <= PCT_HI)
    log(f"  {keep_pct.sum()}/{len(keep_pct)} genes pass pct filter [{PCT_LO},{PCT_HI}]")

    Y_sub = Y[:, keep_pct]
    X_sub = X[:, keep_pct]
    genes = a.var_names[keep_pct].tolist()

    log(f"fitting logistic regressions (n={len(genes)} genes)")
    rows = []
    for i, gene in enumerate(genes):
        result = fit_one(Y_sub[:, i].astype(np.int32), pt_rank)
        if result is None: continue
        slope, intercept, R2_mc, R2_nag, switch_pt = result
        if not (0 <= switch_pt <= 1): continue
        if R2_nag < R2_MIN: continue
        rows.append({
            "gene": gene,
            "switch_pseudotime_rank": switch_pt,
            "direction": "switch_up" if slope > 0 else "switch_down",
            "slope": slope,
            "intercept": intercept,
            "mcfadden_R2": R2_mc,
            "nagelkerke_R2": R2_nag,
            "pct_expressing": float((X_sub[:, i] > 0).mean()),
            "pct_on_thresholded": float(Y_sub[:, i].mean()),
        })
        if (i+1) % 1000 == 0:
            log(f"  {i+1}/{len(genes)}; passed Nagelkerke R²>{R2_MIN}: {len(rows)}")

    df = pd.DataFrame(rows).sort_values("switch_pseudotime_rank").reset_index(drop=True)
    df.to_csv(FIG3/"geneswitches_results.csv", index=False)
    log(f"  geneswitches_results.csv ({df.shape})")
    log(f"  Nagelkerke R²: max={df['nagelkerke_R2'].max():.3f}  "
        f"q90={df['nagelkerke_R2'].quantile(0.9):.3f}  "
        f"q50={df['nagelkerke_R2'].quantile(0.5):.3f}")
    log(f"  McFadden  R²: max={df['mcfadden_R2'].max():.3f}  "
        f"q50={df['mcfadden_R2'].quantile(0.5):.3f}")
    log(f"  Nagelkerke > 0.5: {(df['nagelkerke_R2'] > 0.5).sum()}; "
        f"> 0.3: {(df['nagelkerke_R2'] > 0.3).sum()}; "
        f"> 0.1: {(df['nagelkerke_R2'] > 0.1).sum()}")

    top_up = df[df["direction"]=="switch_up"].sort_values("nagelkerke_R2",
                                                            ascending=False).head(100)
    top_dn = df[df["direction"]=="switch_down"].sort_values("nagelkerke_R2",
                                                              ascending=False).head(100)
    top = pd.concat([top_up, top_dn], ignore_index=True)
    top.to_csv(FIG3/"geneswitches_top.csv", index=False)
    log(f"  geneswitches_top.csv ({top.shape})")
    log("\nTop-15 switch-UP (highest Nagelkerke R²):")
    log(top_up.head(15)[["gene","switch_pseudotime_rank","nagelkerke_R2","mcfadden_R2"]]
        .round(3).to_string(index=False))
    log("\nTop-15 switch-DOWN (highest Nagelkerke R²):")
    log(top_dn.head(15)[["gene","switch_pseudotime_rank","nagelkerke_R2","mcfadden_R2"]]
        .round(3).to_string(index=False))

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
