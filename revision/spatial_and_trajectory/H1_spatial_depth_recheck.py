#!/usr/bin/env python
"""
H1 — Recheck whether the Fig. 8 spatial endpoint (SEC61G enrichment in the tumour-intrinsic ROI) is a sequencing-depth artefact.

Archive 3.1 already reports raw Δ=+0.354 → −0.068 after depth adjustment; recompute independently here
and add one quantity not in the archive that best settles the claim:

  **At what percentile of the random-gene distribution does raw SEC61G Δ fall?**
  If before adjustment it already sits near the median of random genes, it was never special,
  independent of whether the adjustment itself is correct — cleaner than debating the adjustment method.

Report three nested models so the favourable one is not singled out:
  M0  raw      Δ = mean in ROI − mean outside ROI
  M1  + log total counts
  M2  + log total counts + genes detected (most complete; also most prone to over-adjustment)

Over-adjustment controls (must be inspected together; otherwise true signal may have been wiped):
  Positive: KRT19 / MUC1 / KRT7 / KRT18 / EPCAM — epithelial markers; ROI is malignant, should stay positive
  Negative: COL1A1 / PTPRC / CD3D — immune/stromal markers; should stay negative

Usage: python H1_spatial_depth_recheck.py [--n-random 500]
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm

ROI_H5 = "${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad"
C2L_H5 = "${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad"
OUT = "${PROJECT_ROOT}/results/spatial_depth"
TARGETS = ["SEC61G", "SRSF9", "ANGPTL4"]
POS_CTRL = ["KRT19", "MUC1", "KRT7", "KRT18", "EPCAM", "NAPSA"]
NEG_CTRL = ["COL1A1", "PTPRC", "CD3D"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-random", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260805)
    a = ap.parse_args()
    import os
    os.makedirs(OUT, exist_ok=True)

    print("Loading spatial object ...", flush=True)
    ad = sc.read_h5ad(ROI_H5)
    print(f"  {ad.shape}", flush=True)
    print(f"obs columns: {[c for c in ad.obs.columns][:25]}", flush=True)

    # Critical: do not use the existing `roi` column — that is the Fig. 7 stromal-immune ROI
    # (z(NFkB)>0.5 and z(Neutrophil)>0.5, 1,400 spots; malignant cells excluded).
    # Fig. 8 tumour-intrinsic ROI is z(Malignant)>0.5 and z(MP3)>0.5; malignant abundance is from a separate file.
    print("Loading cell2location malignant abundance ...", flush=True)
    c2l = sc.read_h5ad(C2L_H5)
    common = ad.obs_names.intersection(c2l.obs_names)
    ad = ad[common].copy()
    c2l = c2l[common].copy()
    print(f"shared spots {len(common)}", flush=True)

    def z(v):
        v = np.asarray(v, float)
        return (v - v.mean()) / v.std()

    mp3 = ad.obs["MP3_score"].values.astype(float)
    mal = c2l.obsm["q05_cell_abundance_w_sf"][
        "q05cell_abundance_w_sf_Malignant"].values.astype(float)
    roi = (z(mp3) > 0.5) & (z(mal) > 0.5)
    orig = ad.obs["roi"].values.astype(bool)
    print(f"tumour-intrinsic ROI (z(MP3)>0.5 and z(Mal)>0.5): {roi.sum()} / {len(roi)} spots",
          flush=True)
    print(f"(control: Fig. 7 stromal-immune ROI {orig.sum()} spots; overlap {(roi & orig).sum()})",
          flush=True)

    X = ad.X
    if not isinstance(X, np.ndarray):
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    genes = np.array(ad.var_names)

    tot = np.asarray(ad.obs["total_counts"], float) if "total_counts" in ad.obs \
        else np.expm1(X).sum(1)
    ngene = np.asarray(ad.obs["n_genes_by_counts"], float) \
        if "n_genes_by_counts" in ad.obs else (X > 0).sum(1).astype(float)
    slice_ = ad.obs[[c for c in ad.obs.columns
                     if c.lower() in ("sample", "section", "slice", "library_id")][0]] \
        if any(c.lower() in ("sample", "section", "slice", "library_id")
               for c in ad.obs.columns) else None

    print(f"\nMedian depth in ROI {np.median(tot[roi]):.0f}, outside ROI {np.median(tot[~roi]):.0f},"
          f"ratio {np.median(tot[roi])/np.median(tot[~roi]):.2f}", flush=True)
    print(f"Median genes detected in ROI {np.median(ngene[roi]):.0f},"
          f"outside ROI {np.median(ngene[~roi]):.0f},"
          f"ratio {np.median(ngene[roi])/np.median(ngene[~roi]):.2f}", flush=True)

    ltot, lgene = np.log1p(tot), np.log1p(ngene)
    D = pd.get_dummies(slice_, drop_first=True).values.astype(float) \
        if slice_ is not None else None

    def fit(gi):
        y = X[:, gi].astype(float)
        m0 = y[roi].mean() - y[~roi].mean()
        out = {"delta_M0": m0}
        for name, cols in [("M1", [ltot]), ("M2", [ltot, lgene])]:
            Z = [roi.astype(float)] + cols
            Z = np.column_stack(Z + ([D] if D is not None else []))
            Z = sm.add_constant(Z)
            r = sm.OLS(y, Z).fit()
            out[f"beta_{name}"] = r.params[1]
            out[f"p_{name}"] = r.pvalues[1]
        return out

    rows = []
    for g in TARGETS + POS_CTRL + NEG_CTRL:
        w = np.where(genes == g)[0]
        if not len(w):
            print(f"{g} not in matrix"); continue
        r = fit(int(w[0])); r["gene"] = g
        r["kind"] = ("target" if g in TARGETS else
                     "positive control" if g in POS_CTRL else "negative control")
        rows.append(r)
    R = pd.DataFrame(rows)[["gene", "kind", "delta_M0", "beta_M1", "p_M1",
                            "beta_M2", "p_M2"]]
    R.to_csv(f"{OUT}/targets_and_controls.csv", index=False)
    print("\n=== Target genes and controls ===")
    print(R.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    rng = np.random.default_rng(a.seed)
    expressed = np.where((X > 0).mean(0) > 0.05)[0]      # detected in at least 5% of spots
    pick = rng.choice(expressed, min(a.n_random, len(expressed)), replace=False)
    rnd = []
    for i, gi in enumerate(pick, 1):
        rnd.append(fit(int(gi)))
        if i % 100 == 0:
            print(f"random genes {i}/{len(pick)}", flush=True)
    Q = pd.DataFrame(rnd)
    Q["gene"] = genes[pick]
    Q.to_csv(f"{OUT}/random_genes.csv", index=False)

    print(f"\n=== Distribution over {len(Q)} randomly expressed genes ===")
    for col, lab in [("delta_M0", "M0 raw Δ"), ("beta_M1", "M1 adjusted β"),
                     ("beta_M2", "M2 fully adjusted β")]:
        v = Q[col].values
        print(f"{lab:16s} median {np.median(v):+.4f}"
              f"fraction positive {100*(v>0).mean():.1f}%"
              f"2.5–97.5 percentiles [{np.percentile(v,2.5):+.4f}, {np.percentile(v,97.5):+.4f}]")

    print("\n=== Target-gene percentiles within the random distribution ===")
    for _, r in R[R.kind != "negative control"].iterrows():
        p0 = 100 * (Q.delta_M0.values < r.delta_M0).mean()
        p1 = 100 * (Q.beta_M1.values < r.beta_M1).mean()
        p2 = 100 * (Q.beta_M2.values < r.beta_M2).mean()
        print(f"  {r.gene:8s} M0 {p0:5.1f}%   M1 {p1:5.1f}%   M2 {p2:5.1f}%")
    print(f"\nWrote {OUT}/", flush=True)


if __name__ == "__main__":
    main()
