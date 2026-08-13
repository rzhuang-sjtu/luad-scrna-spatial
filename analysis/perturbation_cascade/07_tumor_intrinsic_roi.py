"""Complementary tumor-intrinsic ROI analysis.

Original ROI = z(NFkB)>0.5 AND z(Neutrophil_total)>0.5  →  stromal-immune niche
              (Malignant cells are excluded — Δ Malignant = -3.52 in ROI).
This biases against tumor-cell-intrinsic candidates like SEC61G.

New ROI = z(Malignant)>0.5 AND z(MP3_score)>0.5  →  EMT-active tumor core.
Compute per-spot z-scores within each dataset, then test SEC61G / SRSF9 / ANGPTL4
expression: new-ROI vs not-new-ROI.

Inputs:
  - {cohort}/cohort_with_roi.h5ad : per-spot MP3_score + gene expression
  - {cohort}/all_sections_c2l.h5ad : per-spot Malignant abundance (q05_w_sf)

Outputs to ${PROJECT_ROOT}/results/fig8_plot_data/v2_500/:
  - tumor_intrinsic_roi_long.csv : spot-level ROI flag + gene expr (for plotting)
  - tumor_intrinsic_roi_stats.csv : per-gene per-cohort summary
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/v2_500")
GENES = ["SEC61G", "SRSF9", "ANGPTL4"]

DATASETS = {
    "E-MTAB-13530": dict(
        roi="${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad",
        c2l="${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad",
    ),
    "Okamura": dict(
        roi="${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad",
        c2l="${DATA_ROOT}/ST/results/step09_okamura_validation/all_sections_c2l.h5ad",
    ),
}

def zscore(x):
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / (x.std(ddof=0) + 1e-12)

long_rows, stat_rows = [], []
for label, cfg in DATASETS.items():
    print(f"\n=== {label} ===")
    ad_roi = sc.read_h5ad(cfg["roi"])
    ad_c2l = sc.read_h5ad(cfg["c2l"])
    print(f"  cohort_with_roi: {ad_roi.shape}; c2l: {ad_c2l.shape}")
    if ad_roi.shape[0] != ad_c2l.shape[0]:
        print(f"  WARN: spot count mismatch — aligning by obs_names")
    # align by obs_names
    common = ad_roi.obs_names.intersection(ad_c2l.obs_names)
    print(f"  common spots: {len(common)}")
    ad_roi = ad_roi[common].copy()
    ad_c2l = ad_c2l[common].copy()
    mp3 = ad_roi.obs["MP3_score"].values.astype(float)
    mal = ad_c2l.obsm["q05_cell_abundance_w_sf"]["q05cell_abundance_w_sf_Malignant"].values.astype(float)
    z_mp3 = zscore(mp3)
    z_mal = zscore(mal)
    # original ROI from cohort_with_roi (just for reference)
    orig_roi = ad_roi.obs["roi"].values.astype(bool)
    new_roi = (z_mp3 > 0.5) & (z_mal > 0.5)
    print(f"  z(MP3)>0.5: {(z_mp3 > 0.5).sum()}  "
          f"z(Mal)>0.5: {(z_mal > 0.5).sum()}  "
          f"AND (new ROI): {new_roi.sum()}  "
          f"orig ROI: {orig_roi.sum()}  "
          f"overlap orig∩new: {(orig_roi & new_roi).sum()}")
    # gene expression
    for g in GENES:
        if g not in ad_roi.var_names:
            stat_rows.append({"gene": g, "dataset": label, "in_data": False}); continue
        e = ad_roi[:, g].X
        e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
        a = e[new_roi]; b = e[~new_roi]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        stat_rows.append({"gene": g, "dataset": label, "in_data": True,
                          "n_newROI": int(new_roi.sum()),
                          "n_nonROI": int((~new_roi).sum()),
                          "mean_newROI": float(a.mean()),
                          "mean_nonROI": float(b.mean()),
                          "delta_new_minus_non": float(a.mean() - b.mean()),
                          "mw_p": float(p),
                          "n_origROI_overlap": int((new_roi & orig_roi).sum())})
        # spot-level long table for plotting (subsample non-ROI for size)
        # keep all newROI spots; subsample nonROI to 5000 max for plotting size
        rng = np.random.default_rng(42)
        non_idx = np.where(~new_roi)[0]
        if len(non_idx) > 5000:
            non_idx = rng.choice(non_idx, size=5000, replace=False)
        roi_idx = np.where(new_roi)[0]
        keep_idx = np.concatenate([roi_idx, non_idx])
        long_rows.append(pd.DataFrame({
            "gene": g, "dataset": label,
            "expr": e[keep_idx],
            "new_roi": new_roi[keep_idx],
            "z_mp3": z_mp3[keep_idx],
            "z_mal": z_mal[keep_idx],
        }))

stat_df = pd.DataFrame(stat_rows)
stat_df.to_csv(OUT / "tumor_intrinsic_roi_stats.csv", index=False)
pd.concat(long_rows, ignore_index=True).to_csv(OUT / "tumor_intrinsic_roi_long.csv", index=False)

print("\n" + "=" * 80)
print(" Tumor-intrinsic ROI (z_Malignant>0.5 AND z_MP3>0.5) vs non-ROI")
print("=" * 80)
print(f"{'gene':<10} {'dataset':<14} {'n_ROI':>7} {'mean_ROI':>9} {'mean_non':>9} "
      f"{'Δ':>7} {'p':>10}")
for _, r in stat_df.iterrows():
    if not r.get("in_data", False):
        print(f"{r['gene']:<10} {r['dataset']:<14}  not in data"); continue
    print(f"{r['gene']:<10} {r['dataset']:<14} {r['n_newROI']:>7} "
          f"{r['mean_newROI']:>9.3f} {r['mean_nonROI']:>9.3f} "
          f"{r['delta_new_minus_non']:>+7.3f} {r['mw_p']:>10.2e}")
print(f"\nWrote: {OUT/'tumor_intrinsic_roi_stats.csv'}, {OUT/'tumor_intrinsic_roi_long.csv'}")
