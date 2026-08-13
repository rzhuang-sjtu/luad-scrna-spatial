"""
Fig 8 / S11 — Per-section Mann-Whitney significance for the 3 lead genes
(SEC61G / SRSF9 / ANGPTL4) inside vs outside the tumor-intrinsic ROI
(z(MP3)>0.5 AND z(Malignant)>0.5).

Inputs (per cohort):
  - cohort_with_roi.h5ad  : per-spot MP3_score, sample
  - all_sections_c2l.h5ad : Malignant abundance from c2l q05_w_sf

Outputs (in fig8/v2_500/data/):
  - per_section_lead_gene_sig.csv      : per-(cohort,sample,gene) MW test + FDR
  - per_cohort_lead_gene_sig.csv       : pooled per-(cohort,gene) MW test + FDR
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

DATASETS = {
    "E-MTAB-13530": {
        "roi": "${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad",
        "c2l": "${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad",
    },
    "Okamura": {
        "roi": "${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad",
        "c2l": "${DATA_ROOT}/ST/results/step09_okamura_validation/all_sections_c2l.h5ad",
    },
}
GENES = ["SEC61G", "SRSF9", "ANGPTL4"]
OUT   = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data")


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def stars(p):
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


def zscore(x):
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / (x.std(ddof=0) + 1e-12)


per_section_rows = []
per_cohort_rows  = []

for cohort, cfg in DATASETS.items():
    print(f"\n=== {cohort} ===")
    if not Path(cfg["roi"]).exists() or not Path(cfg["c2l"]).exists():
        print(f"  [SKIP] missing h5ad")
        continue
    ad_roi = sc.read_h5ad(cfg["roi"])
    ad_c2l = sc.read_h5ad(cfg["c2l"])
    common = ad_roi.obs_names.intersection(ad_c2l.obs_names)
    ad_roi = ad_roi[common].copy()
    ad_c2l = ad_c2l[common].copy()
    samples = ad_roi.obs["sample"].astype(str).values
    mp3 = ad_roi.obs["MP3_score"].astype(float).values
    mal = ad_c2l.obsm["q05_cell_abundance_w_sf"]["q05cell_abundance_w_sf_Malignant"].astype(float).values
    # global z-score (matches 07_tumor_intrinsic_roi.py convention)
    z_mp3 = zscore(mp3)
    z_mal = zscore(mal)
    new_roi = (z_mp3 > 0.5) & (z_mal > 0.5)
    print(f"  spots: {len(common)}  ROI={int(new_roi.sum())}  non-ROI={int((~new_roi).sum())}")

    # ----- pooled per-cohort test -----
    cohort_pvals = []
    cohort_idx   = []
    for g in GENES:
        if g not in ad_roi.var_names:
            print(f"  [SKIP] {g} not in var_names"); continue
        e = ad_roi[:, g].X
        e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
        a = e[new_roi]; b = e[~new_roi]
        if len(a) < 5 or len(b) < 5:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        per_cohort_rows.append({
            "cohort": cohort, "gene": g,
            "n_roi": int(len(a)), "n_nonroi": int(len(b)),
            "mean_roi": float(a.mean()), "mean_nonroi": float(b.mean()),
            "delta": float(a.mean() - b.mean()),
            "U_stat": float(u), "p_raw": float(p),
        })
        cohort_pvals.append(p); cohort_idx.append(g)

    # ----- per-section test -----
    for sec in np.unique(samples):
        m = samples == sec
        roi_s = new_roi[m]
        if roi_s.sum() < 3 or (~roi_s).sum() < 3:
            continue
        sec_buf = []; sec_p = []
        for g in GENES:
            if g not in ad_roi.var_names:
                continue
            e = ad_roi[:, g].X
            e = e.toarray().flatten() if hasattr(e, "toarray") else np.asarray(e).flatten()
            e_s = e[m]
            a = e_s[roi_s]; b = e_s[~roi_s]
            if len(a) < 3 or len(b) < 3:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            sec_buf.append({
                "cohort": cohort, "sample": sec, "gene": g,
                "n_roi": int(len(a)), "n_nonroi": int(len(b)),
                "mean_roi": float(a.mean()), "mean_nonroi": float(b.mean()),
                "delta": float(a.mean() - b.mean()),
                "U_stat": float(u), "p_raw": float(p),
            })
            sec_p.append(p)
        if sec_buf:
            q_sec = bh_fdr(np.asarray(sec_p))
            for row, q in zip(sec_buf, q_sec):
                row["p_fdr"] = float(q); row["sig"] = stars(q)
                per_section_rows.append(row)

# Per-cohort BH-FDR within cohort × 3 genes
if per_cohort_rows:
    cdf = pd.DataFrame(per_cohort_rows)
    out_pieces = []
    for c, g in cdf.groupby("cohort", sort=False):
        q = bh_fdr(g["p_raw"].to_numpy())
        tmp = g.copy(); tmp["p_fdr"] = q; tmp["sig"] = [stars(p) for p in q]
        out_pieces.append(tmp)
    cdf = pd.concat(out_pieces, axis=0, ignore_index=True)
    cdf.to_csv(OUT / "per_cohort_lead_gene_sig.csv", index=False)
    print(f"\nwrote {OUT / 'per_cohort_lead_gene_sig.csv'} ({len(cdf)} rows)")
    print(cdf[["cohort","gene","mean_roi","mean_nonroi","delta","p_fdr","sig"]].to_string(index=False))

if per_section_rows:
    sdf = pd.DataFrame(per_section_rows)
    sdf.to_csv(OUT / "per_section_lead_gene_sig.csv", index=False)
    print(f"\nwrote {OUT / 'per_section_lead_gene_sig.csv'} ({len(sdf)} rows)")
    n_sig = sdf.groupby(["cohort","gene"])["sig"].apply(
        lambda s: int((s != "ns").sum())).reset_index(name="n_sig_sections")
    n_tot = sdf.groupby(["cohort","gene"]).size().reset_index(name="n_sections")
    summary = n_sig.merge(n_tot, on=["cohort","gene"])
    print("\nper-cohort × gene consistency:")
    print(summary.to_string(index=False))
