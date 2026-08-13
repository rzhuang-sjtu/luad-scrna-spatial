#!/usr/bin/env python
"""
H4b — Whether MP3 rises with pseudotime: multi-method recheck (revision analysis).

H4 used pooled Spearman across all cells and found a slight MP3 decline along pseudotime (ρ=−0.076).
That approach has a clear flaw: **both MP3 and pseudotime vary across patients, so pooled correlation can mask**
**within-patient relationships** (Simpson's paradox). The MP3–SEC61G analysis failed the same way:
pooled analysis was null; significance appeared only after patient fixed effects.

Recompute with six methods; report all, no cherry-picking:

  M1 pooled Spearman            H4 baseline
  M2 per-patient Spearman + sign test / Wilcoxon   within each patient, then test direction consistency
  M3 patient fixed-effect OLS          MP3 ~ pseudotime + patient dummies
  M4 mixed effects (patient random intercept)
  M5 binned curve shape              not just monotonicity; actual shape (may rise then fall)
  M6 within-branch analysis                opposing PAGA/DPT branches can cancel when pooled

Criterion: if all six methods agree in direction, the conclusion is credible; if only the pooled analysis is negative while within-patient is positive,
the H4 conclusion is a Simpson artefact and must be overturned.

Usage: python H4b_mp3_pseudotime_multitest.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats as sps

H5 = "${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad"
PT = "${PROJECT_ROOT}/results/dpt_root_sensitivity/pseudotime_by_root.csv.gz"
OUT = "${PROJECT_ROOT}/results/dpt_root_sensitivity"
MIN_CELLS = 50          # minimum cells for within-patient analysis


def main():
    ad = sc.read_h5ad(H5, backed="r")
    obs = ad.obs
    P = pd.read_csv(PT)
    ref = P.columns[0]
    d = pd.DataFrame({
        "pt": P[ref].values,
        "MP3": obs["MP3_score"].values.astype(float),
        "MP4": obs["MP4_score"].values.astype(float),
        "patient": obs["patient_id"].astype(str).values,
        "dataset": obs["dataset"].astype(str).values,
    })
    if "dpt_groups" in obs.columns:
        d["branch"] = obs["dpt_groups"].astype(str).values
    d = d[d.pt.notna()].copy()
    print(f"cells {len(d):,}, patients {d.patient.nunique()}, root '{ref}'", flush=True)

    r1 = sps.spearmanr(d.pt, d.MP3)
    r1b = sps.spearmanr(d.pt, d.MP4)
    print(f"\n【M1 pooled Spearman (H4 baseline)】")
    print(f"  MP3 ρ={r1.statistic:+.4f}  p={r1.pvalue:.3g}")
    print(f"MP4 ρ={r1b.statistic:+.4f}  p={r1b.pvalue:.3g} (control: root is the max-MP4 cell, should be negative)")

    rows = []
    for p, g in d.groupby("patient"):
        if len(g) < MIN_CELLS or g.pt.nunique() < 10:
            continue
        rows.append(dict(patient=p, n=len(g),
                         rho_MP3=sps.spearmanr(g.pt, g.MP3).statistic,
                         rho_MP4=sps.spearmanr(g.pt, g.MP4).statistic))
    W = pd.DataFrame(rows)
    W.to_csv(f"{OUT}/per_patient_rho.csv", index=False)
    v = W.rho_MP3.dropna()
    npos = int((v > 0).sum())
    sign_p = sps.binomtest(npos, len(v), 0.5).pvalue
    wil = sps.wilcoxon(v)
    print(f"\n【M2 per-patient Spearman】{len(W)} patients (≥{MIN_CELLS} cells)")
    print(f"MP3 ρ median {v.median():+.4f}   patients positive {npos}/{len(v)}"
          f"（{100*npos/len(v):.1f}%）")
    print(f"sign test p={sign_p:.4g}   Wilcoxon p={wil.pvalue:.4g}")
    print(f"control MP4 ρ median {W.rho_MP4.median():+.4f},"
          f"positive {int((W.rho_MP4>0).sum())}/{len(W)}")

    sub = d.copy()
    sub["zpt"] = (sub.pt - sub.pt.mean()) / sub.pt.std()
    m3 = smf.ols("MP3 ~ zpt + C(patient)", data=sub).fit()
    print(f"\n【M3 patient fixed-effect OLS】MP3 ~ pseudotime + patient")
    print(f"  β={m3.params['zpt']:+.5f}  p={m3.pvalues['zpt']:.3g}  "
          f"95%CI [{m3.conf_int().loc['zpt',0]:+.5f}, {m3.conf_int().loc['zpt',1]:+.5f}]")

    try:
        m4 = smf.mixedlm("MP3 ~ zpt", sub, groups=sub.patient).fit(method="lbfgs")
        print(f"\n【M4 mixed effects (patient random intercept)】")
        print(f"  β={m4.params['zpt']:+.5f}  p={m4.pvalues['zpt']:.3g}")
    except Exception as e:
        print(f"\n【M4 mixed effects】did not converge: {e}")

    sub["bin"] = pd.qcut(sub.pt, 20, labels=False, duplicates="drop")
    curve = sub.groupby("bin").agg(MP3=("MP3", "mean"), MP4=("MP4", "mean"),
                                   n=("MP3", "size"))
    curve.to_csv(f"{OUT}/binned_curve.csv")
    print(f"\n【M5 ventile curve shape】(not just monotonicity; actual shape)")
    print(f"MP3: first {curve.MP3.iloc[0]:.4f} → peak {curve.MP3.max():.4f}"
          f"(ventile {int(curve.MP3.idxmax())+1}) → last {curve.MP3.iloc[-1]:.4f}")
    print(f"MP4: first {curve.MP4.iloc[0]:.4f} → last {curve.MP4.iloc[-1]:.4f}")
    up = int((curve.MP3.diff().dropna() > 0).sum())
    print(f"adjacent ventiles rising {up}/{len(curve)-1} times")
    # Linear trend test (weighted regression of bin means)
    x = np.arange(len(curve))
    sl, ic, rr, pp, se = sps.linregress(x, curve.MP3.values)
    print(f"linear trend of bin means: slope {sl:+.5f}  p={pp:.3g}  R²={rr**2:.3f}")

    if "branch" in d.columns:
        print(f"\n【M6 within-branch analysis】(opposing branches cancel when pooled)")
        for b, g in d.groupby("branch"):
            if len(g) < 200:
                continue
            r = sps.spearmanr(g.pt, g.MP3)
            print(f"branch {b:<6} n={len(g):>6,}  MP3 ρ={r.statistic:+.4f}  p={r.pvalue:.3g}")
    else:
        print("\n【M6】no dpt_groups branch column in object; skip")

    print("\n" + "=" * 60)
    print("Criterion: whether the six methods agree in direction")
    print("=" * 60)
    print(f"M1 pooled Spearman        {r1.statistic:+.4f}")
    print(f"M2 per-patient ρ median    {v.median():+.4f} ({100*npos/len(v):.0f}% patients positive)")
    print(f"M3 patient FE β            {m3.params['zpt']:+.5f}")
    print(f"M5 linear slope of bin means {sl:+.5f}")
    print(f"\nWrote {OUT}/", flush=True)


if __name__ == "__main__":
    main()
