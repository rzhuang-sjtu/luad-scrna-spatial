#!/usr/bin/env python
"""H11 — Correct statistics for adjustment for spot sequencing depth: estimate per section, then test across sections.

Why not H1/H10: those scripts treat >20,000 spots as independent OLS observations,
but adjacent Visium spots are strongly autocorrelated, so effective n is far below spot count —
P values such as 1.2e-40 are inflated by autocorrelation and cannot support cross-cohort comparison.

Treat sections as independent units (discovery 12, validation 8):
  1. Fit y ~ ROI + log(total_counts) + log(n_genes) within each section and take β;
  2. Across sections, one-sample Wilcoxon and t-tests, and report the count of sections with β>0;
  3. Run random genes the same way per section; rank SEC61G by section-mean β
     — a relative measure that cancels overall inflation shared by both cohorts.
ROI definition matches the main text (within-cohort z(MP3)>0.5 and z(Malignant)>0.5); only the statistical unit changes.

Output: results/spatial_depth_per_section/
"""
import glob
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.api as sm
from scipy import stats

OUT = "${PROJECT_ROOT}/results/spatial_depth_per_section"
DISC = "${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad"
C2L = "${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad"
VAL = "${DATA_ROOT}/ST/results/step09_okamura_validation"
TARGETS = ["SEC61G", "SRSF9", "ANGPTL4"]
POS = ["KRT19", "MUC1", "KRT7", "KRT18", "EPCAM", "NAPSA"]
NEG = ["COL1A1", "PTPRC", "CD3D"]
N_RANDOM = 200
MIN_ROI = 30


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std()


def load_discovery():
    ad = sc.read_h5ad(DISC)
    c2l = sc.read_h5ad(C2L)
    common = ad.obs_names.intersection(c2l.obs_names)
    ad, c2l = ad[common].copy(), c2l[common].copy()
    mal = c2l.obsm["q05_cell_abundance_w_sf"][
        "q05cell_abundance_w_sf_Malignant"].values.astype(float)
    roi = (z(ad.obs["MP3_score"].values.astype(float)) > 0.5) & (z(mal) > 0.5)
    X = sp.csr_matrix(ad.X)
    O = pd.DataFrame(dict(section=ad.obs["sample"].astype(str).values,
                          roi=roi,
                          total=np.asarray(ad.obs["total_counts"], float),
                          ngene=np.asarray(ad.obs["n_genes_by_counts"], float)))
    return X, O, np.array(ad.var_names)


def load_validation():
    secs = sorted({os.path.basename(f).replace("_mp.csv", "")
                   for f in glob.glob(f"{VAL}/misty_data/*_mp.csv")})
    Xs, obs, names = [], [], None
    for s in secs:
        h5 = f"{VAL}/section_h5ad/{s}.h5ad"
        if not os.path.exists(h5):
            continue
        a = sc.read_h5ad(h5)
        mp = pd.read_csv(f"{VAL}/misty_data/{s}_mp.csv", index_col=0)
        it = pd.read_csv(f"{VAL}/misty_data/{s}_intra.csv", index_col=0)
        mp.index = mp.index.str.replace(f"-{s}$", "", regex=True)
        it.index = it.index.str.replace(f"-{s}$", "", regex=True)
        keep = a.obs_names.intersection(mp.index).intersection(it.index)
        a = a[keep].copy()
        tot = np.asarray(a.X.sum(1)).ravel()
        ng = np.asarray((a.X > 0).sum(1)).ravel().astype(float)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        Xs.append(sp.csr_matrix(a.X))
        obs.append(pd.DataFrame(dict(section=s, total=tot, ngene=ng,
                                     MP3=mp.loc[keep, "MP3_score"].values,
                                     Mal=it.loc[keep, "Malignant"].values)))
        names = a.var_names if names is None else names
    X = sp.vstack(Xs).tocsr()
    O = pd.concat(obs, ignore_index=True)
    O["roi"] = (z(O.MP3.values) > 0.5) & (z(O.Mal.values) > 0.5)
    return X, O, np.array(names)


def per_section(X, O, gi):
    """Per-section β (depth-adjusted) and Δ (unadjusted). Drop sections with too few ROI spots."""
    y_all = np.asarray(X[:, gi].todense()).ravel().astype(float)
    out = []
    for s, idx in O.groupby("section").groups.items():
        idx = np.asarray(idx)
        r = O.roi.values[idx]
        if r.sum() < MIN_ROI or (~r).sum() < MIN_ROI:
            continue
        y = y_all[idx]
        d = y[r].mean() - y[~r].mean()
        Z = sm.add_constant(np.column_stack(
            [r.astype(float), np.log1p(O.total.values[idx]),
             np.log1p(O.ngene.values[idx])]))
        b = sm.OLS(y, Z).fit().params[1]
        out.append((s, d, b, int(r.sum()), len(idx)))
    return pd.DataFrame(out, columns=["section", "delta", "beta",
                                      "n_roi", "n_spot"])


def summarise(tag, X, O, genes, res):
    print(f"\n{'═'*66}\n{tag}\n{'═'*66}", flush=True)
    print(f"sections {O.section.nunique()}, spots {len(O)},"
          f"ROI {int(O.roi.sum())}", flush=True)
    dep = O.groupby("section").apply(
        lambda g: np.median(g.total[g.roi]) / np.median(g.total[~g.roi]),
        include_groups=False)
    print(f"median per-section ROI/non-ROI depth ratio {dep.median():.2f}"
          f"（{dep.min():.2f}–{dep.max():.2f}）", flush=True)

    rows = []
    for g in TARGETS + POS + NEG:
        w = np.where(genes == g)[0]
        if not len(w):
            print(f"{g} not in matrix"); continue
        t = per_section(X, O, int(w[0]))
        wil = stats.wilcoxon(t.beta)[1] if len(t) >= 6 else np.nan
        tt = stats.ttest_1samp(t.beta, 0)
        rows.append(dict(gene=g,
                         kind=("target" if g in TARGETS else
                               "positive control" if g in POS else "negative control"),
                         n_sec=len(t), delta_med=t.delta.median(),
                         beta_med=t.beta.median(),
                         beta_pos=f"{(t.beta > 0).sum()}/{len(t)}",
                         p_wilcoxon=wil, p_t=tt.pvalue))
        res[(tag, g)] = t
    R = pd.DataFrame(rows)
    print("\nPer-section estimates, cross-section tests (section = independent unit)")
    print(R.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # Random genes: section-median of per-section β, used to percentile-rank SEC61G
    det = np.asarray((X > 0).mean(0)).ravel()
    pool = np.where(det > 0.05)[0]
    rng = np.random.default_rng(20260805)
    pick = rng.choice(pool, min(N_RANDOM, len(pool)), replace=False)
    med_b, med_d = [], []
    for i, gi in enumerate(pick, 1):
        t = per_section(X, O, int(gi))
        med_b.append(t.beta.median())
        med_d.append(t.delta.median())
        if i % 50 == 0:
            print(f"random genes {i}/{len(pick)}", flush=True)
    Q = pd.DataFrame(dict(gene=genes[pick], beta_med=med_b, delta_med=med_d))
    s = R[R.gene == "SEC61G"].iloc[0]
    pct_b = 100 * (Q.beta_med < s.beta_med).mean()
    pct_d = 100 * (Q.delta_med < s.delta_med).mean()
    print(f"\nSEC61G percentile among {len(Q)} randomly expressed genes:"
          f"unadjusted {pct_d:.1f}%, adjusted {pct_b:.1f}%")
    print(f"Fraction of random genes with section-median β>0: {(Q.beta_med > 0).mean()*100:.1f}%"
          f"(unadjusted {(Q.delta_med > 0).mean()*100:.1f}%)")
    Q.to_csv(f"{OUT}/{tag}_random.csv", index=False)
    R.to_csv(f"{OUT}/{tag}_targets.csv", index=False)
    return R, Q, pct_b, pct_d


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    Xd, Od, gd = load_discovery()
    Rd, Qd, pbd, pdd = summarise("discovery", Xd, Od, gd, res)
    Xv, Ov, gv = load_validation()
    Rv, Qv, pbv, pdv = summarise("validation", Xv, Ov, gv, res)

    for tag in ("discovery", "validation"):
        res[(tag, "SEC61G")].assign(cohort=tag).to_csv(
            f"{OUT}/{tag}_SEC61G_by_section.csv", index=False)
    print("\n\n【Conclusion comparison】")
    for tag, R, pb, pd_ in (("discovery cohort", Rd, pbd, pdd),
                            ("validation cohort", Rv, pbv, pdv)):
        s = R[R.gene == "SEC61G"].iloc[0]
        print(f"{tag}: median Δ {s.delta_med:+.3f} → median β {s.beta_med:+.3f};"
              f"sections with β>0 {s.beta_pos}; Wilcoxon p={s.p_wilcoxon:.3g};"
              f"random-gene percentile {pd_:.0f}% → {pb:.0f}%")


if __name__ == "__main__":
    main()
