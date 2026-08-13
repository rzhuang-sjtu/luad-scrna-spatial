#!/usr/bin/env python
"""H10 — Apply H1 sequencing-depth adjustment on the validation cohort (Takano 2024).

H1 covered only the discovery cohort E-MTAB-13530 (SEC61G Δ=+0.354 lost after adjustment);
the main-text spatial endpoint claims enrichment in both cohorts (Takano Δ=+0.532). Withdrawing from one cohort alone
leaves no answer if asked about the other; re-run the validation cohort under the same specification.

The validation objects lack precomputed MP scores and cell2location abundances (those live in misty_data/
per-section CSVs), so rebuild ROI: z(MP3)>0.5 and z(Malignant)>0.5, matching Fig. 8.
Per-section dummy variables, as in H1.

Output: results/spatial_depth_validation/
"""
import argparse
import glob
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.api as sm

BASE = "${DATA_ROOT}/ST/results/step09_okamura_validation"
OUT = "${PROJECT_ROOT}/results/spatial_depth_validation"
TARGETS = ["SEC61G", "SRSF9", "ANGPTL4"]
POS = ["KRT19", "MUC1", "KRT7", "KRT18", "EPCAM", "NAPSA"]
NEG = ["COL1A1", "PTPRC", "CD3D"]


def load():
    secs = sorted({os.path.basename(f).replace("_mp.csv", "")
                   for f in glob.glob(f"{BASE}/misty_data/*_mp.csv")})
    print(f"Validation cohort: {len(secs)} sections: {secs}", flush=True)
    Xs, obs, names = [], [], None
    for s in secs:
        h5 = f"{BASE}/section_h5ad/{s}.h5ad"
        if not os.path.exists(h5):
            print(f"skip {s} (no h5ad)"); continue
        a = sc.read_h5ad(h5)
        mp = pd.read_csv(f"{BASE}/misty_data/{s}_mp.csv", index_col=0)
        it = pd.read_csv(f"{BASE}/misty_data/{s}_intra.csv", index_col=0)
        # CSV barcodes carry a section suffix (AAAC...-1-LUAD_No_1); h5ad barcodes do not — align first
        mp.index = mp.index.str.replace(f"-{s}$", "", regex=True)
        it.index = it.index.str.replace(f"-{s}$", "", regex=True)
        keep = a.obs_names.intersection(mp.index).intersection(it.index)
        assert len(keep) > 0.9 * a.n_obs, f"{s} barcode mismatch: {len(keep)}/{a.n_obs}"
        a = a[keep].copy()
        tot = np.asarray(a.X.sum(1)).ravel()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        o = pd.DataFrame(dict(
            section=s, total_counts=tot,
            n_genes=np.asarray((a.X > 0).sum(1)).ravel().astype(float),
            MP3=mp.loc[keep, "MP3_score"].values if "MP3_score" in mp.columns
            else mp.loc[keep, "MP3"].values,
            Malignant=it.loc[keep, "Malignant"].values), index=keep)
        Xs.append(sp.csr_matrix(a.X))
        obs.append(o)
        names = a.var_names if names is None else names
        print(f"  {s}: {a.shape[0]} spot", flush=True)
    X = sp.vstack(Xs).tocsr()
    O = pd.concat(obs)
    print(f"total {X.shape[0]} spots × {X.shape[1]} genes", flush=True)
    return X, O, np.array(names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-random", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260805)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    X, O, genes = load()

    z = lambda v: (v - np.mean(v)) / np.std(v)
    roi = (z(O.MP3.values) > 0.5) & (z(O.Malignant.values) > 0.5)
    print(f"\nTumour-intrinsic ROI: {roi.sum()} / {len(roi)} spots"
          f"(main text reports 3,504)", flush=True)

    tot, ng = O.total_counts.values, O.n_genes.values
    print(f"Median depth in ROI {np.median(tot[roi]):.0f}, outside ROI {np.median(tot[~roi]):.0f},"
          f"ratio {np.median(tot[roi])/np.median(tot[~roi]):.2f}", flush=True)
    ltot, lgene = np.log1p(tot), np.log1p(ng)
    D = pd.get_dummies(O.section, drop_first=True).values.astype(float)

    def fit(gi):
        y = np.asarray(X[:, gi].todense()).ravel().astype(float)
        out = {"delta_M0": y[roi].mean() - y[~roi].mean()}
        for name, cols in [("M1", [ltot]), ("M2", [ltot, lgene])]:
            Z = sm.add_constant(np.column_stack([roi.astype(float)] + cols + [D]))
            r = sm.OLS(y, Z).fit()
            out[f"beta_{name}"] = r.params[1]
            out[f"p_{name}"] = r.pvalues[1]
        return out

    rows = []
    for g in TARGETS + POS + NEG:
        w = np.where(genes == g)[0]
        if not len(w):
            print(f"{g} not in matrix"); continue
        r = fit(int(w[0]))
        r["gene"] = g
        r["kind"] = ("target" if g in TARGETS else
                     "positive control" if g in POS else "negative control")
        rows.append(r)
    R = pd.DataFrame(rows)[["gene", "kind", "delta_M0", "beta_M1", "p_M1",
                            "beta_M2", "p_M2"]]
    R.to_csv(f"{OUT}/targets_and_controls.csv", index=False)
    print("\n=== Targets and controls (validation cohort) ===")
    print(R.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    det = np.asarray((X > 0).mean(0)).ravel()
    pool = np.where(det > 0.05)[0]
    rng = np.random.default_rng(a.seed)
    pick = rng.choice(pool, min(a.n_random, len(pool)), replace=False)
    rnd = []
    for i, gi in enumerate(pick, 1):
        rnd.append(fit(int(gi)))
        if i % 100 == 0:
            print(f"random genes {i}/{len(pick)}", flush=True)
    Q = pd.DataFrame(rnd)
    Q["gene"] = genes[pick]
    Q.to_csv(f"{OUT}/random_genes.csv", index=False)
    s = R[R.gene == "SEC61G"].iloc[0]
    print(f"\nSEC61G percentile among random genes:"
          f"before adjustment {(Q.delta_M0 < s.delta_M0).mean()*100:.1f}%,"
          f"after adjustment {(Q.beta_M2 < s.beta_M2).mean()*100:.1f}%")
    print(f"Fraction of random genes positive before adjustment {(Q.delta_M0 > 0).mean()*100:.1f}%,"
          f"after adjustment {(Q.beta_M2 > 0).mean()*100:.1f}%")


if __name__ == "__main__":
    main()
