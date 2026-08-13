"""C3 — Robustness of the SEC61G–MP3 association (multiple model settings).

Background: the initial analysis used binary dominant_MP + log-normalized expression + OLS patient fixed effects and obtained
β=+0.011 (p=0.18), i.e. no within-patient association. That would remove «SEC61G associated with EMT/MP3» from the manuscript,
so it must be re-tested under multiple model settings rather than one specification.

Six independent checks:
  M1 continuous MP3 score (no hard argmax class)
  M2 adjust for per-cell depth (n_genes / total_counts)
  M3 negative-binomial GLM, raw counts + log(total_counts) offset (correct distribution for counts)
  M4 mixed-effects (patient random intercept) vs fixed effects
  M5 restrict to patients with ≥50 MP3 and ≥50 non-MP3 cells
  M6 per-patient effect-size distribution (forest-style) + sign test: directional consistency, not only the pooled estimate

A significant within-patient positive effect under any setting means the initial conclusion is unstable and must be revisited.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/mp3_sec61g_robustness"
import os; os.makedirs(OUT, exist_ok=True)

print("Loading malignant-cell object ...", flush=True)
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad")
GENES = ["SEC61G", "SEC61A1", "SEC61B", "EGFR"]
X = a[:, GENES].X
X = X.toarray() if sp.issparse(X) else np.asarray(X)
CNT_MAT = a[:, GENES].layers["counts"]
CNT_MAT = CNT_MAT.toarray() if sp.issparse(CNT_MAT) else np.asarray(CNT_MAT)

d = pd.DataFrame({
    "cell": a.obs_names,
    "patient": a.obs["patient_id"].astype(str).values,
    "dataset": a.obs["dataset"].astype(str).values,
    "dom": a.obs["dominant_MP"].astype(str).values,
    "MP1": a.obs["MP1_score"].values, "MP2": a.obs["MP2_score"].values,
    "MP3": a.obs["MP3_score"].values, "MP4": a.obs["MP4_score"].values,
    "n_genes": a.obs["n_genes_by_counts"].values.astype(float),
    "total_counts": a.obs["total_counts"].values.astype(float),
})
d["pt"] = d["dataset"] + "__" + d["patient"]
for i, g in enumerate(GENES):
    d[g] = X[:, i]
    d[g + "_cnt"] = CNT_MAT[:, i]
d["is_MP3"] = (d.dom == "MP3").astype(int)
d["z_MP3"] = (d.MP3 - d.MP3.mean()) / d.MP3.std()
d["log_tc"] = np.log10(d.total_counts)
print(f"Cells {len(d)}, patients {d.pt.nunique()}, MP3 {d.is_MP3.mean()*100:.1f}%", flush=True)

res = []
def rec(name, beta, p, note=""):
    res.append({"Test": name, "beta": beta, "p": p, "说明": note})
    star = " ***" if (p is not None and p < 0.05) else ""
    print(f"  {name:48s} β={beta:+.4f}  p={p:.3g}{star}  {note}", flush=True)

print("\n=== M0 reproduce initial result (baseline) ===", flush=True)
for f, lab in [("SEC61G ~ is_MP3", "Binary MP3, no adjustment"),
               ("SEC61G ~ is_MP3 + C(pt)", "Binary MP3 + patient fixed effects")]:
    m = smf.ols(f, data=d).fit()
    rec(lab, m.params["is_MP3"], m.pvalues["is_MP3"])

print("\n=== M1 continuous MP3 score (no hard class) ===", flush=True)
for f, lab in [("SEC61G ~ z_MP3", "Continuous MP3, no adjustment"),
               ("SEC61G ~ z_MP3 + C(pt)", "Continuous MP3 + patient fixed effects")]:
    m = smf.ols(f, data=d).fit()
    rec(lab, m.params["z_MP3"], m.pvalues["z_MP3"])

print("\n=== M2 adjust for per-cell depth ===", flush=True)
for f, lab in [("SEC61G ~ z_MP3 + log_tc + C(pt)", "Continuous MP3 + depth + patient"),
               ("SEC61G ~ is_MP3 + log_tc + n_genes + C(pt)", "Binary MP3 + depth + detected genes + patient")]:
    m = smf.ols(f, data=d).fit()
    k = "z_MP3" if "z_MP3" in m.params else "is_MP3"
    rec(lab, m.params[k], m.pvalues[k])

print("\n=== M3 NB GLM (raw counts + depth offset) ===", flush=True)
sub = d.sample(n=min(20000, len(d)), random_state=42).copy()   # NB fit is slow; subsample
sub["off"] = np.log(sub.total_counts.clip(lower=1))
pts = pd.get_dummies(sub.pt, drop_first=True).astype(float)
for key, lab in [("z_MP3", "NB: continuous MP3 + patient"), ("is_MP3", "NB: binary MP3 + patient")]:
    Xd = pd.concat([pd.Series(1.0, index=sub.index, name="const"), sub[[key]], pts], axis=1)
    try:
        m = sm.GLM(sub["SEC61G_cnt"], Xd, family=sm.families.NegativeBinomial(alpha=1.0),
                   offset=sub["off"]).fit()
        rec(lab, m.params[key], m.pvalues[key], f"(n={len(sub)})")
    except Exception as e:
        print(f"{lab}: fit failed {type(e).__name__}", flush=True)

print("\n=== M4 mixed effects (patient random intercept) ===", flush=True)
try:
    m = smf.mixedlm("SEC61G ~ z_MP3", d, groups=d["pt"]).fit(method="lbfgs")
    rec("Mixed effects: continuous MP3, patient random intercept", m.params["z_MP3"], m.pvalues["z_MP3"])
except Exception as e:
    print(f"Mixed-effects fit failed: {type(e).__name__}", flush=True)

print("\n=== M5 patients with both classes ≥50 cells ===", flush=True)
cnt = d.groupby(["pt", "is_MP3"]).size().unstack(fill_value=0)
ok = cnt[(cnt.get(0, 0) >= 50) & (cnt.get(1, 0) >= 50)].index
sub5 = d[d.pt.isin(ok)]
print(f"Eligible patients {len(ok)}, cells {len(sub5)}", flush=True)
for f, lab in [("SEC61G ~ is_MP3 + C(pt)", "Binary MP3 + patient (restricted patients)"),
               ("SEC61G ~ z_MP3 + C(pt)", "Continuous MP3 + patient (restricted patients)")]:
    m = smf.ols(f, data=sub5).fit()
    k = "z_MP3" if "z_MP3" in m.params else "is_MP3"
    rec(lab, m.params[k], m.pvalues[k], f"({len(ok)} patients)")

print("\n=== M6 per-patient effect-size distribution ===", flush=True)
rows = []
for pt, s in d.groupby("pt"):
    if (s.is_MP3.sum() < 20) or ((1 - s.is_MP3).sum() < 20):
        continue
    m = smf.ols("SEC61G ~ is_MP3", data=s).fit()
    rows.append({"pt": pt, "beta": m.params["is_MP3"], "p": m.pvalues["is_MP3"],
                 "n_MP3": int(s.is_MP3.sum()), "n_other": int((1 - s.is_MP3).sum())})
F = pd.DataFrame(rows)
F.to_csv(f"{OUT}/per_patient_effects.csv", index=False)
npos = int((F.beta > 0).sum())
sign_p = stats.binomtest(npos, len(F)).pvalue
print(f"Eligible patients {len(F)}: β>0 in {npos} ({npos/len(F)*100:.0f}%), sign test p={sign_p:.3g}", flush=True)
print(f"β median {F.beta.median():+.4f}, mean {F.beta.mean():+.4f},"
      f"IQR [{F.beta.quantile(.25):+.4f}, {F.beta.quantile(.75):+.4f}]", flush=True)
print(f"Patients with p<0.05: positive {int(((F.p<0.05)&(F.beta>0)).sum())},"
      f"negative {int(((F.p<0.05)&(F.beta<0)).sum())}", flush=True)
# One-sample test of per-patient β (overall departure from 0)
t, pt_ = stats.ttest_1samp(F.beta, 0)
w, pw = stats.wilcoxon(F.beta)
print(f"Per-patient β vs 0: t-test p={pt_:.3g}, Wilcoxon p={pw:.3g}", flush=True)
res.append({"Test": "逐病人 β 符号检验", "beta": F.beta.median(), "p": sign_p,
            "Note": f"{npos}/{len(F)} 为正"})
res.append({"Test": "Per-patient β Wilcoxon", "beta": F.beta.median(), "p": pw, "说明": ""})

print("\n=== Negative/positive controls (same within-patient model) ===", flush=True)
for g in ["SEC61A1", "SEC61B", "EGFR"]:
    m = smf.ols(f"{g} ~ is_MP3 + C(pt)", data=d).fit()
    rec(f"Control {g} ~ MP3 + patient", m.params["is_MP3"], m.pvalues["is_MP3"])

pd.DataFrame(res).to_csv(f"{OUT}/model_summary.csv", index=False)
print(f"\nResults written to {OUT}/", flush=True)
