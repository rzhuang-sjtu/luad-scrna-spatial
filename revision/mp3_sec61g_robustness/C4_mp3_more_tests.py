"""C4 — Second battery of SEC61G × MP3 tests: five methods independent of C3.

T1 Genome-wide benchmark: run the same within-patient model on all genes; report SEC61G β percentile.
   Uses within-transformation (demean within patient), algebraically equivalent to OLS + patient fixed effects,
   but vectorized over all genes in seconds. Decides specificity:
   if SEC61G is mid-ranked, +0.034 is a genome-wide pattern, not a property of SEC61G.

T2 Permutation test: shuffle MP3 scores within each patient, recompute β, build the null distribution.
   No normality/homoscedasticity/linearity assumptions; exact p-value.

T3 Pseudobulk pairing: within each patient, bin cells by MP3 quartiles, take mean expression per quartile,
   then paired Wilcoxon Q4 vs Q1 across patients. Standard robust single-cell DE approach;
   removes pseudo-replication (one observation per patient, not per cell).

T4 Within-patient Spearman + combine: fully non-parametric; per-patient rank correlation then combine (sign test + Fisher).

T5 Zero-inflated hurdle: split into detection (logistic) vs expression given detection (OLS);
   ask whether MP3 affects detection rate or expression intensity.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
import statsmodels.formula.api as smf
from scipy import stats
import warnings, os
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/mp3_sec61g_robustness"
os.makedirs(OUT, exist_ok=True)
TARGET = "SEC61G"
PANEL = ["SEC61G", "EGFR", "VOPP1", "LANCL2", "SEC61A1", "SEC61B",
         "KRT7", "S100A10", "ANXA2", "VIM", "B2M", "GAPDH", "NAPSA"]

print("Loading ...", flush=True)
a = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_malignant_scored.h5ad")
pt = (a.obs["dataset"].astype(str) + "__" + a.obs["patient_id"].astype(str)).values
mp3 = a.obs["MP3_score"].values.astype(float)
X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
genes = np.array(a.var_names)
gi = {g: i for i, g in enumerate(genes)}
n, G = X.shape
print(f"Cells {n}, genes {G}, patients {len(np.unique(pt))}", flush=True)

codes, uniq = pd.factorize(pt)
K = len(uniq)
cnt = np.bincount(codes, minlength=K).astype(float)
def demean(v):                       # v: (n,) or (n,G)
    if v.ndim == 1:
        s = np.bincount(codes, weights=v, minlength=K)
        return v - (s / cnt)[codes]
    s = np.zeros((K, v.shape[1]))
    np.add.at(s, codes, v)
    return v - (s / cnt[:, None])[codes]

x_dm = demean(mp3)
sx = x_dm.std()
x_dm_z = x_dm / sx                   # β unit = per 1 SD of within-patient MP3 variation
X_dm = demean(X)
denom = (x_dm_z ** 2).sum()
beta_all = (X_dm * x_dm_z[:, None]).sum(axis=0) / denom

print("\n=== T1 genome-wide benchmark (within-patient β, all genes) ===", flush=True)
b_t = beta_all[gi[TARGET]]
pct = (beta_all < b_t).mean() * 100
print(f"{TARGET} β = {b_t:+.4f}, percentile {pct:.1f} among {G} genes", flush=True)
print(f"Genome-wide β: median {np.median(beta_all):+.4f}, mean {beta_all.mean():+.4f},"
      f"IQR [{np.percentile(beta_all,25):+.4f}, {np.percentile(beta_all,75):+.4f}]", flush=True)
print(f"Genes with larger β: {int((beta_all > b_t).sum())} ({100-pct:.1f}%)", flush=True)
print("\n  Panel-gene percentiles:", flush=True)
rows = []
for g in PANEL:
    if g not in gi: continue
    b = beta_all[gi[g]]; p_ = (beta_all < b).mean() * 100
    print(f"{g:<9} β={b:+.4f}  percentile {p_:5.1f}", flush=True)
    rows.append(dict(gene=g, beta=b, percentile=p_))
pd.DataFrame({"gene": genes, "beta_within_patient": beta_all}).to_csv(
    f"{OUT}/T1_genomewide_beta.csv.gz", index=False)

print("\n=== T2 permutation test (shuffle MP3 within patient, 2000×) ===", flush=True)
rng = np.random.default_rng(20260804)
y_dm = X_dm[:, gi[TARGET]]
order = np.argsort(codes, kind="stable")
starts = np.r_[0, np.cumsum(np.bincount(codes, minlength=K))]
null = np.empty(2000)
for it in range(2000):
    xp = np.empty(n)
    for k in range(K):
        idx = order[starts[k]:starts[k + 1]]
        xp[idx] = rng.permutation(mp3[idx])
    xp_dm = demean(xp) / sx
    null[it] = (y_dm * xp_dm).sum() / (xp_dm ** 2).sum()
p_perm = ((np.abs(null) >= abs(b_t)).sum() + 1) / (len(null) + 1)
print(f"Observed β={b_t:+.4f}; null mean {null.mean():+.5f}, SD {null.std():.5f}", flush=True)
print(f"Two-sided permutation p = {p_perm:.4g}", flush=True)

print("\n=== T3 pseudobulk pairing (per-patient MP3 quartiles, Q4 vs Q1) ===", flush=True)
df = pd.DataFrame({"pt": pt, "mp3": mp3})
res3 = []
for g in PANEL:
    if g not in gi: continue
    df["y"] = X[:, gi[g]]
    q4, q1, keep = [], [], 0
    for p_, s in df.groupby("pt"):
        if len(s) < 40: continue
        lo, hi = s.mp3.quantile(0.25), s.mp3.quantile(0.75)
        a1 = s.loc[s.mp3 <= lo, "y"]; a4 = s.loc[s.mp3 >= hi, "y"]
        if len(a1) < 10 or len(a4) < 10: continue
        q1.append(a1.mean()); q4.append(a4.mean()); keep += 1
    q1, q4 = np.array(q1), np.array(q4)
    w, pw = stats.wilcoxon(q4, q1)
    d_ = (q4 - q1)
    res3.append(dict(gene=g, n_pt=keep, delta_median=np.median(d_),
                     frac_pos=(d_ > 0).mean(), p_wilcoxon=pw))
    print(f"{g:<9} {keep} patients  Δ(Q4−Q1) median {np.median(d_):+.4f}"
          f"{(d_>0).mean()*100:3.0f}% positive  Wilcoxon p={pw:.3g}", flush=True)
pd.DataFrame(res3).to_csv(f"{OUT}/T3_pseudobulk_quartile.csv", index=False)

print("\n=== T4 within-patient Spearman (non-parametric) + combine ===", flush=True)
res4 = []
for g in PANEL:
    if g not in gi: continue
    y = X[:, gi[g]]
    rhos = []
    for k in range(K):
        idx = order[starts[k]:starts[k + 1]]
        if len(idx) < 30: continue
        r = stats.spearmanr(mp3[idx], y[idx])[0]
        if np.isfinite(r): rhos.append(r)
    rhos = np.array(rhos)
    npos = int((rhos > 0).sum())
    sp_ = stats.binomtest(npos, len(rhos)).pvalue
    wp = stats.wilcoxon(rhos)[1]
    res4.append(dict(gene=g, n_pt=len(rhos), median_rho=np.median(rhos),
                     frac_pos=npos / len(rhos), p_sign=sp_, p_wilcoxon=wp))
    print(f"{g:<9} {len(rhos)} patients  median ρ={np.median(rhos):+.4f}"
          f"{npos/len(rhos)*100:3.0f}% positive  sign p={sp_:.3g}  Wilcoxon p={wp:.3g}", flush=True)
pd.DataFrame(res4).to_csv(f"{OUT}/T4_within_patient_spearman.csv", index=False)

print("\n=== T5 zero-inflated split (detection rate vs expression given detection) ===", flush=True)
dd = pd.DataFrame({"pt": pt, "z": (mp3 - mp3.mean()) / mp3.std()})
for g in ["SEC61G", "EGFR", "B2M", "KRT7"]:
    if g not in gi: continue
    y = X[:, gi[g]]
    dd["det"] = (y > 0).astype(int); dd["val"] = y
    m1 = smf.logit("det ~ z + C(pt)", data=dd).fit(disp=0)
    sub = dd[dd.det == 1]
    m2 = smf.ols("val ~ z + C(pt)", data=sub).fit()
    print(f"{g:<9} detection rate {dd.det.mean()*100:4.1f}%"
          f"logistic β={m1.params['z']:+.4f} (p={m1.pvalues['z']:.2g})   "
          f"OLS given detection β={m2.params['z']:+.4f} (p={m2.pvalues['z']:.2g})", flush=True)

print(f"\nResults written to {OUT}/", flush=True)
