"""C5 — Pre-check for ambient RNA contamination (assess risk without running decontX).

Motivation: OSM, IL1B and CXCL8 are among the most common myeloid ambient
transcripts, and the whole neutrophil axis rests on them. Ambient RNA has
to be ruled out before the OSM-priming states are re-derived.

Physical basis of ambient contamination: free mRNA from lysed cells is packaged into droplets at random, forming an approximately
uniform additive background within a sample, with each cell receiving an amount roughly proportional to its total counts. This can create false-positive detection,
but is unlikely to create within-sample differences between subtypes unless sequencing depth differs systematically between them.

This script checks three points to predict what decontX will show:

  Q1 Depth confounding: do OSM_priming / OSM_low have higher total counts and genes detected than other subtypes?
     If not, ambient cannot explain the subtype differences.
  Q2 After adjusting for depth, does the OSM subtype difference remain? (sample fixed effects + log total counts)
  Q3 Within the same sample, is myeloid OSM abundance (proxy for ambient source) correlated with neutrophil OSM
     detection rate? Strong correlation → background carry-over; weak → cell-intrinsic expression.

Plus Q4: compare OSM with neutrophil-intrinsic markers that cannot be ambient-derived,
     and with truly high-abundance myeloid genes (LYZ / S100A8), to see which class OSM behaviour matches.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
import statsmodels.formula.api as smf
from scipy import stats
import warnings, os
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/ambient_precheck"
os.makedirs(OUT, exist_ok=True)
KEY = ["OSM", "IL1B", "CXCL8"]
CTRL_NEU = ["FCGR3B", "CSF3R", "S100A8", "S100A9"]      # neutrophil-intrinsic markers
CTRL_HI = ["LYZ", "B2M", "ACTB"]                         # high-abundance genes prone to ambient contamination

print("Loading neutrophil object ...", flush=True)
neu = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad")
ann = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad")
assert (neu.obs_names == ann.obs_names).all()
neu.obs["sub"] = ann.obs["neu_subtype"].values
X = sp.csr_matrix(neu.X)
tc = np.asarray(X.sum(1)).ravel()
ng = np.asarray((X > 0).sum(1)).ravel()
gi = {g: i for i, g in enumerate(neu.var_names)}
d = pd.DataFrame({
    "sub": neu.obs["sub"].astype(str).values,
    "sample": neu.obs["sample_id"].astype(str).values,
    "dataset": neu.obs["dataset"].astype(str).values,
    "tc": tc, "ng": ng, "log_tc": np.log10(tc.clip(1)),
})
for g in KEY + CTRL_NEU + CTRL_HI:
    if g in gi:
        v = np.asarray(X[:, gi[g]].todense()).ravel()
        d[g] = v
        d[g + "_det"] = (v > 0).astype(int)
print(f"neutrophils {len(d)}, samples {d['sample'].nunique()}, subtypes {d['sub'].nunique()}\n", flush=True)

print("=== Q1 Sequencing depth by subtype (ambient uptake scales with this) ===", flush=True)
q1 = d.groupby("sub").agg(n=("tc", "size"), median_total_counts=("tc", "median"),
                          中位检出基因=("ng", "median")).sort_values("median_total_counts", ascending=False)
print(q1.to_string(), flush=True)
osm_subs = ["Neu_OSM_priming", "Neu_OSM_low"]
a_ = d[d["sub"].isin(osm_subs)].tc; b_ = d[~d["sub"].isin(osm_subs)].tc
print(f"\n  OSM-related subtypes median total counts {a_.median():.0f} vs others {b_.median():.0f}"
      f"MW p={stats.mannwhitneyu(a_, b_).pvalue:.3g}", flush=True)
print("→ if OSM subtypes are not deeper, ambient cannot explain their OSM enrichment", flush=True)
q1.to_csv(f"{OUT}/Q1_depth_by_subtype.csv")

print("\n=== Q2 Gene enrichment in OSM subtypes after adjusting for sample and depth ===", flush=True)
d["is_osm_sub"] = d["sub"].isin(osm_subs).astype(int)
rows = []
for g in KEY + CTRL_NEU + CTRL_HI:
    if g not in d.columns: continue
    m0 = smf.ols(f"{g} ~ is_osm_sub", data=d).fit()
    m1 = smf.ols(f"{g} ~ is_osm_sub + log_tc + C(sample)", data=d).fit()
    rows.append(dict(gene=g, beta_raw=m0.params["is_osm_sub"], p_raw=m0.pvalues["is_osm_sub"],
                     beta_adj=m1.params["is_osm_sub"], p_adj=m1.pvalues["is_osm_sub"],
                     det_rate=d[g + "_det"].mean()))
    print(f"{g:<8} detection rate {d[g+'_det'].mean()*100:5.1f}%"
          f"unadjusted β={m0.params['is_osm_sub']:+.4f} (p={m0.pvalues['is_osm_sub']:.2g})"
          f"adjusted for sample+depth β={m1.params['is_osm_sub']:+.4f} (p={m1.pvalues['is_osm_sub']:.2g})", flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/Q2_subtype_enrichment_adjusted.csv", index=False)

print("\n=== Q3 Same-sample myeloid OSM abundance (ambient source proxy) vs neutrophil OSM detection rate ===", flush=True)
mye = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad", backed="r")
mv = {g: i for i, g in enumerate(mye.var_names)}
if "OSM" in mv:
    sub_ids = sorted(set(d["sample"]) & set(mye.obs["sample_id"].astype(str)))
    print(f"samples shared by both sides: {len(sub_ids)}", flush=True)
    mye_mem = mye.to_memory()
    MX = sp.csr_matrix(mye_mem.X)
    msamp = mye_mem.obs["sample_id"].astype(str).values
    rec = []
    for s in sub_ids:
        mm = msamp == s
        if mm.sum() < 20: continue
        osm_mye = float(np.asarray(MX[mm, mv["OSM"]].todense()).mean())
        nn = d["sample"] == s
        if nn.sum() < 20: continue
        rec.append(dict(sample=s, n_neu=int(nn.sum()), n_mye=int(mm.sum()),
                        mye_OSM_mean=osm_mye,
                        neu_OSM_detrate=float(d.loc[nn, "OSM_det"].mean()),
                        neu_OSM_mean=float(d.loc[nn, "OSM"].mean()),
                        neu_S100A8_detrate=float(d.loc[nn, "S100A8_det"].mean()) if "S100A8_det" in d else np.nan))
    R = pd.DataFrame(rec)
    R.to_csv(f"{OUT}/Q3_myeloid_vs_neutrophil_OSM.csv", index=False)
    if len(R) > 5:
        r1, p1 = stats.spearmanr(R.mye_OSM_mean, R.neu_OSM_detrate)
        r2, p2 = stats.spearmanr(R.mye_OSM_mean, R.neu_OSM_mean)
        print(f"n samples {len(R)}", flush=True)
        print(f"myeloid mean OSM vs neutrophil OSM detection rate: Spearman ρ={r1:+.3f} (p={p1:.3g})", flush=True)
        print(f"myeloid mean OSM vs neutrophil mean OSM expression: Spearman ρ={r2:+.3f} (p={p2:.3g})", flush=True)
        print("→ high ρ (>0.6) suggests background carry-over; low ρ suggests neutrophil-intrinsic expression", flush=True)
else:
    print("OSM absent from myeloid object; skip", flush=True)

print(f"\nResults written to {OUT}/", flush=True)
