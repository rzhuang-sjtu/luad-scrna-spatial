"""D3 — Robustness of OSM subtype differences (multiple methods; ambient-corrected data).

An initial test with raw counts + OLS + sample fixed effects + log depth gave OSM p=0.74,
suggesting no subtype difference. That test has four weaknesses addressed here:

  Weakness 1 power: median 33 neutrophils per sample across 8 subtypes leaves almost no within-sample contrast,
        yet 131 sample fixed effects exhaust degrees of freedom.
        Approach: quantify per-sample comparability, then restrict analysis to samples with enough cells in both groups.
  Weakness 2 distribution: raw counts in OLS are inappropriate → negative binomial GLM + log-depth offset.
  Weakness 3 over-adjustment: OSM subtypes have systematically lower depth; depth may be a mediator rather than a confounder
        → report both depth-unadjusted and depth-adjusted versions.
  Weakness 4 no specificity baseline → genome-wide β distribution for this contrast; report OSM percentile rank.

Also: pseudobulk pairing (per-sample group means, paired Wilcoxon across samples); fully non-parametric, no pseudo-replication.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, glob, os
import anndata as ad
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/decontx"
OSM_SUBS = ["Neu_OSM_priming", "Neu_OSM_low"]
PANEL = ["OSM", "IL1B", "CXCL8", "CSF3R", "FCGR3B", "S100A8", "S100A9", "VEGFA", "MMP9"]

print("Loading decontaminated matrix ...", flush=True)
mats, bcs = [], []
for d in sorted(glob.glob(f"{OUT}/output/*/")):
    f, b = os.path.join(d, "decontaminated.mtx"), os.path.join(d, "neu_barcodes.tsv")
    if os.path.exists(f) and os.path.exists(b):
        mats.append(sp.csr_matrix(sio.mmread(f)).T)
        bcs += [l.strip() for l in open(b)]
D = sp.vstack(mats).tocsr()

raw = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad", backed="r")
ann = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad", backed="r")
mg = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_merged_annotated.h5ad", backed="r")
genes = np.array(mg.var_names); gi = {g: i for i, g in enumerate(genes)}
idx = pd.Index(raw.obs_names).get_indexer(bcs); ok = idx >= 0
D = D[ok]; idx = idx[ok]
sub = ann.obs["neu_subtype"].astype(str).values[idx]
samp = raw.obs["sample_id"].astype(str).values[idx]
tc = np.asarray(D.sum(1)).ravel()
is_osm = np.isin(sub, OSM_SUBS).astype(int)
print(f"cells {D.shape[0]}, samples {len(np.unique(samp))}, OSM subtype {is_osm.sum()}\n", flush=True)

print("=== Per-sample comparability (within-sample contrast needs cells in both groups) ===", flush=True)
t = pd.DataFrame({"samp": samp, "is_osm": is_osm})
ct = t.groupby(["samp", "is_osm"]).size().unstack(fill_value=0)
ct.columns = ["n_other", "n_osm"] if 0 in ct.columns else ct.columns
print(f"total samples {len(ct)}", flush=True)
for k in [1, 3, 5, 10, 20]:
    n = int(((ct.iloc[:, 0] >= k) & (ct.iloc[:, 1] >= k)).sum())
    cells = int(ct[(ct.iloc[:, 0] >= k) & (ct.iloc[:, 1] >= k)].sum().sum())
    print(f"samples with >={k:<2} cells in each group: {n:>3} (covering {cells} cells)", flush=True)
good = ct[(ct.iloc[:, 0] >= 5) & (ct.iloc[:, 1] >= 5)].index
print(f"→ subsequent restricted-sample analyses use the {len(good)} samples with >=5 cells per group", flush=True)

d = pd.DataFrame({"samp": samp, "is_osm": is_osm, "tc": tc,
                  "log_tc": np.log10(np.clip(tc, 1, None)), "off": np.log(np.clip(tc, 1, None))})
sel = d.samp.isin(good).values

def nb_beta(y, X, off):
    try:
        m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=1.0), offset=off).fit()
        return m.params["is_osm"], m.pvalues["is_osm"]
    except Exception:
        return np.nan, np.nan

print("\n=== Weaknesses 2+3: negative binomial GLM (correct count distribution), four settings ===", flush=True)
print(f"{'gene':<8}{'all samples no depth':>18}{'all samples w/ depth':>18}{'restricted no depth':>20}{'restricted w/ depth':>20}", flush=True)
rows = []
dum_all = pd.get_dummies(d.samp, drop_first=True).astype(float)
dum_sel = pd.get_dummies(d.loc[sel, "samp"], drop_first=True).astype(float)
for g in PANEL:
    if g not in gi: continue
    y = np.asarray(D[:, gi[g]].todense()).ravel()
    out = []
    for use_sel in [False, True]:
        m_ = sel if use_sel else np.ones(len(d), bool)
        dum = dum_sel if use_sel else dum_all
        base = pd.concat([pd.Series(1.0, index=dum.index, name="const"),
                          pd.Series(d.loc[m_, "is_osm"].values, index=dum.index, name="is_osm"),
                          dum], axis=1)
        for use_off in [False, True]:
            off = d.loc[m_, "off"].values if use_off else np.zeros(m_.sum())
            out.append(nb_beta(y[m_], base, off))
    print(f"{g:<8}" + "".join(f"{b:+9.3f}(p={p:.1g})" if np.isfinite(b) else f"{'NA':>18}"
                              for b, p in [out[0], out[1], out[2], out[3]]), flush=True)
    rows.append(dict(gene=g, b_all_nooff=out[0][0], p_all_nooff=out[0][1],
                     b_all_off=out[1][0], p_all_off=out[1][1],
                     b_sel_nooff=out[2][0], p_sel_nooff=out[2][1],
                     b_sel_off=out[3][0], p_sel_off=out[3][1]))
pd.DataFrame(rows).to_csv(f"{OUT}/D3_nb_models.csv", index=False)

print("\n=== Pseudobulk pairing (per-sample group means, Wilcoxon across samples) ===", flush=True)
cpm = D.multiply(1e4 / np.clip(tc, 1, None)[:, None]).tocsr()
res = []
for g in PANEL:
    if g not in gi: continue
    y = np.log1p(np.asarray(cpm[:, gi[g]].todense()).ravel())
    a_, b_ = [], []
    for s in good:
        m_ = samp == s
        a_.append(y[m_ & (is_osm == 1)].mean()); b_.append(y[m_ & (is_osm == 0)].mean())
    a_, b_ = np.array(a_), np.array(b_)
    w, p = stats.wilcoxon(a_, b_)
    print(f"{g:<8} {len(a_)} samples  Δ median {np.median(a_-b_):+.4f}"
          f"{(a_>b_).mean()*100:3.0f}% positive  Wilcoxon p={p:.3g}", flush=True)
    res.append(dict(gene=g, n=len(a_), delta=np.median(a_ - b_),
                    frac_pos=(a_ > b_).mean(), p=p))
pd.DataFrame(res).to_csv(f"{OUT}/D3_pseudobulk.csv", index=False)

print("\n=== Weakness 4: genome-wide specificity baseline (within-group demeaning, restricted samples) ===", flush=True)
Ds = D[sel]; ss = samp[sel]; xs = is_osm[sel].astype(float)
codes, uq = pd.factorize(ss); K = len(uq)
cntk = np.bincount(codes, minlength=K).astype(float)
def dm(v):
    if v.ndim == 1:
        return v - (np.bincount(codes, weights=v, minlength=K) / cntk)[codes]
    s = np.zeros((K, v.shape[1])); np.add.at(s, codes, v)
    return v - (s / cntk[:, None])[codes]
Y = np.log1p(Ds.multiply(1e4 / np.clip(np.asarray(Ds.sum(1)).ravel(), 1, None)[:, None]).toarray())
xd = dm(xs); xd /= xd.std()
Yd = dm(Y)
beta = (Yd * xd[:, None]).sum(0) / (xd ** 2).sum()
for g in PANEL:
    if g not in gi: continue
    b = beta[gi[g]]; pc = (beta < b).mean() * 100
    print(f"{g:<8} β={b:+.4f}  genome-wide percentile {pc:5.1f}", flush=True)
pd.DataFrame({"gene": genes, "beta": beta}).to_csv(f"{OUT}/D3_genomewide_beta.csv.gz", index=False)
print(f"\nResults written to {OUT}/", flush=True)
