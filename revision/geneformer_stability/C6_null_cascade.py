"""C6 — Expression-matched random genes through a five-endpoint cascade (rebuild; write to a durable directory).

Motivation: the cascade is only informative if a random gene set of the
same size and expression range does not pass it. Several sets of ~27
expression-matched random genes are put through the identical cascade.

Improvements over the first version:
  1. Write results to results/null_cascade/ (previously a temp dir that was cleaned up)
  2. Besides the binary «passed all five» criterion, add a **continuous** metric — number of endpoints passed (0–5).
     Binary has very low power with only 27 candidates and 1 full pass (Fisher p=0.113),
     discarding much information; the continuous metric is far more powerful.
  3. Add per-endpoint statistic-level comparisons (not only threshold pass/fail) and a permutation distribution for random 27-gene sets.

Five endpoints use the manuscript criteria:
  E1 DepMap 24Q2: LUAD mean Chronos < -0.5 and one-sided MW p<0.05
  E2 TCGA tumour vs normal: log2FC>0 and Wilcoxon p<0.05
  E3 TCGA survival: univariate Cox HR>1 and p<0.05
  E4 Spatial tumour-intrinsic ROI: Δ>0 and MW p<0.05 (discovery cohort E-MTAB-13530)
  E5 Immunotherapy: higher in non-responders; ≥1 of 3 cohorts p<0.05
"""
import numpy as np, pandas as pd, os, time, warnings
import anndata as ad, scipy.sparse as sp
from scipy import stats
from lifelines import CoxPHFitter
warnings.filterwarnings("ignore")

OUT = "${PROJECT_ROOT}/results/null_cascade"
os.makedirs(OUT, exist_ok=True)
GF = "${PROJECT_ROOT}/results/fig8_geneformer"
N_CTRL = 100
t0 = time.time()
def log(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

cand = pd.read_csv("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn_500_diff/8A_500_candidate_pool.csv")
CAND = sorted(cand.Gene_name.unique()); log(f"Candidates {len(CAND)}")

tpm = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv", index_col=0)
tpm = tpm[~tpm.index.duplicated()].astype(np.float32)
tum = [c for c in tpm.columns if c[13:15] == "01"]
nor = [c for c in tpm.columns if c[13:15] == "11"]
expr_mean = np.log2(tpm[tum] + 1).mean(axis=1)
log(f"TCGA {tpm.shape} tumour {len(tum)} normal {len(nor)}")

D = "${DATA_ROOT}/depmap/24Q2/"
mod = pd.read_csv(D + "Model.csv", low_memory=False)
luad = set(mod.loc[mod.OncotreeCode == "LUAD", "ModelID"])
cr = pd.read_csv(D + "CRISPRGeneEffect.csv", index_col=0)
cr.columns = [c.split(" (")[0] for c in cr.columns]
cr = cr.loc[:, ~cr.columns.duplicated()]
is_luad = cr.index.isin(luad); log(f"DepMap {cr.shape} LUAD {is_luad.sum()}")

a_roi = ad.read_h5ad("${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad")
a_c2l = ad.read_h5ad("${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad")
common = a_roi.obs_names.intersection(a_c2l.obs_names)
a_roi = a_roi[common]; a_c2l = a_c2l[common]
z = lambda v: (np.asarray(v, float) - np.asarray(v, float).mean()) / (np.asarray(v, float).std(ddof=0) + 1e-12)
roi = (z(a_roi.obs["MP3_score"].values) > 0.5) & \
      (z(a_c2l.obsm["q05_cell_abundance_w_sf"]["q05cell_abundance_w_sf_Malignant"].values) > 0.5)
SP_G = set(a_roi.var_names); log(f"Spatial {len(common)} spots, ROI {roi.sum()}")

PD_ = "${WORK_ROOT}/luad_figures/fig_treatment/"
c1 = pd.read_csv("${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
c1 = c1[~c1.index.duplicated()]
s1 = pd.read_csv(PD_ + "gse207422_mp_scores.csv").set_index("Sample")
d2 = pd.read_csv("${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz", sep="\t", index_col=0).groupby(level=0).max()
c2 = np.log2(d2.div(d2.sum(axis=0), axis=1) * 1e6 + 1)
s2 = pd.read_csv(PD_ + "gse126044_mp_scores.csv").set_index("Sample")
d3 = pd.read_csv("${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", sep="\t", index_col=0)
d3.index = d3.index.astype(str).str.split(".").str[0]; d3 = d3.groupby(level=0).max()
c3 = np.log2(d3 + 1)
s3 = pd.read_csv(PD_ + "gse135222_mp_scores.csv").set_index("Sample")
mm = pd.read_csv("${PROJECT_ROOT}/data/external/gse135222_ensg2symbol.csv")
mm["e"] = mm["ensg"].astype(str).str.split(".").str[0]
sym2ensg = dict(zip(mm["symbol"], mm["e"]))
IO = [("GSE207422", c1, s1, "MPR", "NMPR", None), ("GSE126044", c2, s2, "R", "NR", None),
      ("GSE135222", c3, s3, "R", "NR", sym2ensg)]
log("IO cohorts ready")

clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin = clin[clin.sample_type == "Primary Tumor"].drop_duplicates("sample_barcode").set_index("sample_barcode")
sur = clin.loc[[s for s in tum if s in clin.index]].copy()
sur["event"] = (sur.vital_status.str.strip().str.lower() == "dead").astype(int)
sur["time"] = np.where(sur.event == 1, sur.days_to_death, sur.days_to_last_follow_up)
sur = sur[sur.time.notna() & (sur.time > 0)]
surv_s = list(sur.index); log(f"Survival n={len(sur)} events {int(sur.event.sum())}")

universe = sorted(set(tpm.index) & set(cr.columns) & SP_G)
cand_in = [g for g in CAND if g in universe]
pool = pd.Series(expr_mean[universe]).sort_values(); order = list(pool.index)
pos = {g: i for i, g in enumerate(order)}
ctrl = set()
for g in cand_in:
    i = pos[g]; lo, hi = max(0, i - 400), min(len(order), i + 400)
    nb = [x for x in order[lo:hi] if x not in CAND and x not in ctrl]
    nb.sort(key=lambda x: abs(expr_mean[x] - expr_mean[g]))
    ctrl.update(nb[:N_CTRL])
TEST = cand_in + sorted(ctrl)
log(f"Universe {len(universe)}, candidates {len(cand_in)}, controls {len(ctrl)}, total {len(TEST)}")

T = np.log2(tpm[tum] + 1); N = np.log2(tpm[nor] + 1)
sp_genes = [g for g in TEST if g in SP_G]
SPX = a_roi[:, sp_genes].X
SPX = SPX.toarray() if sp.issparse(SPX) else np.asarray(SPX)
sp_col = {g: i for i, g in enumerate(sp_genes)}
surv_expr = np.log2(tpm.loc[TEST, surv_s] + 1)

rows = []
for k, g in enumerate(TEST):
    if k % 400 == 0: log(f"  {k}/{len(TEST)}")
    r = {"gene": g, "is_candidate": g in CAND, "expr_mean": float(expr_mean[g])}
    v = cr[g].values; la, lb = v[is_luad], v[~is_luad]
    la, lb = la[~np.isnan(la)], lb[~np.isnan(lb)]
    r["dep_luad_mean"] = float(la.mean()); r["dep_delta"] = float(la.mean() - lb.mean())
    r["dep_p"] = float(stats.mannwhitneyu(la, lb, alternative="less").pvalue)
    r["E1"] = bool(la.mean() < -0.5 and r["dep_p"] < 0.05)
    x, y = T.loc[g].values, N.loc[g].values
    r["tn_log2fc"] = float(x.mean() - y.mean())
    r["tn_p"] = float(stats.mannwhitneyu(x, y).pvalue)
    r["E2"] = bool(r["tn_log2fc"] > 0 and r["tn_p"] < 0.05)
    e = surv_expr.loc[g].values
    dd = pd.DataFrame({"time": sur.time.values, "event": sur.event.values,
                       "x": (e - e.mean()) / (e.std() + 1e-12)})
    try:
        c = CoxPHFitter().fit(dd, duration_col="time", event_col="event")
        r["os_hr"] = float(c.summary.loc["x", "exp(coef)"]); r["os_p"] = float(c.summary.loc["x", "p"])
    except Exception:
        r["os_hr"], r["os_p"] = np.nan, np.nan
    r["E3"] = bool(r["os_hr"] > 1 and r["os_p"] < 0.05)
    if g in sp_col:
        xx = SPX[:, sp_col[g]]; aa, bb = xx[roi], xx[~roi]
        r["sp_delta"] = float(aa.mean() - bb.mean())
        r["sp_p"] = float(stats.mannwhitneyu(aa, bb).pvalue)
        r["E4"] = bool(r["sp_delta"] > 0 and r["sp_p"] < 0.05)
    else:
        r["sp_delta"], r["sp_p"], r["E4"] = np.nan, np.nan, False
    nsig = ntest = 0; io_min_p = 1.0
    for name, mat, sc, pos_l, neg_l, mp_ in IO:
        key = mp_.get(g) if mp_ else g
        if key is None or key not in mat.index: continue
        smp = [s for s in mat.columns if s in sc.index]
        if len(smp) < 6: continue
        resp = sc.loc[smp, "response_group"]
        vals = mat.loc[key, smp].astype(float)
        gp, gn = vals[resp.values == pos_l], vals[resp.values == neg_l]
        if len(gp) < 3 or len(gn) < 3: continue
        ntest += 1
        p_ = stats.mannwhitneyu(gn, gp, alternative="greater").pvalue
        io_min_p = min(io_min_p, p_)
        if p_ < 0.05: nsig += 1
    r["io_ntested"] = ntest; r["io_nsig"] = nsig; r["io_min_p"] = io_min_p
    r["E5"] = bool(nsig >= 1)
    r["n_pass"] = sum(int(r[f"E{i}"]) for i in range(1, 6))
    r["pass_all5"] = r["n_pass"] == 5
    rows.append(r)

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/null_cascade_results.csv", index=False)
C, K = res[res.is_candidate], res[~res.is_candidate]
log(f"Done, candidates {len(C)} controls {len(K)}")

print("\n" + "=" * 72)
print("Expression-matched null: binary vs continuous criteria")
print("=" * 72)
a, b = int(C.pass_all5.sum()), len(C) - int(C.pass_all5.sum())
c_, d_ = int(K.pass_all5.sum()), len(K) - int(K.pass_all5.sum())
print(f"\n[Binary] Passed all five: candidates {a}/{len(C)} ({a/len(C)*100:.1f}%)"
      f"controls {c_}/{len(K)} ({c_/len(K)*100:.2f}%)"
      f"Fisher p={stats.fisher_exact([[a,b],[c_,d_]],alternative='greater')[1]:.4f}")
print(f"\n[Continuous] Endpoints passed (0–5):")
print(f"Candidates  mean {C.n_pass.mean():.3f} ± {C.n_pass.std():.3f}  median {C.n_pass.median():.0f}")
print(f"Controls  mean {K.n_pass.mean():.3f} ± {K.n_pass.std():.3f}  median {K.n_pass.median():.0f}")
u, p = stats.mannwhitneyu(C.n_pass, K.n_pass, alternative="greater")
dcoh = (C.n_pass.mean() - K.n_pass.mean()) / np.sqrt((C.n_pass.var() + K.n_pass.var()) / 2)
print(f"Mann-Whitney one-sided p={p:.3g}   Cohen's d={dcoh:.3f}")
print(f"\n  {'n_pass':<8}{'cand':>9}{'ctrl':>9}")
for k in range(6):
    print(f"  {k:<8}{(C.n_pass==k).mean()*100:>8.1f}%{(K.n_pass==k).mean()*100:>8.1f}%")
rng = np.random.default_rng(7); v = K.n_pass.values
sim = np.array([v[rng.choice(len(v), len(C), replace=False)].mean() for _ in range(20000)])
print(f"\n  Mean endpoints passed over 20,000 random {len(C)}-gene sets: {sim.mean():.3f},"
      f"95% interval [{np.percentile(sim,2.5):.3f}, {np.percentile(sim,97.5):.3f}]")
print(f"Observed candidate mean {C.n_pass.mean():.3f}, empirical p={(sim>=C.n_pass.mean()).mean():.4g}")
print(f"\n[Per endpoint] pass rate candidates vs controls")
for e, lab in [("E1","DepMap essentiality"),("E2","TCGA tumour high expression"),("E3","TCGA 生存"),
               ("E4","Spatial tumour-intrinsic ROI"),("E5","免疫治疗非响应")]:
    aa,bb = int(C[e].sum()), len(C)-int(C[e].sum())
    cc,dd2 = int(K[e].sum()), len(K)-int(K[e].sum())
    print(f"  {lab:<18} {C[e].mean()*100:5.1f}%  {K[e].mean()*100:5.1f}%  "
          f"Fisher p={stats.fisher_exact([[aa,bb],[cc,dd2]],alternative='greater')[1]:.4f}")
print(f"\nResults written to {OUT}/null_cascade_results.csv")
