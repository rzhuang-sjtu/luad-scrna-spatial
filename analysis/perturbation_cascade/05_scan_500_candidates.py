"""Lightweight DepMap + TCGA scan of 27 new 500-cell candidates (no plots).

For each gene:
  - DepMap CRISPR: LUAD mean effect, non-LUAD mean, MW one-sided p (LUAD<other)
  - TCGA T vs N: log2FC, Wilcoxon p
  - TCGA OS: KM logrank (median split) + Cox HR/p (continuous log2 expr)

Output: combined ranked table.
"""
import re
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
warnings.filterwarnings("ignore")

DEPMAP = Path("${DATA_ROOT}/depmap/24Q2")
TCGA   = Path("${DATA_ROOT}/TCGA_LUAD_analysis")
CAND   = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn_500_diff/8A_500_candidate_pool.csv")
OUT    = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn_500_diff")

PRIORITY = {"GABARAP", "CXCL8", "ANGPTL4", "BSG", "HK2", "SLC2A1", "SLC25A5", "SRSF9"}

cands = pd.read_csv(CAND)["Gene_name"].tolist()
print(f"loaded {len(cands)} candidates from 500-cell pool")

# ----- DepMap -----
print("\n[A] DepMap 24Q2")
model = pd.read_csv(DEPMAP / "Model.csv", low_memory=False)
luad_models = set(model.loc[model["OncotreeCode"] == "LUAD", "ModelID"])
print(f"  LUAD lines: {len(luad_models)}")
crispr_header = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", nrows=0).columns.tolist()
def name_of(c):
    m = re.match(r"^([A-Za-z0-9\-]+)\s*\(\d+\)$", c)
    return m.group(1) if m else c
sym2col = {name_of(c): c for c in crispr_header[1:]}
sel = [crispr_header[0]] + [sym2col[g] for g in cands if g in sym2col]
crispr = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", usecols=sel)
crispr = crispr.rename(columns={crispr_header[0]: "ModelID"})
crispr = crispr.rename(columns={sym2col[g]: g for g in cands if g in sym2col})
crispr["is_LUAD"] = crispr["ModelID"].isin(luad_models)

dep_rows = []
for g in cands:
    if g not in crispr.columns:
        dep_rows.append({"gene": g, "depmap_in_data": False}); continue
    a = crispr.loc[crispr["is_LUAD"], g].dropna()
    b = crispr.loc[~crispr["is_LUAD"], g].dropna()
    p = stats.mannwhitneyu(a, b, alternative="less").pvalue if len(a) >= 5 and len(b) >= 5 else np.nan
    dep_rows.append({"gene": g, "depmap_in_data": True,
                     "luad_mean": a.mean(), "luad_n": len(a),
                     "other_mean": b.mean(), "other_n": len(b),
                     "delta_LUAD_minus_other": a.mean() - b.mean(),
                     "mw_p_LUAD_lt_other": p})
dep_df = pd.DataFrame(dep_rows)

# ----- TCGA -----
print("\n[B] TCGA-LUAD T vs N + OS")
cln = pd.read_csv(TCGA / "TCGA_LUAD_clinical.csv").rename(columns={"sample_barcode": "sample"})
tpm_rows = []
for chunk in pd.read_csv(TCGA / "TCGA_LUAD_TPM_matrix.csv", chunksize=5000):
    chunk = chunk.rename(columns={chunk.columns[0]: "gene"})
    sub = chunk[chunk["gene"].isin(cands)]
    if len(sub): tpm_rows.append(sub)
tpm = pd.concat(tpm_rows, ignore_index=True).set_index("gene")
print(f"  TPM matrix shape: {tpm.shape}; genes found: {tpm.index.tolist()}")

st_map = dict(zip(cln["sample"], cln["sample_type"]))
tcga_rows = []
for g in cands:
    row = {"gene": g, "tcga_in_data": g in tpm.index}
    if g not in tpm.index:
        tcga_rows.append(row); continue
    expr = tpm.loc[g]
    t_vals, n_vals = [], []
    for s, v in expr.items():
        st = st_map.get(s, "Unknown")
        if st == "Primary Tumor":      t_vals.append(np.log2(v + 1))
        elif st == "Solid Tissue Normal": n_vals.append(np.log2(v + 1))
    t_arr, n_arr = np.array(t_vals), np.array(n_vals)
    p_tn = stats.mannwhitneyu(t_arr, n_arr, alternative="two-sided").pvalue
    row.update(tumor_n=len(t_arr), normal_n=len(n_arr),
               log2FC_T_minus_N=t_arr.mean() - n_arr.mean(),
               wilcoxon_p_TvN=p_tn)
    tcga_rows.append(row)
tn_df = pd.DataFrame(tcga_rows)

# ----- OS -----
tumor_samples = [s for s in cln.loc[cln["sample_type"] == "Primary Tumor", "sample"]
                 if s in tpm.columns]
surv = cln[cln["sample"].isin(tumor_samples)].copy()
surv["event"] = (surv["vital_status"] == "Dead").astype(int)
surv["time"] = np.where(surv["event"] == 1,
                        surv["days_to_death"], surv["days_to_last_follow_up"])
surv = surv.dropna(subset=["time"])
surv = surv[surv["time"] > 0].copy()
print(f"  OS cohort: n={len(surv)} events={int(surv['event'].sum())}")

os_rows = []
for g in cands:
    if g not in tpm.index:
        os_rows.append({"gene": g, "os_in_data": False}); continue
    expr_g = tpm.loc[g, surv["sample"].values].values
    df = pd.DataFrame({"time": surv["time"].values,
                       "event": surv["event"].values,
                       "expr": np.log2(expr_g + 1)})
    med = df["expr"].median()
    high, low = df[df["expr"] >= med], df[df["expr"] < med]
    lr_p = logrank_test(high["time"], low["time"],
                        high["event"], low["event"]).p_value
    cph = CoxPHFitter(); cph.fit(df, duration_col="time", event_col="event")
    s = cph.summary.loc["expr"]
    os_rows.append({"gene": g, "os_in_data": True,
                    "km_logrank_p": lr_p,
                    "cox_HR": np.exp(s["coef"]),
                    "cox_HR_low": np.exp(s["coef lower 95%"]),
                    "cox_HR_high": np.exp(s["coef upper 95%"]),
                    "cox_p": s["p"]})
os_df = pd.DataFrame(os_rows)

# ----- merge -----
out = dep_df.merge(tn_df, on="gene", how="outer").merge(os_df, on="gene", how="outer")
out["is_priority"] = out["gene"].isin(PRIORITY)

# composite score: count of significant + favorable signals
def score(r):
    s = 0
    if pd.notna(r.get("luad_mean")) and r["luad_mean"] < -0.3: s += 1
    if pd.notna(r.get("mw_p_LUAD_lt_other")) and r["mw_p_LUAD_lt_other"] < 0.05: s += 1
    if pd.notna(r.get("log2FC_T_minus_N")) and r["log2FC_T_minus_N"] > 0.5 and r.get("wilcoxon_p_TvN", 1) < 0.05: s += 1
    if pd.notna(r.get("km_logrank_p")) and r["km_logrank_p"] < 0.05: s += 1
    if pd.notna(r.get("cox_p")) and r["cox_p"] < 0.05 and r.get("cox_HR", 1) > 1: s += 1
    return s
out["score"] = out.apply(score, axis=1)
out = out.sort_values(["is_priority", "score", "gene"], ascending=[False, False, True])
out.to_csv(OUT / "scan_27_summary.csv", index=False)

# ----- print compact ranked table -----
print("\n" + "=" * 110)
print(" Ranked candidate scan (=priority pick)")
print("=" * 110)
hdr = f"{'gene':<10} {'':<2} {'DM_LUAD':>8} {'DM_p':>8} {'log2FC':>7} {'TvN_p':>8} {'KM_p':>8} {'HR':>6} {'Cox_p':>8} {'score':>5}"
print(hdr); print("-" * len(hdr))
for _, r in out.iterrows():
    star = "" if r["is_priority"] else " "
    def f(x, fmt=".2f"):
        if pd.isna(x): return "  —"
        return format(x, fmt)
    print(f"{r['gene']:<10} {star:<2} {f(r.get('luad_mean')):>8} {f(r.get('mw_p_LUAD_lt_other'),'.1e'):>8} "
          f"{f(r.get('log2FC_T_minus_N')):>7} {f(r.get('wilcoxon_p_TvN'),'.1e'):>8} "
          f"{f(r.get('km_logrank_p'),'.1e'):>8} {f(r.get('cox_HR')):>6} {f(r.get('cox_p'),'.1e'):>8} {int(r['score']):>5}")
print(f"\nWrote: {OUT / 'scan_27_summary.csv'}")
print("DONE.")
