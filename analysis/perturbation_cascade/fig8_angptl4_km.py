"""Compute ANGPTL4 KM (Q1 vs Q4 + multivariate Cox) for TCGA-LUAD."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

OUT = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data")
TPM = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv", index_col=0)
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin["event"] = (clin["vital_status"].str.lower() == "dead").astype(int)
clin["time"] = np.where(clin["event"] == 1,
                         clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin.dropna(subset=["time"]); clin = clin[clin["time"] > 0]


def stage_to_num(s):
    if not isinstance(s, str): return np.nan
    s = s.strip().lower().replace("stage ","")
    return {"i":1,"ia":1,"ib":1,"ii":2,"iia":2,"iib":2,
            "iii":3,"iiia":3,"iiib":3,"iv":4,"iva":4,"ivb":4}.get(s, np.nan)
clin["stage_num"] = clin["ajcc_stage"].apply(stage_to_num)

GENE = "ANGPTL4"
if GENE not in TPM.index:
    raise SystemExit(f"{GENE} not in TPM matrix")

# log2(TPM+1) per sample (column = sample)
expr = np.log2(TPM.loc[GENE].astype(float) + 1)
expr.index.name = "sample_barcode"
df_e = expr.reset_index().rename(columns={GENE: "expr"})
df_e["case_id"] = df_e["sample_barcode"].str[:12]
df = df_e.merge(clin[["case_id","time","event","age_at_diagnosis","stage_num"]],
                 on="case_id").drop_duplicates("case_id").reset_index(drop=True)
print(f"merged: {len(df)} patients, {df['event'].sum()} events")

# Q1 vs Q4
q1 = df["expr"].quantile(0.25); q3 = df["expr"].quantile(0.75)
df["grp"] = np.where(df["expr"] >= q3, "High",
              np.where(df["expr"] <= q1, "Low", "drop"))
sub = df[df["grp"].isin(["High","Low"])]
r = logrank_test(sub.loc[sub.grp=="High","time"],
                 sub.loc[sub.grp=="Low","time"],
                 event_observed_A=sub.loc[sub.grp=="High","event"],
                 event_observed_B=sub.loc[sub.grp=="Low","event"])
n_hi = int((sub.grp=="High").sum()); e_hi = int(sub.loc[sub.grp=="High","event"].sum())
n_lo = int((sub.grp=="Low").sum());  e_lo = int(sub.loc[sub.grp=="Low","event"].sum())

# Multivariate Cox with Age + Stage
mv = df[["expr","age_at_diagnosis","stage_num","time","event"]].dropna()
cph = CoxPHFitter(); cph.fit(mv, duration_col="time", event_col="event")
s = cph.summary.loc["expr"]

stats_row = {
    "gene": GENE,
    "cutoff_low": float(q1), "cutoff_high": float(q3),
    "n_high": n_hi, "n_low": n_lo,
    "events_high": e_hi, "events_low": e_lo,
    "logrank_chi2": float(r.test_statistic),
    "logrank_p":    float(r.p_value),
    "cox_mv_HR":    float(s["exp(coef)"]),
    "cox_mv_HR_lo": float(s["exp(coef) lower 95%"]),
    "cox_mv_HR_hi": float(s["exp(coef) upper 95%"]),
    "cox_mv_p":     float(s["p"]),
}
pd.DataFrame([stats_row]).to_csv(OUT / "8Q_km_ANGPTL4_stats.csv", index=False)

# KM curve points for plotting
curve_rows = []
for grp_name in ("High","Low"):
    sl = sub[sub["grp"] == grp_name]
    kmf = KaplanMeierFitter(); kmf.fit(sl["time"], sl["event"])
    sf = kmf.survival_function_.reset_index()
    sf.columns = ["time","surv_prob"]
    sf["n_at_risk"] = [int((sl["time"] >= t).sum()) for t in sf["time"]]
    sf["gene"]  = GENE
    sf["group"] = grp_name
    curve_rows.append(sf)

# Append to existing 8OP_km_long.csv (which has SRSF9 + SEC61G)
existing = pd.read_csv(OUT / "8OP_km_long.csv")
new = pd.concat(curve_rows, ignore_index=True)
all_curves = pd.concat([existing, new], ignore_index=True)
all_curves.to_csv(OUT / "8OPQ_km_long.csv", index=False)

print(f"\nANGPTL4: Q4 (n={n_hi}, {e_hi}d) vs Q1 (n={n_lo}, {e_lo}d)")
print(f"  Q1/Q4 log-rank p = {r.p_value:.4f}")
print(f"  Multivariate Cox (adj. age + stage): HR = {s['exp(coef)']:.2f} "
      f"({s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}), p = {s['p']:.4f}")
print(f"\nwrote {OUT/'8Q_km_ANGPTL4_stats.csv'}")
print(f"wrote {OUT/'8OPQ_km_long.csv'}  ({len(all_curves)} total curve rows)")
