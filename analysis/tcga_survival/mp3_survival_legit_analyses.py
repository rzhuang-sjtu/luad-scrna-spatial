"""Legitimate prognostic analyses for MP3 (no data-dredging).

  1. Tertile KM (predefined 33/66 cutoffs)
  2. Quartile KM (Q1 vs Q4 only, max-contrast but predefined)
  3. Continuous Cox (univariate)
  4. Multivariate Cox (already in tcga_luad_mp_cox_multivariate.csv)
  5. Predicted-risk KM (median split of Cox linear predictor)
  6. Subgroup KM: MP3 split within Neutrophil-high half
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

OUT = Path("${WORK_ROOT}/luad_figures/fig3")
SIG = pd.read_csv(OUT / "tcga_luad_mp_ssgsea.csv.gz")
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin["event"] = (clin["vital_status"].str.lower() == "dead").astype(int)
clin["time"] = np.where(clin["event"] == 1,
                         clin["days_to_death"],
                         clin["days_to_last_follow_up"])
clin = clin.dropna(subset=["time"]); clin = clin[clin["time"] > 0]

SIG["case_id"] = SIG["sample_barcode"].str[:12]
df = SIG.merge(clin[["case_id","time","event","age_at_diagnosis","ajcc_stage"]],
                on="case_id", how="inner").drop_duplicates("case_id").reset_index(drop=True)
print(f"merged: {len(df)} patients, {df['event'].sum()} events")

# stage numeric
def stage_to_num(s):
    if not isinstance(s, str): return np.nan
    s = s.strip().lower().replace("stage ","")
    m = {"i":1,"ia":1,"ib":1,"ii":2,"iia":2,"iib":2,
         "iii":3,"iiia":3,"iiib":3,"iv":4,"iva":4,"ivb":4}
    return m.get(s, np.nan)
df["stage_num"] = df["ajcc_stage"].apply(stage_to_num)

results = []

# Helper
def km_from_groups(g, label):
    """log-rank statistic + p across all groups."""
    r = multivariate_logrank_test(df["time"], g, df["event"])
    return {"analysis": label,
            "groups": int(len(pd.Series(g).dropna().unique())),
            "n_total": int(len(df)),
            "n_events": int(df["event"].sum()),
            "chi2": float(r.test_statistic),
            "p": float(r.p_value)}

# 1. Tertile KM
mp3 = df["MP3"].to_numpy()
tert = pd.qcut(mp3, 3, labels=["T1_low","T2_mid","T3_high"])
results.append({**km_from_groups(tert, "MP3 tertile log-rank (T1 vs T2 vs T3)")})

# 2. Q1 vs Q4 KM
q = pd.qcut(mp3, 4, labels=["Q1","Q2","Q3","Q4"])
mask_q14 = q.isin(["Q1","Q4"])
r_q = logrank_test(df.loc[mask_q14 & (q == "Q4"), "time"],
                    df.loc[mask_q14 & (q == "Q1"), "time"],
                    event_observed_A=df.loc[mask_q14 & (q == "Q4"), "event"],
                    event_observed_B=df.loc[mask_q14 & (q == "Q1"), "event"])
results.append({"analysis": "MP3 Q1 vs Q4 log-rank",
                "groups": 2,
                "n_total": int(mask_q14.sum()),
                "n_events": int(df.loc[mask_q14, "event"].sum()),
                "chi2": float(r_q.test_statistic),
                "p": float(r_q.p_value)})

# 3. Continuous Cox (univariate, single test, no multiplicity)
cph = CoxPHFitter()
cph.fit(df[["MP3","time","event"]].dropna(),
        duration_col="time", event_col="event")
mp3_uni = cph.summary.loc["MP3"]
results.append({"analysis": "MP3 continuous Cox (univariate)",
                "groups": np.nan,
                "n_total": int(len(df)),
                "n_events": int(df["event"].sum()),
                "chi2": float(mp3_uni["z"] ** 2),
                "p": float(mp3_uni["p"]),
                "HR": float(mp3_uni["exp(coef)"]),
                "HR_lo": float(mp3_uni["exp(coef) lower 95%"]),
                "HR_hi": float(mp3_uni["exp(coef) upper 95%"])})

# 4. Multivariate Cox (MP1-4 + age + stage)
mv = df[["MP1","MP2","MP3","MP4","age_at_diagnosis","stage_num","time","event"]].dropna()
cph_mv = CoxPHFitter(); cph_mv.fit(mv, duration_col="time", event_col="event")
mp3_mv = cph_mv.summary.loc["MP3"]
results.append({"analysis": "MP3 multivariate Cox (adj. age, stage, MP1/2/4)",
                "groups": np.nan,
                "n_total": int(len(mv)),
                "n_events": int(mv["event"].sum()),
                "chi2": float(mp3_mv["z"] ** 2),
                "p": float(mp3_mv["p"]),
                "HR": float(mp3_mv["exp(coef)"]),
                "HR_lo": float(mp3_mv["exp(coef) lower 95%"]),
                "HR_hi": float(mp3_mv["exp(coef) upper 95%"])})

# 5. Predicted-risk KM from multivariate Cox
mv2 = mv.copy()
mv2["lp"] = cph_mv.predict_log_partial_hazard(mv2)
mv2["risk_grp"] = (mv2["lp"] > mv2["lp"].median()).map({True:"high", False:"low"})
r_lp = logrank_test(mv2.loc[mv2["risk_grp"]=="high","time"],
                    mv2.loc[mv2["risk_grp"]=="low","time"],
                    event_observed_A=mv2.loc[mv2["risk_grp"]=="high","event"],
                    event_observed_B=mv2.loc[mv2["risk_grp"]=="low","event"])
results.append({"analysis": "Predicted-risk KM (multi-Cox linear predictor median)",
                "groups": 2,
                "n_total": int(len(mv2)),
                "n_events": int(mv2["event"].sum()),
                "chi2": float(r_lp.test_statistic),
                "p": float(r_lp.p_value)})

# 6. Subgroup KM: MP3 high vs low within Neu-high half
neu = pd.read_csv(OUT / "tcga_neutrophil_scores.csv")
if "Neutrophil_core" in neu.columns:
    nc = "Neutrophil_core"
elif "neutrophil_score" in neu.columns:
    nc = "neutrophil_score"
else:
    nc = neu.columns[1]
neu["case_id"] = neu.iloc[:,0].astype(str).str[:12]
df_n = df.merge(neu[["case_id", nc]], on="case_id").drop_duplicates("case_id")
neu_hi = df_n[nc] >= df_n[nc].median()
sub = df_n[neu_hi].copy()
sub["mp3_grp"] = (sub["MP3"] >= sub["MP3"].median()).map({True:"hi", False:"lo"})
r_sub = logrank_test(sub.loc[sub["mp3_grp"]=="hi","time"],
                     sub.loc[sub["mp3_grp"]=="lo","time"],
                     event_observed_A=sub.loc[sub["mp3_grp"]=="hi","event"],
                     event_observed_B=sub.loc[sub["mp3_grp"]=="lo","event"])
results.append({"analysis": "MP3 high vs low KM, restricted to Neu-high half",
                "groups": 2,
                "n_total": int(len(sub)),
                "n_events": int(sub["event"].sum()),
                "chi2": float(r_sub.test_statistic),
                "p": float(r_sub.p_value)})

out_df = pd.DataFrame(results)
out_df.to_csv(OUT / "mp3_legit_survival_summary.csv", index=False)
print("\n=== MP3 prognostic evidence (no multiplicity) ===")
for _, r in out_df.iterrows():
    star = "***" if r["p"] < 1e-3 else "**" if r["p"] < 1e-2 else "*" if r["p"] < 0.05 else "ns"
    print(f"  {r['analysis']:60s}  p={r['p']:.4f}  {star}")
print("\nwrote", OUT / "mp3_legit_survival_summary.csv")
