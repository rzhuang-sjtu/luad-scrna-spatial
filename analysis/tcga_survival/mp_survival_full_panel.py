"""All-MP legitimate prognostic analysis (Q1 vs Q4 + Cox uni + Cox multivar)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

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


def stage_to_num(s):
    if not isinstance(s, str): return np.nan
    s = s.strip().lower().replace("stage ","")
    return {"i":1,"ia":1,"ib":1,"ii":2,"iia":2,"iib":2,
            "iii":3,"iiia":3,"iiib":3,"iv":4,"iva":4,"ivb":4}.get(s, np.nan)
df["stage_num"] = df["ajcc_stage"].apply(stage_to_num)

# Multivariate Cox shared model (single fit) — all MPs + age + stage
mv_input = df[["MP1","MP2","MP3","MP4","age_at_diagnosis","stage_num","time","event"]].dropna()
cph_mv = CoxPHFitter(); cph_mv.fit(mv_input, duration_col="time", event_col="event")

rows = []
for mp in ("MP1","MP2","MP3","MP4"):
    score = df[mp].to_numpy()

    # 1) Q1 vs Q4 KM (predefined cutoffs)
    q = pd.qcut(score, 4, labels=["Q1","Q2","Q3","Q4"])
    m14 = q.isin(["Q1","Q4"])
    r_q = logrank_test(df.loc[m14 & (q == "Q4"), "time"],
                        df.loc[m14 & (q == "Q1"), "time"],
                        event_observed_A=df.loc[m14 & (q == "Q4"), "event"],
                        event_observed_B=df.loc[m14 & (q == "Q1"), "event"])
    n_q4 = int((m14 & (q == "Q4")).sum())
    n_q1 = int((m14 & (q == "Q1")).sum())
    e_q4 = int(df.loc[m14 & (q == "Q4"), "event"].sum())
    e_q1 = int(df.loc[m14 & (q == "Q1"), "event"].sum())

    # 2) Univariate continuous Cox
    cu = CoxPHFitter(); cu.fit(df[[mp,"time","event"]].dropna(),
                                duration_col="time", event_col="event")
    s_u = cu.summary.loc[mp]

    # 3) Multivariate (from shared fit)
    s_m = cph_mv.summary.loc[mp]

    rows.append({
        "MP": mp,
        "Q1_n": n_q1, "Q1_events": e_q1,
        "Q4_n": n_q4, "Q4_events": e_q4,
        "Q1Q4_chi2": float(r_q.test_statistic),
        "Q1Q4_p":    float(r_q.p_value),
        "Cox_uni_HR":   float(s_u["exp(coef)"]),
        "Cox_uni_HR_lo": float(s_u["exp(coef) lower 95%"]),
        "Cox_uni_HR_hi": float(s_u["exp(coef) upper 95%"]),
        "Cox_uni_p":    float(s_u["p"]),
        "Cox_mv_HR":    float(s_m["exp(coef)"]),
        "Cox_mv_HR_lo": float(s_m["exp(coef) lower 95%"]),
        "Cox_mv_HR_hi": float(s_m["exp(coef) upper 95%"]),
        "Cox_mv_p":     float(s_m["p"]),
    })

out = pd.DataFrame(rows)
out.to_csv(OUT / "mp_full_panel_survival_summary.csv", index=False)

print(f"merged: {len(df)} patients, {df['event'].sum()} events")
print()
print(f"{'MP':<5} {'Q1vQ4_p':>9} {'CoxUni_p':>9} {'CoxMV_p':>9}  "
      f"{'CoxUni_HR':>10}  {'CoxMV_HR':>10}")
for _, r in out.iterrows():
    sig_q = "***" if r["Q1Q4_p"]<1e-3 else "**" if r["Q1Q4_p"]<1e-2 else "*" if r["Q1Q4_p"]<5e-2 else "ns"
    sig_u = "***" if r["Cox_uni_p"]<1e-3 else "**" if r["Cox_uni_p"]<1e-2 else "*" if r["Cox_uni_p"]<5e-2 else "ns"
    sig_m = "***" if r["Cox_mv_p"]<1e-3 else "**" if r["Cox_mv_p"]<1e-2 else "*" if r["Cox_mv_p"]<5e-2 else "ns"
    print(f"{r['MP']:<5} {r['Q1Q4_p']:>9.3g}{sig_q:>3}  "
          f"{r['Cox_uni_p']:>9.3g}{sig_u:>3}  "
          f"{r['Cox_mv_p']:>9.3g}{sig_m:>3}  "
          f"{r['Cox_uni_HR']:>5.2f}({r['Cox_uni_HR_lo']:.2f}-{r['Cox_uni_HR_hi']:.2f})  "
          f"{r['Cox_mv_HR']:>5.2f}({r['Cox_mv_HR_lo']:.2f}-{r['Cox_mv_HR_hi']:.2f})")
print()
print(f"wrote {OUT/'mp_full_panel_survival_summary.csv'}")
