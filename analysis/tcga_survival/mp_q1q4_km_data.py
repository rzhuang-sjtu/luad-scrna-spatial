"""Compute Q1 vs Q4 KM curves for all 4 MPs and write CSV for R plotting."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = Path("${WORK_ROOT}/luad_figures/fig3")
SIG = pd.read_csv(OUT / "tcga_luad_mp_ssgsea.csv.gz")
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin["event"] = (clin["vital_status"].str.lower() == "dead").astype(int)
clin["time"] = np.where(clin["event"] == 1,
                         clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin.dropna(subset=["time"]); clin = clin[clin["time"] > 0]
SIG["case_id"] = SIG["sample_barcode"].str[:12]
df = SIG.merge(clin[["case_id","time","event"]], on="case_id", how="inner")\
        .drop_duplicates("case_id").reset_index(drop=True)

curve_rows = []
stat_rows  = []
for mp in ("MP1","MP2","MP3","MP4"):
    s = df[mp].to_numpy()
    q1 = np.quantile(s, 0.25); q3 = np.quantile(s, 0.75)
    grp = pd.Series(np.where(s >= q3, "High",
                     np.where(s <= q1, "Low", "drop")),
                    index=df.index)
    sub = df[grp.isin(["High","Low"])].copy()
    sub["grp"] = grp[sub.index]
    r = logrank_test(sub.loc[sub.grp=="High","time"],
                     sub.loc[sub.grp=="Low","time"],
                     event_observed_A=sub.loc[sub.grp=="High","event"],
                     event_observed_B=sub.loc[sub.grp=="Low","event"])
    n_hi = int((sub.grp=="High").sum()); e_hi = int(sub.loc[sub.grp=="High","event"].sum())
    n_lo = int((sub.grp=="Low").sum());  e_lo = int(sub.loc[sub.grp=="Low","event"].sum())
    stat_rows.append({"MP": mp, "n_high": n_hi, "events_high": e_hi,
                      "n_low": n_lo, "events_low": e_lo,
                      "log_rank_chi2": float(r.test_statistic),
                      "log_rank_p":    float(r.p_value)})

    for grp_name in ("High","Low"):
        sl = sub[sub.grp == grp_name]
        kmf = KaplanMeierFitter()
        kmf.fit(sl["time"], sl["event"])
        sf = kmf.survival_function_.reset_index()
        sf.columns = ["time","surv_prob"]
        n_remain = []
        for t in sf["time"]:
            n_remain.append(int((sl["time"] >= t).sum()))
        sf["n_at_risk"] = n_remain
        sf["MP"] = mp; sf["group"] = grp_name
        curve_rows.append(sf)

curves = pd.concat(curve_rows, ignore_index=True)
stats = pd.DataFrame(stat_rows)
curves.to_csv(OUT / "tcga_luad_mp_q1q4_km_curves.csv.gz",
              index=False, compression="gzip")
stats.to_csv(OUT / "tcga_luad_mp_q1q4_km_logrank.csv", index=False)
print(f"merged: {len(df)} patients, {df['event'].sum()} events")
print(f"\nQ1 (bottom 25%) vs Q4 (top 25%) log-rank per MP:")
for _, r in stats.iterrows():
    s = "***" if r['log_rank_p']<1e-3 else "**" if r['log_rank_p']<1e-2 else "*" if r['log_rank_p']<5e-2 else "ns"
    print(f"  {r['MP']:5s}  Q4 n={r['n_high']:3d}/{r['events_high']:3d}d  "
          f"Q1 n={r['n_low']:3d}/{r['events_low']:3d}d  p={r['log_rank_p']:.4f} {s}")
print(f"\nwrote {OUT/'tcga_luad_mp_q1q4_km_curves.csv.gz'}")
print(f"wrote {OUT/'tcga_luad_mp_q1q4_km_logrank.csv'}")
