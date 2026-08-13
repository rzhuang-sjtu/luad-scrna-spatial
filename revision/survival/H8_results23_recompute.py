#!/usr/bin/env python
"""H8 — Recompute manuscript §2.3 (Fig 3 survival) end-to-end on **one documented case set**.

Why recompute rather than only fix units: in that paragraph, multivariable HRs are per 1-SD,
whereas univariate HRs (4.94 / 0.19 / 15.71) were on the raw ssGSEA scale, mixed in one sentence,
while the sentence opens with "all hazard ratios are per 1-SD". Fixing units requires per-1-SD univariate values
on the same sample set; main-text n=458/188 is not reproducible from the current clinical table (gives 443/182),
so recompute every number in the paragraph on this one case set and document inclusion criteria in the text.

Output: results/results23_recompute/ for line-by-line replacement in the main-text revision checklist.
"""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

SIG = "${PROJECT_ROOT}/results/step28_tcga_combined_scores.csv.gz"
CLIN = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv"
OUT = "${PROJECT_ROOT}/results/results23_recompute"
MPS = ["MP1", "MP2", "MP3", "MP4"]
CLINCOV = ["age", "stage_late", "male"]


def build():
    sig = pd.read_csv(SIG)
    clin = pd.read_csv(CLIN)
    sig["patient"] = sig.sample_barcode.str[:12]
    clin["patient"] = clin["sample_barcode"].astype(str).str[:12]
    d = sig.merge(clin, on="patient", how="inner", suffixes=("", "_c"))
    d = d.drop_duplicates(subset="patient", keep="first").reset_index(drop=True)
    n_merge = len(d)
    d["event"] = (d.vital_status.astype(str).str.lower() == "dead").astype(int)
    d["time"] = np.where(d.event == 1,
                         pd.to_numeric(d.days_to_death, errors="coerce"),
                         pd.to_numeric(d.days_to_last_follow_up, errors="coerce"))
    n_notime = int(d.time.isna().sum())
    n_nonpos = int((d.time <= 0).sum())
    d = d[d.time.notna() & (d.time > 0)].copy()
    age = pd.to_numeric(d.age_at_diagnosis, errors="coerce")
    d["age"] = age / 365.25 if age.median() > 200 else age
    d["male"] = (d.gender.astype(str).str.lower() == "male").astype(int)
    st = (d.ajcc_stage.astype(str).str.upper()
          .str.replace("STAGE ", "", regex=False)
          .str.replace(r"[AB]$", "", regex=True).str.strip())
    d["stage_late"] = st.str.startswith(("III", "IV")).astype(int)
    d["stage_ord"] = st.map({"I": 1, "II": 2, "III": 3, "IV": 4})
    for m in MPS:
        d[f"z{m}"] = (d[m] - d[m].mean()) / d[m].std()
    d["z_MP3_minus_MP4"] = d.zMP3 - d.zMP4
    print(f"after merge {n_merge} cases; follow-up missing {n_notime}, non-positive {n_nonpos} →"
          f"analysis set {len(d)} cases, {int(d.event.sum())} events")
    return d


def cox(df, cols):
    x = df[["time", "event"] + cols].dropna()
    c = CoxPHFitter().fit(x, "time", "event")
    s = c.summary
    return c, x, {k: dict(HR=s.loc[k, "exp(coef)"],
                          lo=s.loc[k, "exp(coef) lower 95%"],
                          hi=s.loc[k, "exp(coef) upper 95%"],
                          p=s.loc[k, "p"]) for k in cols}


def main():
    import os
    os.makedirs(OUT, exist_ok=True)
    d = build()
    rows = []

    # ── Univariate: raw scale vs z-scored (p must match exactly; proves unit issue only)
    for m in MPS:
        _, x_raw, r_raw = cox(d, [m])
        _, x_z, r_z = cox(d, [f"z{m}"])
        rows.append(dict(块="univariate", 变量=m,
                         HR原始尺度=r_raw[m]["HR"], HR_per1SD=r_z[f"z{m}"]["HR"],
                         CI_low=r_z[f"z{m}"]["lo"], CI_high=r_z[f"z{m}"]["hi"],
                         p=r_z[f"z{m}"]["p"], n=len(x_z),
                         events=int(x_z.event.sum()), SD=d[m].std()))

    # ── Multivariable: four MPs + clinical (main-text model)
    cols = [f"z{m}" for m in MPS] + CLINCOV
    _, xm, rm = cox(d, cols)
    for m in MPS:
        rows.append(dict(块="multivariable (four MPs + clinical)", 变量=m,
                         HR原始尺度=np.nan, HR_per1SD=rm[f"z{m}"]["HR"],
                         CI_low=rm[f"z{m}"]["lo"], CI_high=rm[f"z{m}"]["hi"],
                         p=rm[f"z{m}"]["p"], n=len(xm),
                         events=int(xm.event.sum()), SD=d[m].std()))
    # Alternative: ordinal stage
    _, xo, ro = cox(d, [f"z{m}" for m in MPS] + ["age", "stage_ord", "male"])
    for m in MPS:
        rows.append(dict(块="multivariable (ordinal stage)", 变量=m, HR原始尺度=np.nan,
                         HR_per1SD=ro[f"z{m}"]["HR"], CI_low=ro[f"z{m}"]["lo"],
                         CI_high=ro[f"z{m}"]["hi"], p=ro[f"z{m}"]["p"],
                         n=len(xo), events=int(xo.event.sum()), SD=d[m].std()))

    # ── Composite score
    _, xc, rc = cox(d, ["z_MP3_minus_MP4"] + CLINCOV)
    rows.append(dict(块="multivariable (composite + clinical)", 变量="z(MP3)-z(MP4)",
                     HR原始尺度=np.nan, HR_per1SD=rc["z_MP3_minus_MP4"]["HR"],
                     CI_low=rc["z_MP3_minus_MP4"]["lo"],
                     CI_high=rc["z_MP3_minus_MP4"]["hi"],
                     p=rc["z_MP3_minus_MP4"]["p"], n=len(xc),
                     events=int(xc.event.sum()), SD=np.nan))
    tab = pd.DataFrame(rows)
    tab.to_csv(f"{OUT}/cox_table.csv", index=False)

    # ── Q1 vs Q4 KM
    km = []
    for m in MPS:
        q = d[m].quantile([.25, .75])
        a = d[d[m] <= q.iloc[0]]
        b = d[d[m] >= q.iloc[1]]
        r = logrank_test(a.time, b.time, a.event, b.event)
        km.append(dict(MP=m, n_Q1=len(a), n_Q4=len(b),
                       events_Q1=int(a.event.sum()), events_Q4=int(b.event.sum()),
                       logrank_p=r.p_value))
    q = d.z_MP3_minus_MP4.quantile([.25, .75])
    a = d[d.z_MP3_minus_MP4 <= q.iloc[0]]
    b = d[d.z_MP3_minus_MP4 >= q.iloc[1]]
    r = logrank_test(a.time, b.time, a.event, b.event)
    km.append(dict(MP="z(MP3)-z(MP4)", n_Q1=len(a), n_Q4=len(b),
                   events_Q1=int(a.event.sum()), events_Q4=int(b.event.sum()),
                   logrank_p=r.p_value))
    pd.DataFrame(km).to_csv(f"{OUT}/km_q1q4.csv", index=False)

    # ── ΔC-index: stage-only model vs stage + composite
    base = d[["time", "event"] + CLINCOV + ["z_MP3_minus_MP4"]].dropna()
    c0 = CoxPHFitter().fit(base[["time", "event"] + CLINCOV], "time", "event")
    c1 = CoxPHFitter().fit(base, "time", "event")
    ci0 = concordance_index(base.time, -c0.predict_partial_hazard(base), base.event)
    ci1 = concordance_index(base.time, -c1.predict_partial_hazard(base), base.event)
    lr = 2 * (c1.log_likelihood_ - c0.log_likelihood_)
    from scipy.stats import chi2
    lr_p = chi2.sf(lr, 1)
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(1000):
        s = base.sample(len(base), replace=True, random_state=int(rng.integers(1e9)))
        try:
            m0 = CoxPHFitter().fit(s[["time", "event"] + CLINCOV], "time", "event")
            m1 = CoxPHFitter().fit(s, "time", "event")
            boots.append(
                concordance_index(s.time, -m1.predict_partial_hazard(s), s.event) -
                concordance_index(s.time, -m0.predict_partial_hazard(s), s.event))
        except Exception:
            pass
    lo, hi = np.percentile(boots, [2.5, 97.5])
    pd.DataFrame([dict(C_clinical=ci0, C_plus_composite=ci1, delta=ci1 - ci0,
                       boot_mean=float(np.mean(boots)), boot_lo=lo, boot_hi=hi,
                       LR_chi2=lr, LR_p=lr_p, n=len(base),
                       events=int(base.event.sum()))]
                 ).to_csv(f"{OUT}/c_index.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n── Cox (all HRs per 1-SD; raw scale for comparison only)")
    print(tab.round(4).to_string(index=False))
    print("\n── Q1 vs Q4 KM")
    print(pd.DataFrame(km).round(5).to_string(index=False))
    print("\n── ΔC-index")
    print(pd.read_csv(f"{OUT}/c_index.csv").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
