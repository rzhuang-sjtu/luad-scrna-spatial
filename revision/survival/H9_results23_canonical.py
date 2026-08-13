#!/usr/bin/env python
"""H9 — Reproduce the **original case set** for manuscript §2.3 survival and fill missing numbers on that set.

  · KM Q1-vs-Q4 p and 111 cases per arm  ← `mp_survival_full_panel.py` (patient-deduplicated, n=443)
  · univariate / multivariable Cox            ← `09_tcga_survival.py` (primary-tumour samples, n=458)
and within the same script univariate ran on the **raw ssGSEA scale**, multivariable on the **z-scored** scale,
while the text uniformly says "all hazard ratios are per 1-SD" — archive 3.4 unit inconsistency.

Strictly reproduce script 09 merge rules (clin keeps Primary Tumor rows only; join on sample_barcode;
time>0) and recompute on this single case set:
  1. Univariate per 1-SD for all four MPs (p must match published values exactly → unit conversion only)
  2. Published multivariable (should fully recover 1.43 / 0.72 / 1.02 / 0.82)
  3. MP3 nested-model ladder and composite score (needed for response letter item 9; same case set as main text)
Output: results/results23_canonical/
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

SIG = "${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz"
CLIN = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv"
OUT = "${PROJECT_ROOT}/results/results23_canonical"
MPS = ["MP1", "MP2", "MP3", "MP4"]


def build():
    sc = pd.read_csv(SIG)
    clin = pd.read_csv(CLIN)
    clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()
    df = clin_pt.set_index("sample_barcode").join(
        sc.set_index("sample_barcode"), how="inner").reset_index()
    df["event"] = (df.vital_status.str.strip().str.lower() == "dead").astype(int)
    df["time"] = np.where(df.event == 1, df.days_to_death, df.days_to_last_follow_up)
    df = df[df.time.notna() & (df.time > 0)].copy()
    df["age"] = pd.to_numeric(df.age_at_diagnosis, errors="coerce") / 365.25
    st = df.ajcc_stage.fillna("Unknown").str.extract(r"(Stage [IV]+)", expand=False)
    df["stage_num"] = st.fillna("Unknown").map(
        {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
    df["male"] = (df.gender.str.lower() == "male").astype(int)
    print(f"case set: n={len(df)}  events={int(df.event.sum())}"
          f"with stage {int(df.stage_num.notna().sum())}")
    return df


def fit(df, cols, keep):
    x = df[["time", "event"] + cols].dropna()
    s = CoxPHFitter().fit(x, "time", "event").summary
    return dict(HR=s.loc[keep, "exp(coef)"], lo=s.loc[keep, "exp(coef) lower 95%"],
                hi=s.loc[keep, "exp(coef) upper 95%"], p=s.loc[keep, "p"],
                n=len(x), events=int(x.event.sum()))


def main():
    os.makedirs(OUT, exist_ok=True)
    df = build()

    # ── 1. Univariate: raw scale vs per 1-SD (same samples, same model; unit only)
    uni = []
    for m in MPS:
        x = df[["time", "event", m]].dropna()
        sd = x[m].std()
        raw = fit(df, [m], m)
        z = x.copy()
        z[m] = (z[m] - z[m].mean()) / sd
        s = CoxPHFitter().fit(z, "time", "event").summary
        uni.append(dict(MP=m, SD=sd, HR_raw=raw["HR"], p_raw=raw["p"],
                        HR_per1SD=s.loc[m, "exp(coef)"],
                        CI_lo=s.loc[m, "exp(coef) lower 95%"],
                        CI_hi=s.loc[m, "exp(coef) upper 95%"],
                        p=s.loc[m, "p"], n=len(x), events=int(x.event.sum())))
    uni = pd.DataFrame(uni)
    uni.to_csv(f"{OUT}/univariate_per1SD.csv", index=False)

    # ── 2. Published multivariable (four MPs z-scored + age + ordinal stage + sex)
    mv = df[["time", "event", "age", "stage_num", "male"] + MPS].dropna()
    z = mv.copy()
    for m in MPS:
        z[m] = (z[m] - z[m].mean()) / z[m].std()
    s = CoxPHFitter().fit(z, "time", "event").summary
    mvt = s[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].copy()
    mvt.columns = ["HR", "CI_lo", "CI_hi", "p"]
    mvt["n"] = len(mv)
    mvt["events"] = int(mv.event.sum())
    mvt.to_csv(f"{OUT}/multivariate_published.csv")

    # ── 3. MP3 nested ladder + composite (all on the same z-scored case set)
    d = df.copy()
    for m in MPS:
        d["z" + m] = (d[m] - d[m].mean()) / d[m].std()
    d["comp"] = d.zMP3 - d.zMP4
    CLIN3 = ["age", "stage_num", "male"]
    ladder = [
        ("MP3 alone", ["zMP3"], "zMP3"),
        ("MP3 + age, stage, sex", ["zMP3"] + CLIN3, "zMP3"),
        ("MP3 + the other three MPs", ["z" + m for m in MPS], "zMP3"),
        ("MP3 + three MPs + clinical (published model)",
         ["z" + m for m in MPS] + CLIN3, "zMP3"),
        ("MP3 + MP4 only", ["zMP3", "zMP4"], "zMP3"),
        ("Composite z(MP3)-z(MP4), alone", ["comp"], "comp"),
        ("Composite + clinical", ["comp"] + CLIN3, "comp"),
    ]
    rows = [dict(model=name, **fit(d, cols, keep)) for name, cols, keep in ladder]
    lad = pd.DataFrame(rows)
    lad.to_csv(f"{OUT}/mp3_ladder.csv", index=False)

    # Rerun univariate on the multivariable subset to address case-exclusion as an explanation
    sub = d[d[["age", "stage_num", "male"]].notna().all(axis=1)]
    same = fit(sub, ["zMP3"], "zMP3")
    pd.DataFrame([same]).to_csv(f"{OUT}/mp3_univariate_on_mv_subset.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n── Univariate (published used the HR_raw column; text claims per 1-SD)")
    print(uni.round(4).to_string(index=False))
    print("\n── Multivariable (should recover published MP3 1.43 / MP4 0.72 / MP2 1.02 / MP1 0.82)")
    print(mvt.round(4).to_string())
    print("\n── MP3 nested ladder")
    print(lad.round(4).to_string(index=False))
    print(f"\n── MP3 univariate on the same multivariable subset (rules out case exclusion):"
          f"HR={same['HR']:.3f} ({same['lo']:.3f}–{same['hi']:.3f}) "
          f"p={same['p']:.4f} n={same['n']}/{same['events']}")


if __name__ == "__main__":
    main()
