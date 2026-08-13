#!/usr/bin/env python
"""
H5 — Independent recheck of archive 3.4 (inconsistent HR units) and 3.5 (MP3 "stage-independent" wording).

Both are marked must-fix for the response letter; recompute rather than cite the 08-04 record alone.

3.4 claim: main text states all HRs are "per 1-SD increase in z-scored MP score",
but the univariate column was run on the raw ssGSEA scale.
  Check: same samples, univariate Cox on raw scale vs z-scored,
  see which set reproduces the reported numbers. p-values should match exactly (scale does not change the test);
  only HR changes — key evidence that this is a scale issue, not a different model.

3.5 claim: MP3 significance comes from mutual adjustment with the other three MPs (suppression),
not stage-independence; adding stage alone actually weakens it.
  Check: four nested models + collinearity diagnostics (correlation matrix, VIF) + alternative stage coding.
  If conclusions flip under alternative stage coding, the original claim is unstable.

Output: results/survival_recheck/
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm

SIG = "${PROJECT_ROOT}/results/step28_tcga_combined_scores.csv.gz"
CLIN = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv"
OUT = "${PROJECT_ROOT}/results/survival_recheck"
MPS = ["MP1", "MP2", "MP3", "MP4"]


def build():
    sig = pd.read_csv(SIG)
    clin = pd.read_csv(CLIN)
    sig["patient"] = sig.sample_barcode.str[:12]
    clin["patient"] = clin["sample_barcode"].astype(str).str[:12]
    d = sig.merge(clin, on="patient", how="inner", suffixes=("", "_clin"))
    d = d.drop_duplicates(subset="patient", keep="first").reset_index(drop=True)
    d["event"] = (d.vital_status.astype(str).str.lower() == "dead").astype(int)
    d["time"] = np.where(d.event == 1,
                         pd.to_numeric(d.days_to_death, errors="coerce"),
                         pd.to_numeric(d.days_to_last_follow_up, errors="coerce"))
    d = d[d.time.notna() & (d.time > 0)].copy()
    age = pd.to_numeric(d.age_at_diagnosis, errors="coerce")
    d["age"] = age / 365.25 if age.median() > 200 else age
    d["male"] = (d.gender.astype(str).str.lower() == "male").astype(int)
    st = d.ajcc_stage.astype(str).str.upper().str.replace("STAGE ", "", regex=False)
    st = st.str.replace(r"[AB]$", "", regex=True).str.strip()
    d["stage_late"] = st.str.startswith(("III", "IV")).astype(int)
    d["stage_ord"] = st.map({"I": 1, "II": 2, "III": 3, "IV": 4})
    d["stage_known"] = d.stage_ord.notna()
    for m in MPS:
        d[f"z{m}"] = (d[m] - d[m].mean()) / d[m].std()
    d["z_MP3_minus_MP4"] = d.zMP3 - d.zMP4
    return d


def cox(df, cols, dur="time", ev="event"):
    c = CoxPHFitter().fit(df[[dur, ev] + cols].dropna(), dur, ev)
    return c


def main():
    os.makedirs(OUT, exist_ok=True)
    d = build()
    print(f"after merge, analysable {len(d)} cases, {int(d.event.sum())} events", flush=True)

    print("\n" + "=" * 68)
    print("3.4 recheck: univariate Cox, raw ssGSEA scale vs z-scored")
    print("=" * 68)
    rows = []
    for m in MPS:
        craw = cox(d, [m])
        cz = cox(d, [f"z{m}"])
        rows.append(dict(
            MP=m,
            HR_原始尺度=float(np.exp(craw.params_[m])),
            HR_z标准化=float(np.exp(cz.params_[f"z{m}"])),
            CI_low=float(np.exp(cz.confidence_intervals_.loc[f"z{m}"].iloc[0])),
            CI_high=float(np.exp(cz.confidence_intervals_.loc[f"z{m}"].iloc[1])),
            p_原始=float(craw.summary.loc[m, "p"]),
            p_z=float(cz.summary.loc[f"z{m}", "p"]),
            SD=float(d[m].std()),
            n=int(len(d)), events=int(d.event.sum())))
    U = pd.DataFrame(rows)
    U.to_csv(f"{OUT}/univariate_scale.csv", index=False)
    print(U.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    same_p = np.allclose(U.p_原始, U.p_z, rtol=1e-6)
    print(f"\n  p-values identical across scales: {same_p}")
    print("→ same p with different HR is the criterion for 'same model, different scale';"
          "if the model itself differed, p would change too.")

    print("\n" + "=" * 68)
    print("3.5 recheck: source of MP3 significance")
    print("=" * 68)
    ds = d[d.stage_known].copy()
    print(f"stage-available subset n={len(ds)}, events {int(ds.event.sum())}")

    specs = [
        ("MP3 alone", ["zMP3"]),
        ("MP3 + clinical (age/binary stage/sex)", ["zMP3", "age", "stage_late", "male"]),
        ("MP3 + other three MPs (no clinical)", ["zMP1", "zMP2", "zMP3", "zMP4"]),
        ("MP3 + other three MPs + clinical (main-text model)",
         ["zMP1", "zMP2", "zMP3", "zMP4", "age", "stage_late", "male"]),
        ("alternative: stage as ordinal I–IV",
         ["zMP1", "zMP2", "zMP3", "zMP4", "age", "stage_ord", "male"]),
        ("alternative: MP3 + MP4 only", ["zMP3", "zMP4"]),
        ("composite z(MP3)−z(MP4) alone", ["z_MP3_minus_MP4"]),
        ("composite + clinical", ["z_MP3_minus_MP4", "age", "stage_late", "male"]),
    ]
    rows = []
    for name, cols in specs:
        c = cox(ds, cols)
        key = "zMP3" if "zMP3" in cols else "z_MP3_minus_MP4"
        ci = c.confidence_intervals_.loc[key]
        rows.append(dict(模型=name, 关注变量=key,
                         HR=float(np.exp(c.params_[key])),
                         CI_low=float(np.exp(ci.iloc[0])),
                         CI_high=float(np.exp(ci.iloc[1])),
                         p=float(c.summary.loc[key, "p"]),
                         n=int(c._n_examples), events=int(c.event_observed.sum())))
    M = pd.DataFrame(rows)
    M.to_csv(f"{OUT}/mp3_model_ladder.csv", index=False)
    print(M.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # Collinearity and suppression diagnostics
    print("\n【Inter-MP Spearman correlations and MP3 vs stage】")
    cm = ds[[f"z{m}" for m in MPS]].corr(method="spearman")
    print(cm.to_string(float_format=lambda x: f"{x:+.3f}"))
    print(f"\n  MP3 vs ordinal stage Spearman ="
          f"{ds[['zMP3','stage_ord']].corr(method='spearman').iloc[0,1]:+.4f}")

    X = ds[[f"z{m}" for m in MPS] + ["age", "stage_late", "male"]].dropna()
    X = sm.add_constant(X)
    vif = pd.DataFrame({"variable": X.columns[1:],
                        "VIF": [variance_inflation_factor(X.values, i)
                                for i in range(1, X.shape[1])]})
    print("\n【Variance inflation factor】VIF > 5 suggests substantial collinearity")
    print(vif.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    vif.to_csv(f"{OUT}/vif.csv", index=False)
    print(f"\nWrote {OUT}/", flush=True)


if __name__ == "__main__":
    main()
