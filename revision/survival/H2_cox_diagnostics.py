#!/usr/bin/env python
"""
H2 — Proportional hazards assumption and tied event times for the Cox model (revision analysis).

Motivation: the handling of tied event times and the proportional-hazards assumption
for the Cox models were not stated anywhere.

Neither item was in the archive; this script supplies both:

1. **Tied event times**: report how many ties exist (multiple events on the same day),
   and compare HRs under Efron vs Breslow approximations — if nearly identical,
   the tie-handling choice does not affect conclusions.
2. **Proportional hazards**: global and per-covariate Schoenfeld residual tests.
   p < 0.05 indicates time-varying effect; proportional hazards does not hold for that covariate.
3. If any covariate violates PH, report remedies: time stratification or time-interaction results.

Output: results/cox_diagnostics/
"""
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

SIG = "${PROJECT_ROOT}/results/step28_tcga_combined_scores.csv.gz"
CLIN = "${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv"
OUT = "${PROJECT_ROOT}/results/cox_diagnostics"


def main():
    os.makedirs(OUT, exist_ok=True)
    sig = pd.read_csv(SIG)
    clin = pd.read_csv(CLIN)
    print(f"signature table {sig.shape}, clinical table {clin.shape}", flush=True)
    print(f"clinical columns: {[c for c in clin.columns][:20]}", flush=True)

    # Align sample barcodes (clinical is usually patient-level TCGA-XX-XXXX; signatures are sample-level …-01A)
    sig["patient"] = sig.sample_barcode.str[:12]
    pid = [c for c in clin.columns
           if c.lower() in ("bcr_patient_barcode", "submitter_id", "patient",
                            "case_submitter_id", "sample_barcode")]
    if not pid:
        raise SystemExit(f"patient ID column not found; clinical columns: {list(clin.columns)}")
    clin["patient"] = clin[pid[0]].astype(str).str[:12]

    d = sig.merge(clin, on="patient", how="inner")
    # A patient may have multiple samples/aliquots; collapse to patient level or the same person re-enters the risk set
    n0 = len(d)
    d = d.drop_duplicates(subset="patient", keep="first").reset_index(drop=True)
    print(f"after merge {n0} rows → patient-level {len(d)} cases", flush=True)

    def pick(*names):
        for n in names:
            for c in d.columns:
                if c.lower() == n:
                    return c
        return None

    c_vital = pick("vital_status")
    c_death = pick("days_to_death")
    c_fu = pick("days_to_last_follow_up", "days_to_last_followup")
    c_stage = pick("ajcc_stage", "ajcc_pathologic_stage", "tumor_stage", "stage")
    c_age = pick("age_at_diagnosis", "age_at_index",
                 "age_at_initial_pathologic_diagnosis", "age")
    c_sex = pick("gender", "sex")
    print(f"using columns: vital={c_vital} death={c_death} fu={c_fu}"
          f"stage={c_stage} age={c_age} sex={c_sex}", flush=True)

    d["event"] = (d[c_vital].astype(str).str.lower() == "dead").astype(int)
    d["time"] = np.where(d.event == 1,
                         pd.to_numeric(d[c_death], errors="coerce"),
                         pd.to_numeric(d[c_fu], errors="coerce"))
    d = d[d.time.notna() & (d.time > 0)].copy()

    # covariates
    age = pd.to_numeric(d[c_age], errors="coerce")
    # TCGA age_at_diagnosis is often in days; convert to years if >200
    d["age"] = age / 365.25 if age.median() > 200 else age
    d["male"] = (d[c_sex].astype(str).str.lower() == "male").astype(int)
    st = d[c_stage].astype(str).str.upper().str.replace("STAGE ", "", regex=False)
    d["stage_late"] = st.str.startswith(("III", "IV")).astype(int)
    for m in ["MP1", "MP2", "MP3", "MP4"]:
        d[f"z{m}"] = (d[m] - d[m].mean()) / d[m].std()

    d = d.dropna(subset=["age", "time", "event"]).copy()
    print(f"\nanalysable {len(d)} cases, {int(d.event.sum())} events", flush=True)

    et = d.loc[d.event == 1, "time"]
    vc = et.value_counts()
    n_tied_times = int((vc > 1).sum())
    n_in_ties = int(vc[vc > 1].sum())
    print(f"\n【Tied event times】")
    print(f"distinct event times {et.nunique()}, of which {n_tied_times} are tied")
    print(f"events involved {n_in_ties} / {int(d.event.sum())}"
          f"({100*n_in_ties/d.event.sum():.1f}%), max tie size {int(vc.max())}")

    covs = ["zMP1", "zMP2", "zMP3", "zMP4", "age", "male", "stage_late"]
    sub = d[["time", "event"] + covs].dropna()
    rows = []
    for tie in ["efron", "breslow"]:
        cph = CoxPHFitter()
        cph.fit(sub, "time", "event", strata=None,
                fit_options={"step_size": 0.5}) if False else None
        cph = CoxPHFitter()
        # lifelines defaults to Efron; Breslow is not available via a ties argument — compare manually:
        # lifelines implements only Efron; use statsmodels PHReg for Breslow comparison
        break
    from statsmodels.duration.hazard_regression import PHReg
    print("\n【Effect of tie approximation on results】")
    print(f"{'covariate':<12}{'Efron HR':>10}{'Breslow HR':>12}{'diff%':>8}")
    hrs = {}
    for method in ["efron", "breslow"]:
        m = PHReg(sub.time.values, sub[covs].values, status=sub.event.values,
                  ties=method).fit()
        hrs[method] = np.exp(m.params)
    for i, c in enumerate(covs):
        e, b = hrs["efron"][i], hrs["breslow"][i]
        print(f"  {c:<12}{e:>10.4f}{b:>12.4f}{100*(b-e)/e:>7.2f}%")
    pd.DataFrame({"covariate": covs, "HR_efron": hrs["efron"],
                  "HR_breslow": hrs["breslow"]}).to_csv(
        f"{OUT}/tie_methods.csv", index=False)

    cph = CoxPHFitter().fit(sub, "time", "event")
    print("\n【Model (Efron, lifelines)】")
    print(cph.summary[["exp(coef)", "exp(coef) lower 95%",
                       "exp(coef) upper 95%", "p"]].to_string(
        float_format=lambda x: f"{x:.4g}"))
    cph.summary.to_csv(f"{OUT}/cox_model.csv")

    print("\n【Proportional hazards: Schoenfeld residual test】")
    print("p < 0.05 indicates a time-varying effect; proportional hazards does not hold")
    res = proportional_hazard_test(cph, sub, time_transform=["rank", "km"])
    print(res.summary.to_string(float_format=lambda x: f"{x:.4g}"))
    res.summary.to_csv(f"{OUT}/ph_test.csv")

    bad = res.summary[res.summary.p < 0.05]
    if len(bad):
        print(f"\n   covariates violating proportional hazards: {sorted(set(bad.index.get_level_values(0)))}")
        print("remedy: stratified Cox on violators; re-estimate remaining covariates")
        strat = sorted(set(bad.index.get_level_values(0)))
        strat = [s for s in strat if sub[s].nunique() <= 5]
        if strat:
            cph2 = CoxPHFitter().fit(sub, "time", "event", strata=strat)
            print(f"after stratifying by {strat}:")
            print(cph2.summary[["exp(coef)", "exp(coef) lower 95%",
                                "exp(coef) upper 95%", "p"]].to_string(
                float_format=lambda x: f"{x:.4g}"))
            cph2.summary.to_csv(f"{OUT}/cox_stratified.csv")
        else:
            print("violators are continuous; need time-interaction terms; see ph_test.csv")
    else:
        print("\n   no covariate violates proportional hazards (global and per-term p ≥ 0.05)")
    print(f"\nWrote {OUT}/", flush=True)


if __name__ == "__main__":
    main()
