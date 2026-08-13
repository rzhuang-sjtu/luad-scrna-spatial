"""Step 9c: Enhanced TCGA-LUAD survival — optimal cutoff, stage-stratified, risk score.

Follows Step 9 (median-split) with four additional analyses:
  B1. Optimal cutoff per MP (max log-rank chi2 over 25-75 percentile)
  B2. Stage-stratified median-split KM (Early=I+II, Late=III+IV)
  B3. Combined risk score = z(MP2) - z(MP4)
  B4. Stage-stratified univariate Cox per MP

Outputs → ${WORK_ROOT}/luad_figures/fig3/:
  tcga_optimal_cutoff_km.csv / tcga_optimal_cutoff_summary.csv
  tcga_stage_stratified_km.csv / tcga_stage_stratified_summary.csv
  tcga_stage_stratified_cox.csv
  tcga_risk_score_km.csv / tcga_risk_score_cox.csv
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

MP = Path("${WORK_ROOT}/luad_figures/fig3/tcga_luad_mp_ssgsea.csv.gz")
CLIN = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
FIG = Path("${WORK_ROOT}/luad_figures/fig3")
MPS = ["MP1", "MP2", "MP3", "MP4"]

def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def km_curve(df, group_col, label_fmt="{v}"):
    rows = []
    for g, sub in df.groupby(group_col, observed=True):
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], sub["event"], label=label_fmt.format(v=g))
        tab = kmf.survival_function_.reset_index()
        tab.columns = ["time", "surv_prob"]
        tab["group"] = g
        tab["n_at_risk"] = [(sub["time"] >= t).sum() for t in tab["time"]]
        tab["n_total"] = len(sub)
        tab["events"] = int(sub["event"].sum())
        rows.append(tab)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    t0 = time.time()
    log("loading MP ssGSEA + clinical")
    mp = pd.read_csv(MP, index_col=0)
    clin = pd.read_csv(CLIN)
    clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()

    df = clin_pt.set_index("sample_barcode").join(mp, how="inner").reset_index()
    df["event"] = (df["vital_status"].str.strip().str.lower() == "dead").astype(int)
    df["time"] = np.where(df["event"] == 1, df["days_to_death"], df["days_to_last_follow_up"])
    df = df[df["time"].notna() & (df["time"] > 0)].copy()

    df["stage_simple"] = df["ajcc_stage"].fillna("Unknown") \
        .str.extract(r"(Stage [IV]+)", expand=False).fillna("Unknown")
    df["stage_group"] = df["stage_simple"].map({
        "Stage I": "Early", "Stage II": "Early",
        "Stage III": "Late", "Stage IV": "Late",
    })
    df["age"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce") / 365.25
    df["gender_bin"] = (df["gender"].str.lower() == "male").astype(int)

    log(f"  n={len(df)} events={int(df['event'].sum())}")
    log(f"  stage_group counts:\n{df['stage_group'].value_counts(dropna=False).to_string()}")

    log("B1. optimal cutoff per MP")
    opt_km_rows, opt_sum_rows = [], []
    for mp_col in MPS:
        scores = df[mp_col].values
        q_lo, q_hi = np.quantile(scores, [0.25, 0.75])
        candidates = np.linspace(q_lo, q_hi, 51)
        best = None
        for c in candidates:
            grp = np.where(scores >= c, "High", "Low")
            if sum(grp == "High") < 20 or sum(grp == "Low") < 20:
                continue
            mh = grp == "High"; ml = ~mh
            lr = logrank_test(df.loc[mh, "time"], df.loc[ml, "time"],
                              df.loc[mh, "event"], df.loc[ml, "event"])
            if (best is None) or (lr.test_statistic > best["chi2"]):
                best = {"cutpoint": float(c), "chi2": float(lr.test_statistic),
                        "p": float(lr.p_value),
                        "n_high": int(mh.sum()), "n_low": int(ml.sum()),
                        "events_high": int(df.loc[mh, "event"].sum()),
                        "events_low": int(df.loc[ml, "event"].sum())}
        best["MP"] = mp_col
        best["score_q25"] = float(q_lo); best["score_q75"] = float(q_hi)
        opt_sum_rows.append(best)
        # KM curve at optimal cutpoint
        cp = best["cutpoint"]
        df_t = df.assign(_grp=np.where(df[mp_col] >= cp, "High", "Low"))
        km = km_curve(df_t, "_grp")
        km["MP"] = mp_col; km["cutpoint"] = cp
        opt_km_rows.append(km)

    pd.DataFrame(opt_sum_rows).to_csv(FIG/"tcga_optimal_cutoff_summary.csv", index=False)
    pd.concat(opt_km_rows, ignore_index=True).to_csv(FIG/"tcga_optimal_cutoff_km.csv", index=False)
    log("  optimal cutoff summary:")
    log(pd.DataFrame(opt_sum_rows).round(4).to_string(index=False))

    log("B2. stage-stratified KM (Early=I+II, Late=III+IV)")
    ss_km_rows, ss_sum_rows = [], []
    for stage in ["Early", "Late"]:
        sub_stage = df[df["stage_group"] == stage].copy()
        for mp_col in MPS:
            med = sub_stage[mp_col].median()
            grp = np.where(sub_stage[mp_col] >= med, "High", "Low")
            d = sub_stage.assign(_grp=grp)
            lr = logrank_test(d[grp=="High"]["time"], d[grp=="Low"]["time"],
                              d[grp=="High"]["event"], d[grp=="Low"]["event"])
            ss_sum_rows.append({
                "stage": stage, "MP": mp_col, "median_cut": float(med),
                "n_high": int((grp=="High").sum()), "n_low": int((grp=="Low").sum()),
                "events_high": int(d[grp=="High"]["event"].sum()),
                "events_low": int(d[grp=="Low"]["event"].sum()),
                "chi2": float(lr.test_statistic), "p": float(lr.p_value),
            })
            km = km_curve(d, "_grp")
            km["stage"] = stage; km["MP"] = mp_col
            ss_km_rows.append(km)

    pd.DataFrame(ss_sum_rows).to_csv(FIG/"tcga_stage_stratified_summary.csv", index=False)
    pd.concat(ss_km_rows, ignore_index=True).to_csv(FIG/"tcga_stage_stratified_km.csv", index=False)
    log("  stage-stratified KM summary:")
    log(pd.DataFrame(ss_sum_rows).round(4).to_string(index=False))

    log("B3. combined risk score = z(MP2) - z(MP4)")
    rs_df = df.copy()
    rs_df["risk_score"] = ((rs_df["MP2"] - rs_df["MP2"].mean()) / rs_df["MP2"].std()
                           - (rs_df["MP4"] - rs_df["MP4"].mean()) / rs_df["MP4"].std())

    # median split KM
    med = rs_df["risk_score"].median()
    rs_df["rs_grp_median"] = np.where(rs_df["risk_score"] >= med, "High", "Low")
    lr_med = logrank_test(
        rs_df[rs_df["rs_grp_median"]=="High"]["time"],
        rs_df[rs_df["rs_grp_median"]=="Low"]["time"],
        rs_df[rs_df["rs_grp_median"]=="High"]["event"],
        rs_df[rs_df["rs_grp_median"]=="Low"]["event"],
    )
    # optimal cutpoint on risk_score
    scores = rs_df["risk_score"].values
    q_lo, q_hi = np.quantile(scores, [0.25, 0.75])
    best = None
    for c in np.linspace(q_lo, q_hi, 51):
        grp = np.where(scores >= c, "High", "Low")
        if sum(grp=="High") < 20 or sum(grp=="Low") < 20: continue
        mh = grp=="High"; ml = ~mh
        lr = logrank_test(rs_df.loc[mh,"time"], rs_df.loc[ml,"time"],
                          rs_df.loc[mh,"event"], rs_df.loc[ml,"event"])
        if (best is None) or (lr.test_statistic > best["chi2"]):
            best = {"cutpoint": float(c), "chi2": float(lr.test_statistic),
                    "p": float(lr.p_value),
                    "n_high": int(mh.sum()), "n_low": int(ml.sum())}
    rs_df["rs_grp_optimal"] = np.where(rs_df["risk_score"] >= best["cutpoint"], "High", "Low")

    km_rows = []
    for label, col in [("median", "rs_grp_median"), ("optimal", "rs_grp_optimal")]:
        km = km_curve(rs_df, col)
        km["split"] = label
        km["cutpoint"] = med if label=="median" else best["cutpoint"]
        km_rows.append(km)
    pd.concat(km_rows, ignore_index=True).to_csv(FIG/"tcga_risk_score_km.csv", index=False)

    # univariate Cox on continuous risk_score + covariates
    cox_df = rs_df[["time","event","risk_score","age",
                    "stage_group","gender_bin"]].dropna().copy()
    cox_df["stage_num"] = cox_df["stage_group"].map({"Early":1, "Late":2})
    cox_df = cox_df.drop(columns=["stage_group"]).dropna()
    cox_df["risk_score"] = (cox_df["risk_score"] - cox_df["risk_score"].mean()) \
                            / cox_df["risk_score"].std()
    cph_uni = CoxPHFitter()
    cph_uni.fit(cox_df[["time","event","risk_score"]], duration_col="time",
                event_col="event")
    uni_row = cph_uni.summary.loc["risk_score"]

    cph_mv = CoxPHFitter()
    cph_mv.fit(cox_df, duration_col="time", event_col="event")
    mv_tab = cph_mv.summary[["coef","exp(coef)","exp(coef) lower 95%",
                              "exp(coef) upper 95%","p"]].copy()
    mv_tab.columns = ["coef","HR","HR_lo","HR_hi","p"]
    mv_tab["concordance"] = cph_mv.concordance_index_

    cox_out = pd.DataFrame([
        {"model": "univariate", "covariate": "risk_score",
         "HR": uni_row["exp(coef)"], "HR_lo": uni_row["exp(coef) lower 95%"],
         "HR_hi": uni_row["exp(coef) upper 95%"], "p": uni_row["p"],
         "concordance": cph_uni.concordance_index_}
    ])
    for cov, r in mv_tab.iterrows():
        cox_out = pd.concat([cox_out, pd.DataFrame([{
            "model": "multivariate", "covariate": cov,
            "HR": r["HR"], "HR_lo": r["HR_lo"], "HR_hi": r["HR_hi"],
            "p": r["p"], "concordance": r["concordance"],
        }])], ignore_index=True)

    # Also record log-rank for median + optimal splits
    rs_summary = pd.DataFrame([
        {"method": "median_split", "cutpoint": float(med),
         "chi2": float(lr_med.test_statistic), "p": float(lr_med.p_value),
         "n_high": int((rs_df["rs_grp_median"]=="High").sum()),
         "n_low":  int((rs_df["rs_grp_median"]=="Low").sum())},
        {"method": "optimal_cutoff", **best},
    ])
    rs_summary.to_csv(FIG/"tcga_risk_score_summary.csv", index=False)
    cox_out.to_csv(FIG/"tcga_risk_score_cox.csv", index=False)
    log(f"  risk_score median-split chi2={lr_med.test_statistic:.2f} p={lr_med.p_value:.4f}")
    log(f"  risk_score optimal cutpoint={best['cutpoint']:.3f} chi2={best['chi2']:.2f} p={best['p']:.4f}")
    log("  Cox (multivariate):")
    log(mv_tab.round(4).to_string())

    log("B4. stage-stratified univariate Cox per MP")
    cox_rows = []
    for stage in ["Early", "Late"]:
        sub = df[df["stage_group"] == stage].copy()
        if len(sub) < 30: continue
        for mp_col in MPS:
            ss = sub[["time","event", mp_col]].dropna()
            if len(ss) < 30: continue
            cph = CoxPHFitter()
            cph.fit(ss, duration_col="time", event_col="event")
            s = cph.summary.loc[mp_col]
            cox_rows.append({
                "stage": stage, "MP": mp_col,
                "n": len(ss), "events": int(ss["event"].sum()),
                "coef": s["coef"], "HR": s["exp(coef)"],
                "HR_lo": s["exp(coef) lower 95%"],
                "HR_hi": s["exp(coef) upper 95%"],
                "p": s["p"], "concordance": cph.concordance_index_,
            })
    pd.DataFrame(cox_rows).to_csv(FIG/"tcga_stage_stratified_cox.csv", index=False)
    log("  stage-stratified Cox:")
    log(pd.DataFrame(cox_rows).round(4).to_string(index=False))

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
