"""Step 9: TCGA-LUAD survival analysis of the 4 LUAD MPs.

Pipeline:
  1. Load TCGA-LUAD TPM matrix; filter to Primary Tumor (-01A) samples.
  2. ssGSEA-score each sample against MP1-MP4 top-100 signatures (using
     the rank-based ssGSEA implementation in gseapy).
  3. Merge with clinical. Build OS = days_to_death (if dead) else
     days_to_last_follow_up; event = vital_status == 'Dead'.
  4. Per MP: median split → Kaplan-Meier + log-rank; also continuous Cox.
  5. Multivariate Cox: MP scores + age + stage + gender.

Outputs (${WORK_ROOT}/luad_figures/fig3/):
  - tcga_luad_mp_ssgsea.csv.gz
  - tcga_luad_mp_km_logrank.csv
  - tcga_luad_mp_cox_univariate.csv
  - tcga_luad_mp_cox_multivariate.csv
  - tcga_luad_km_curves_MP{1..4}.csv   (time, group, n_at_risk, surv_prob)
  - tcga_luad_survival_summary.md
"""

from __future__ import annotations
import os, sys, time, json
import numpy as np
import pandas as pd
from pathlib import Path

TPM_CSV = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv")
CLIN_CSV = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
SIG_CSV = Path.home() / "luad/results/step6_mp_signatures_top100.csv"
OUTDIR = Path("${WORK_ROOT}/luad_figures/fig3")
OUTDIR.mkdir(parents=True, exist_ok=True)

VALID_MPS = ["MP1", "MP2", "MP3", "MP4"]
TOPN = 100

print(f"[{time.strftime('%H:%M:%S')}] loading signatures ...")
sig = pd.read_csv(SIG_CSV)
sig = sig[(sig["MP"].isin(VALID_MPS)) & (sig["rank"] <= TOPN)].copy()
gene_sets = {mp: df["gene"].tolist() for mp, df in sig.groupby("MP")}
print("  sig sizes:", {k: len(v) for k, v in gene_sets.items()})

print(f"[{time.strftime('%H:%M:%S')}] loading clinical ...")
clin = pd.read_csv(CLIN_CSV)
print(f"  clinical rows: {len(clin)}")
# Keep primary tumors only
clin_pt = clin[clin["sample_type"] == "Primary Tumor"].copy()
print(f"  primary-tumor rows: {len(clin_pt)}")

print(f"[{time.strftime('%H:%M:%S')}] loading TPM (~202MB) ...")
tpm = pd.read_csv(TPM_CSV, index_col=0)
print(f"  TPM shape: {tpm.shape}")
# columns are sample barcodes; rows are gene symbols
# filter columns to primary tumor samples in clinical
pt_samples = [s for s in tpm.columns if s in set(clin_pt["sample_barcode"])]
print(f"  primary-tumor samples present in TPM: {len(pt_samples)}")
tpm = tpm[pt_samples]

# Collapse duplicate gene symbols by max (TPM)
if not tpm.index.is_unique:
    tpm = tpm.groupby(tpm.index).max()
    print(f"  collapsed duplicate genes → shape: {tpm.shape}")

# log2(TPM+1) for ssGSEA stability (gseapy ssgsea uses rank so log scale
# doesn't change result, but reduces memory pressure during rank)
print(f"[{time.strftime('%H:%M:%S')}] computing ssGSEA ...")
import gseapy as gp

expr_for_gsea = np.log2(tpm + 1.0).astype("float32")
# gseapy ssgsea expects genes as index, samples as columns
ss = gp.ssgsea(
    data=expr_for_gsea,
    gene_sets=gene_sets,
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    min_size=5,
    max_size=5000,
    permutation_num=0,
    seed=0,
    threads=8,
)

scores = ss.res2d.pivot_table(index="Name", columns="Term", values="NES").astype(float)
scores.index.name = "sample_barcode"
scores = scores[VALID_MPS]
print(f"  score shape: {scores.shape}")
scores.to_csv(OUTDIR / "tcga_luad_mp_ssgsea.csv.gz", compression="gzip")

# Merge with clinical
print(f"[{time.strftime('%H:%M:%S')}] merging clinical + scores ...")
df = clin_pt.set_index("sample_barcode").join(scores, how="inner").reset_index()
# Build survival variables
df["event"] = (df["vital_status"].str.strip().str.lower() == "dead").astype(int)
df["time"] = np.where(df["event"] == 1, df["days_to_death"], df["days_to_last_follow_up"])
# Clean
df = df[df["time"].notna() & (df["time"] > 0)].copy()
df["age"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce") / 365.25
df["stage_simple"] = (
    df["ajcc_stage"]
    .fillna("Unknown")
    .str.extract(r"(Stage [IV]+)", expand=False)
    .fillna("Unknown")
)
df["stage_num"] = df["stage_simple"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
df["gender_bin"] = (df["gender"].str.lower() == "male").astype(int)
print(f"  merged n={len(df)}  events={int(df['event'].sum())}")

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

lr_rows, cox_uni_rows = [], []
km_curves = []

for mp in VALID_MPS:
    med = df[mp].median()
    grp = np.where(df[mp] >= med, "High", "Low")
    dh = df[grp == "High"]
    dl = df[grp == "Low"]
    lr = logrank_test(dh["time"], dl["time"], dh["event"], dl["event"])
    lr_rows.append({
        "MP": mp, "median_cut": med, "n_high": len(dh), "n_low": len(dl),
        "events_high": int(dh["event"].sum()), "events_low": int(dl["event"].sum()),
        "logrank_chi2": lr.test_statistic, "logrank_p": lr.p_value,
    })

    # KM curve rows for plotting
    for label, sub in [("High", dh), ("Low", dl)]:
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], sub["event"], label=label)
        tab = kmf.survival_function_.reset_index()
        tab.columns = ["time", "surv_prob"]
        tab["group"] = label
        tab["MP"] = mp
        # append at-risk counts
        at_risk = [(sub["time"] >= t).sum() for t in tab["time"]]
        tab["n_at_risk"] = at_risk
        km_curves.append(tab)

    # Univariate continuous Cox
    cph = CoxPHFitter()
    sub_df = df[["time", "event", mp]].dropna()
    cph.fit(sub_df, duration_col="time", event_col="event")
    s = cph.summary.loc[mp]
    cox_uni_rows.append({
        "MP": mp, "coef": s["coef"], "HR": s["exp(coef)"],
        "HR_lo": s["exp(coef) lower 95%"], "HR_hi": s["exp(coef) upper 95%"],
        "p": s["p"], "n": len(sub_df),
    })

pd.DataFrame(lr_rows).to_csv(OUTDIR / "tcga_luad_mp_km_logrank.csv", index=False)
pd.DataFrame(cox_uni_rows).to_csv(OUTDIR / "tcga_luad_mp_cox_univariate.csv", index=False)
pd.concat(km_curves, ignore_index=True).to_csv(OUTDIR / "tcga_luad_mp_km_curves.csv.gz",
                                                 compression="gzip", index=False)

mv = df[["time", "event", "age", "stage_num", "gender_bin"] + VALID_MPS].dropna()
print(f"  multivariate Cox n={len(mv)} events={int(mv['event'].sum())}")

# Z-score MP covariates so HRs are per-SD
mv_z = mv.copy()
for mp in VALID_MPS:
    mv_z[mp] = (mv_z[mp] - mv_z[mp].mean()) / mv_z[mp].std()

cph_mv = CoxPHFitter()
cph_mv.fit(mv_z, duration_col="time", event_col="event")
mv_out = cph_mv.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].copy()
mv_out.columns = ["coef", "HR", "HR_lo", "HR_hi", "p"]
mv_out["concordance"] = cph_mv.concordance_index_
mv_out.to_csv(OUTDIR / "tcga_luad_mp_cox_multivariate.csv")

print("\n=== Univariate continuous Cox ===")
print(pd.DataFrame(cox_uni_rows).round({"coef": 3, "HR": 3, "HR_lo": 3, "HR_hi": 3, "p": 5}))
print("\n=== Log-rank (median split) ===")
print(pd.DataFrame(lr_rows).round({"logrank_chi2": 2, "logrank_p": 5}))
print("\n=== Multivariate Cox (MP scores z-scored) ===")
print(mv_out.round({"coef": 3, "HR": 3, "HR_lo": 3, "HR_hi": 3, "p": 5}))

# Summary markdown
md = OUTDIR / "tcga_luad_survival_summary.md"
with open(md, "w") as f:
    f.write("# TCGA-LUAD survival summary\n\n")
    f.write(f"- n samples (merged, time>0): {len(df)}\n")
    f.write(f"- events (deaths): {int(df['event'].sum())}\n\n")
    f.write("## Univariate continuous Cox (per-unit ssGSEA NES)\n\n")
    f.write(pd.DataFrame(cox_uni_rows).round(4).to_markdown(index=False) + "\n\n")
    f.write("## Median-split log-rank\n\n")
    f.write(pd.DataFrame(lr_rows).round(4).to_markdown(index=False) + "\n\n")
    f.write("## Multivariate Cox (MP z-scored, + age + stage + gender)\n\n")
    f.write(mv_out.round(4).to_markdown() + "\n\n")
    f.write(f"- concordance index: {cph_mv.concordance_index_:.4f}\n")

print(f"\n[{time.strftime('%H:%M:%S')}] done → {OUTDIR}")
