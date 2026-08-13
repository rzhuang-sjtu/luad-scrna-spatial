"""ANGPTL4 KM with median split (matches SRSF9/SEC61G style in 8M/8N)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = Path("${WORK_ROOT}/luad_figures/fig8/v2_500/data")
TPM = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_TPM_matrix.csv", index_col=0)
clin = pd.read_csv("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin["event"] = (clin["vital_status"].str.lower() == "dead").astype(int)
clin["time"] = np.where(clin["event"] == 1, clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin.dropna(subset=["time"]); clin = clin[clin["time"] > 0]

GENE = "ANGPTL4"
expr = np.log2(TPM.loc[GENE].astype(float) + 1).reset_index().rename(
    columns={"index":"sample_barcode", GENE:"expr"})
expr["case_id"] = expr["sample_barcode"].str[:12]
df = expr.merge(clin[["case_id","time","event"]], on="case_id")\
        .drop_duplicates("case_id").reset_index(drop=True)

cut = df["expr"].median()
df["group"] = np.where(df["expr"] >= cut, "High", "Low")
n_hi = int((df["group"]=="High").sum()); n_lo = int((df["group"]=="Low").sum())
e_hi = int(df.loc[df.group=="High","event"].sum())
e_lo = int(df.loc[df.group=="Low","event"].sum())
r = logrank_test(df.loc[df.group=="High","time"],
                 df.loc[df.group=="Low","time"],
                 event_observed_A=df.loc[df.group=="High","event"],
                 event_observed_B=df.loc[df.group=="Low","event"])

stats_row = {"gene": GENE, "split": "median", "cutoff_logTPM": float(cut),
             "n_high": n_hi, "n_low": n_lo,
             "events_high": e_hi, "events_low": e_lo,
             "logrank_chi2": float(r.test_statistic),
             "logrank_p":    float(r.p_value)}
pd.DataFrame([stats_row]).to_csv(OUT / "8O_km_ANGPTL4_stats.csv", index=False)

# KM curves: replace ANGPTL4 entries in the long table with median-split
existing = pd.read_csv(OUT / "8OPQ_km_long.csv")
existing = existing[existing["gene"] != GENE].copy()
curve_rows = []
for grp in ("High","Low"):
    sl = df[df["group"] == grp]
    kmf = KaplanMeierFitter(); kmf.fit(sl["time"], sl["event"])
    sf = kmf.survival_function_.reset_index()
    sf.columns = ["time","surv_prob"]
    sf["n_at_risk"] = [int((sl["time"] >= t).sum()) for t in sf["time"]]
    sf["gene"] = GENE; sf["group"] = grp
    curve_rows.append(sf)
all_curves = pd.concat([existing, pd.concat(curve_rows, ignore_index=True)],
                       ignore_index=True)
all_curves.to_csv(OUT / "8OPQ_km_long.csv", index=False)

print(f"ANGPTL4 (median split): n_high={n_hi}/{e_hi}d, n_low={n_lo}/{e_lo}d, "
      f"log-rank p = {r.p_value:.4f}")
print(f"wrote {OUT/'8O_km_ANGPTL4_stats.csv'}")
