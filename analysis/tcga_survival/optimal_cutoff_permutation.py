"""Permutation-correct the Fig 3D optimal-cutoff KM p-values."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = Path("${WORK_ROOT}/luad_figures/fig3")
SIG = pd.read_csv(OUT / "tcga_luad_mp_ssgsea.csv.gz")
CLIN_CSV = Path("${DATA_ROOT}/TCGA_LUAD_analysis/TCGA_LUAD_clinical.csv")
clin = pd.read_csv(CLIN_CSV)

# Match script 09_tcga_survival.py: build OS time + event
clin["event"] = (clin["vital_status"].str.lower() == "dead").astype(int)
clin["time"] = np.where(clin["event"] == 1,
                         clin["days_to_death"],
                         clin["days_to_last_follow_up"])
clin = clin.dropna(subset=["time"])
clin = clin[clin["time"] > 0]

# Match SIG sample_barcode to clinical case_id
SIG["case_id"] = SIG["sample_barcode"].str[:12]
clin["case_id"] = clin["case_id"] if "case_id" in clin.columns else clin.iloc[:, 0]
df = SIG.merge(clin[["case_id","time","event"]], on="case_id", how="inner")
df = df.drop_duplicates("case_id").reset_index(drop=True)
print(f"merged: {len(df)} patients")

PCTLS = np.arange(0.30, 0.71, 0.01)
N_PERM = 1000
rng = np.random.default_rng(2026)


def scan_max_chi2(score, time, event):
    """Scan cutoffs in 30-70 percentile, return the maximum log-rank chi2."""
    best = 0.0; best_p = 1.0; best_cut = np.nan
    for q in PCTLS:
        cut = np.quantile(score, q)
        hi = score >= cut
        lo = ~hi
        if hi.sum() < 10 or lo.sum() < 10:
            continue
        try:
            r = logrank_test(time[hi], time[lo],
                             event_observed_A=event[hi],
                             event_observed_B=event[lo])
            if r.test_statistic > best:
                best = r.test_statistic
                best_p = r.p_value
                best_cut = cut
        except Exception:
            continue
    return best, best_p, best_cut


rows = []
for mp in ("MP1", "MP2", "MP3", "MP4"):
    score = df[mp].to_numpy()
    time = df["time"].to_numpy(); event = df["event"].to_numpy()
    obs_chi2, obs_p, cut = scan_max_chi2(score, time, event)

    null = np.empty(N_PERM)
    for i in range(N_PERM):
        idx = rng.permutation(len(score))
        chi2_i, _, _ = scan_max_chi2(score, time[idx], event[idx])
        null[i] = chi2_i
    p_perm = (1 + np.sum(null >= obs_chi2)) / (1 + N_PERM)

    rows.append({
        "MP": mp,
        "cutpoint_score": float(cut),
        "log_rank_chi2": float(obs_chi2),
        "p_uncorrected": float(obs_p),
        "p_permutation": float(p_perm),
        "n_permutations": int(N_PERM),
    })
    print(f"  {mp}: chi2={obs_chi2:.3f}  raw_p={obs_p:.2e}  perm_p={p_perm:.4f}")

out = pd.DataFrame(rows)
out.to_csv(OUT / "tcga_optimal_cutoff_perm_corrected.csv", index=False)
print(f"\nwrote {OUT / 'tcga_optimal_cutoff_perm_corrected.csv'}")
print(out.to_string(index=False))
