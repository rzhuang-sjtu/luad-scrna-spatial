"""
Step 9e: ROI vs non-ROI per-spot significance testing.

Reads ${DATA_ROOT}/ST/results/r_data/per_section/*.csv, pools spots within each
cohort (E-MTAB-13530, Takano/Okamura FF), then for every metric column
(cell-type abundance, MP score, PROGENy pathway, gene expression) runs a
two-sided Mann-Whitney U test of values inside ROI vs outside ROI.

Outputs:
  ${DATA_ROOT}/ST/results/r_data/roi_vs_nonroi_stats_with_pvalues.csv
      Per-section, per-metric raw stats + spot counts.
  ${DATA_ROOT}/ST/results/r_data/roi_vs_nonroi_aggregate_pvalues.csv
      Cohort-pooled Mann-Whitney + BH-FDR.

The aggregate file is the one to cite in the paper.
"""
from __future__ import annotations
import os, glob
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PER = Path("${DATA_ROOT}/ST/results/r_data/per_section")
OUT = Path("${DATA_ROOT}/ST/results/r_data")

# Group columns by metric type (matching S10J factor levels).
def classify(col: str) -> str | None:
    if col.startswith("ct_"):       return "celltype"
    if col.startswith("progeny_"):  return "obs"
    if col.startswith("gex_"):      return "gene"
    if col.endswith("_score") and col.startswith("MP"): return "obs"
    return None


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR — returns adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    # enforce monotonicity from the back
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def stars(p: float) -> str:
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


def main():
    files = sorted(glob.glob(str(PER / "*.csv")))
    print(f"per-section CSVs: {len(files)}")

    # Pool by cohort prefix
    by_cohort: dict[str, list[pd.DataFrame]] = {}
    for f in files:
        name = Path(f).stem            # e.g. EMTAB13530__P10_T1
        cohort = name.split("__")[0]
        df = pd.read_csv(f)
        df["__sample"] = name.split("__", 1)[1]
        df["__cohort"] = cohort
        by_cohort.setdefault(cohort, []).append(df)

    per_section_rows = []
    aggregate_rows   = []

    for cohort, lst in by_cohort.items():
        print(f"\n=== cohort: {cohort} ({len(lst)} sections) ===")
        cohort_df = pd.concat(lst, axis=0, ignore_index=True)
        if "roi" not in cohort_df.columns:
            print(f"  [SKIP] no roi column"); continue
        # Coerce roi to boolean even if stored as 0/1/"True"/"False"
        roi_raw = cohort_df["roi"]
        if roi_raw.dtype == object:
            roi = roi_raw.astype(str).str.lower().isin({"true","1","t"})
        else:
            roi = roi_raw.astype(bool)
        n_roi    = int(roi.sum())
        n_nonroi = int((~roi).sum())
        print(f"  spots: ROI={n_roi}, non-ROI={n_nonroi}")
        if n_roi < 5 or n_nonroi < 5:
            print(f"  [SKIP] too few spots"); continue

        metric_cols = [c for c in cohort_df.columns if classify(c) is not None]
        print(f"  metrics: {len(metric_cols)}")

        cohort_pvals = []
        cohort_idx   = []
        for c in metric_cols:
            v = pd.to_numeric(cohort_df[c], errors="coerce").to_numpy()
            mask = ~np.isnan(v)
            if mask.sum() < 10:
                continue
            v_in  = v[mask & roi.to_numpy()]
            v_out = v[mask & (~roi.to_numpy())]
            if len(v_in) < 5 or len(v_out) < 5:
                continue
            try:
                u, p = mannwhitneyu(v_in, v_out, alternative="two-sided")
            except ValueError:
                continue
            mean_in  = float(np.mean(v_in))
            mean_out = float(np.mean(v_out))
            cohort_pvals.append(p)
            cohort_idx.append(c)
            aggregate_rows.append({
                "cohort": cohort,
                "metric": c,
                "type": classify(c),
                "mean_roi": mean_in,
                "mean_nonroi": mean_out,
                "delta": mean_in - mean_out,
                "U_stat": float(u),
                "p_raw": float(p),
                "n_roi": int(len(v_in)),
                "n_nonroi": int(len(v_out)),
            })

        # Per-section Mann-Whitney + per-section BH-FDR (across metrics).
        # Lets us paint a per-sample × per-metric consistency heatmap.
        for sec_df in lst:
            sample = sec_df["__sample"].iloc[0]
            roi_s_raw = sec_df["roi"]
            if roi_s_raw.dtype == object:
                roi_s = roi_s_raw.astype(str).str.lower().isin({"true","1","t"})
            else:
                roi_s = roi_s_raw.astype(bool)
            if roi_s.sum() < 3 or (~roi_s).sum() < 3:
                continue
            sec_p = []
            sec_idx = []
            sec_rows_buf = []
            for c in metric_cols:
                v = pd.to_numeric(sec_df[c], errors="coerce").to_numpy()
                mask = ~np.isnan(v)
                v_in  = v[mask & roi_s.to_numpy()]
                v_out = v[mask & (~roi_s.to_numpy())]
                if len(v_in) < 3 or len(v_out) < 3:
                    continue
                try:
                    u, p = mannwhitneyu(v_in, v_out, alternative="two-sided")
                except ValueError:
                    continue
                sec_p.append(p)
                sec_idx.append(c)
                sec_rows_buf.append({
                    "cohort": cohort,
                    "sample": sample,
                    "metric": c,
                    "type": classify(c),
                    "mean_roi": float(np.mean(v_in)),
                    "mean_nonroi": float(np.mean(v_out)),
                    "delta": float(np.mean(v_in) - np.mean(v_out)),
                    "U_stat": float(u),
                    "p_raw": float(p),
                    "n_roi": int(len(v_in)),
                    "n_nonroi": int(len(v_out)),
                })
            if not sec_rows_buf:
                continue
            p_fdr_sec = bh_fdr(np.asarray(sec_p))
            for row, q in zip(sec_rows_buf, p_fdr_sec):
                row["p_fdr"] = float(q)
                row["sig"]   = stars(q)
                per_section_rows.append(row)

    # FDR-adjust within each cohort
    agg_df = pd.DataFrame(aggregate_rows)
    if not agg_df.empty:
        adj_list = []
        for cohort, g in agg_df.groupby("cohort", sort=False):
            p_fdr = bh_fdr(g["p_raw"].to_numpy())
            tmp = g.copy()
            tmp["p_fdr"] = p_fdr
            tmp["sig"]    = [stars(p) for p in p_fdr]
            adj_list.append(tmp)
        agg_df = pd.concat(adj_list, axis=0, ignore_index=True)
        agg_df = agg_df.sort_values(["cohort", "p_fdr", "delta"],
                                    ascending=[True, True, False])
        agg_path = OUT / "roi_vs_nonroi_aggregate_pvalues.csv"
        agg_df.to_csv(agg_path, index=False)
        print(f"\nwrote {agg_path}  ({len(agg_df)} rows)")
        # Quick summary
        for cohort, g in agg_df.groupby("cohort", sort=False):
            n_sig = (g["p_fdr"] < 0.05).sum()
            print(f"  {cohort}: {n_sig}/{len(g)} metrics FDR<0.05")

    sec_df_out = pd.DataFrame(per_section_rows)
    if not sec_df_out.empty:
        sec_path = OUT / "roi_vs_nonroi_stats_with_pvalues.csv"
        sec_df_out.to_csv(sec_path, index=False)
        print(f"wrote {sec_path}  ({len(sec_df_out)} rows)")


if __name__ == "__main__":
    main()
