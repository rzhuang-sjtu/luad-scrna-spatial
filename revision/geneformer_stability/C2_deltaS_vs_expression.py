"""C2 — ΔS vs expression / detection frequency, and expression-matched null model. No GPU.

Motivation:
  "correlation of delta_S to goal with mean expression and detection frequency,
   with an expression-matched null"

Core question
----------
If Shift_to_goal_end (ΔS) is driven mainly by expression level or detection frequency, then «deleting this gene pushes
cells toward the target state» is not a biological signal but the mechanical effect of «highly expressed genes cause large rank-encoding
perturbations when removed». Quantify this as requested.

Approach
----
1. For each transition, take the published 500-cell result table and compute Spearman of ΔS with
   (a) mean expression in sender cells  (b) detection frequency = N_Detections / n_cells
   .
2. Expression-matched null: bin genes by expression and assess top-200 enrichment within bins.
   If ΔS is fully expression-driven, top-200 concentrates in the highest-expression bins.
3. Compare expression/detection distributions of top-200 vs all genes (Mann-Whitney).

Write CSV under results/fig8_geneformer/deltaS_expression/
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

ROOT = "${PROJECT_ROOT}/results/fig8_geneformer"
OUT = f"{ROOT}/deltaS_expression"
FDR, N_MIN, TOPN = 0.05, 5, 200

# Published-run inputs (neu used the HVG-truncated version — the bug identified here)
INPUTS = {
    "macro_spp1_to_c1qc":     (f"{ROOT}/inputs/macro_spp1_to_c1qc/data.h5ad", "v1"),
    "mal_mp3_to_mp1":         (f"{ROOT}/inputs/mal_mp3_to_mp1/data.h5ad", "v1"),
    "neu_osm_priming_to_low": (f"{ROOT}/inputs/neu_osm_priming_to_low/data.h5ad", "v1"),
}


def sender_expression(h5ad):
    """Per-gene mean log1p expression and detection frequency in sender cells (indexed by ENSEMBL id)."""
    import anndata as ad
    a = ad.read_h5ad(h5ad)
    m = (a.obs["cell_state"].astype(str) == "sender").values
    X = a.X[m]
    X = sp.csr_matrix(X) if not sp.issparse(X) else X
    n = X.shape[0]
    mean_cnt = np.asarray(X.mean(axis=0)).ravel()
    det = np.asarray((X > 0).sum(axis=0)).ravel() / n
    return pd.DataFrame({"Ensembl_ID": list(a.var_names),
                         "mean_count": mean_cnt,
                         "log_mean": np.log1p(mean_cnt),
                         "detect_frac": det,
                         "n_sender": n})


def load_stats(d, t):
    f = os.path.join(d, t, f"{t}_stats.csv")
    if not os.path.exists(f):
        return None
    s = pd.read_csv(f, index_col=0)
    s["passed"] = ((s.Sig == 1) & (s.Shift_to_goal_end > 0) &
                   (s.Goal_end_FDR < FDR) & (s.N_Detections >= N_MIN))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perturb-dir", default=f"{ROOT}/perturb_500")
    ap.add_argument("--tag", default="published500")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows_corr, rows_strat, merged_all = [], [], []

    for t, (h5, _) in INPUTS.items():
        s = load_stats(a.perturb_dir, t)
        if s is None:
            print(f"Skip {t} (no result table)"); continue
        e = sender_expression(h5)
        m = s.merge(e, on="Ensembl_ID", how="inner")
        m["transition"] = t
        m["is_top200"] = False
        pas = m[m.passed].sort_values("Shift_to_goal_end", ascending=False)
        m.loc[m.index.isin(pas.head(TOPN).index), "is_top200"] = True
        merged_all.append(m)
        print(f"\n=== {t} ===  genes {len(m)}, passed filter {int(m.passed.sum())},"
              f"sender cells {int(m.n_sender.iloc[0])}")

        # (1) Correlation of ΔS with expression / detection frequency
        for sub, lab in [(m, "All genes"), (m[m.passed], "genes that passed the filter")]:
            if len(sub) < 30: continue
            for xcol, xlab in [("log_mean", "log mean expression"), ("detect_frac", "detection frequency"),
                               ("N_Detections", "N_Detections")]:
                rho, p = stats.spearmanr(sub.Shift_to_goal_end, sub[xcol])
                rows_corr.append({"transition": t, "gene_set": lab, "x": xlab,
                                  "spearman_rho": rho, "p": p, "n": len(sub)})
                print(f"  [{lab:14s}] ΔS vs {xlab:12s}  rho={rho:+.3f}  p={p:.2g}  n={len(sub)}")

        # (2) Expression-matched null: 10 expression strata; top-200 distribution across strata
        mm = m[m.passed].copy()
        if len(mm) >= 50:
            mm["expr_decile"] = pd.qcut(mm.log_mean, 10, labels=False, duplicates="drop")
            tot = mm.groupby("expr_decile").size()
            hit = mm[mm.is_top200].groupby("expr_decile").size().reindex(tot.index, fill_value=0)
            for d in tot.index:
                rows_strat.append({"transition": t, "expr_decile": int(d),
                                   "n_genes": int(tot[d]), "n_top200": int(hit[d]),
                                   "frac_top200": float(hit[d] / tot[d])})
            print("Fraction of top-200 in each expression decile:" +
                  " ".join(f"{d}:{hit[d]/tot[d]:.2f}" for d in tot.index))

        # (3) Expression difference: top-200 vs other genes that passed the filter
        if m.is_top200.sum() > 10:
            g1 = m.loc[m.is_top200, "log_mean"]
            g2 = m.loc[m.passed & ~m.is_top200, "log_mean"]
            if len(g2) > 10:
                u, p = stats.mannwhitneyu(g1, g2)
                print(f"top-200 log mean expression {g1.mean():.3f} vs other passing genes {g2.mean():.3f}"
                      f"  MW p={p:.2g}")
                rows_corr.append({"transition": t, "gene_set": "top200 vs rest (expression)",
                                  "x": "log_mean", "spearman_rho": float(g1.mean() - g2.mean()),
                                  "p": p, "n": len(g1)})

    pd.DataFrame(rows_corr).to_csv(f"{OUT}/{a.tag}_correlations.csv", index=False)
    pd.DataFrame(rows_strat).to_csv(f"{OUT}/{a.tag}_expr_decile.csv", index=False)
    if merged_all:
        pd.concat(merged_all, ignore_index=True)[
            ["transition", "Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR",
             "N_Detections", "log_mean", "detect_frac", "passed", "is_top200"]
        ].to_csv(f"{OUT}/{a.tag}_gene_level.csv.gz", index=False)
    print(f"\nResults written to {OUT}/")


if __name__ == "__main__":
    main()
