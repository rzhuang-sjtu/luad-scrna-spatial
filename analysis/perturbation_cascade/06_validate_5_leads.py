"""Validate 5 candidate leads (SEC61G/SLC2A1/EGLN3/SRSF9/ANGPTL4) across:
  1. ST E-MTAB-13530 — ROI vs non-ROI spot expression
  2. ST Okamura 2024 (Takano) — ROI vs non-ROI spot expression
  3. GSE207422 / GSE126044 / GSE135222 — R vs NR sample expression

Output: per-gene per-dataset effect size + p, plus combined wide ranking table.
No plots.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
warnings.filterwarnings("ignore")

GENES = ["SEC61G", "SLC2A1", "EGLN3", "SRSF9", "ANGPTL4"]
ENSG = {"ANGPTL4": "ENSG00000167772", "EGLN3": "ENSG00000129521",
        "SEC61G": "ENSG00000132432", "SLC2A1": "ENSG00000117394",
        "SRSF9":  "ENSG00000111786"}

OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn_500_diff")
OUT.mkdir(parents=True, exist_ok=True)

ST_DATASETS = {
    "E-MTAB-13530": "${DATA_ROOT}/ST/results/step08_roi/cohort_with_roi.h5ad",
    "Okamura":      "${DATA_ROOT}/ST/results/step09_okamura_validation/cohort_with_roi.h5ad",
}

def st_scan(label, path):
    print(f"\n[ST] {label}: {path}")
    ad = sc.read_h5ad(path)
    print(f"  shape={ad.shape}, ROI=True {(ad.obs['roi']==True).sum()}, "
          f"non-ROI {(ad.obs['roi']==False).sum()}, samples={ad.obs['sample'].nunique()}")
    rows = []
    for g in GENES:
        if g not in ad.var_names:
            rows.append({"gene": g, "dataset": label, "in_data": False}); continue
        e = ad[:, g].X
        if hasattr(e, "toarray"): e = e.toarray().flatten()
        else: e = np.asarray(e).flatten()
        roi = ad.obs["roi"].values.astype(bool)
        a, b = e[roi], e[~roi]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"gene": g, "dataset": label, "in_data": True,
                     "n_roi": int(roi.sum()), "n_non": int((~roi).sum()),
                     "mean_roi": float(a.mean()), "mean_non": float(b.mean()),
                     "delta_roi_minus_non": float(a.mean() - b.mean()),
                     "mw_p": float(p)})
    return pd.DataFrame(rows)

st_dfs = [st_scan(k, v) for k, v in ST_DATASETS.items()]
st_df = pd.concat(st_dfs, ignore_index=True)
st_df.to_csv(OUT / "validate5_st.csv", index=False)
print("\n=== ST results ===")
print(st_df[["gene","dataset","mean_roi","mean_non","delta_roi_minus_non","mw_p"]].to_string(index=False))

PD_ = Path("${WORK_ROOT}/luad_figures/fig_treatment")
COHORTS = {
    "GSE207422": dict(expr="${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
                      kind="log2TPM_symbol",
                      scores=PD_ / "gse207422_mp_scores.csv",
                      pos="MPR", neg="NMPR"),
    "GSE126044": dict(expr="${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz",
                      kind="counts_symbol",
                      scores=PD_ / "gse126044_mp_scores.csv",
                      pos="R", neg="NR"),
    "GSE135222": dict(expr="${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
                      kind="TPM_ENSG",
                      scores=PD_ / "gse135222_mp_scores.csv",
                      pos="R", neg="NR"),
}

def load_expr(path, kind):
    if kind == "log2TPM_symbol":
        return pd.read_csv(path, sep="\t", index_col=0)
    elif kind == "counts_symbol":
        df = pd.read_csv(path, sep="\t", index_col=0)
        if not df.index.is_unique: df = df.groupby(level=0).max()
        cpm = df.div(df.sum(axis=0), axis=1) * 1e6
        return np.log2(cpm + 1)
    elif kind == "TPM_ENSG":
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index = df.index.astype(str).str.split(".").str[0]
        if not df.index.is_unique: df = df.groupby(level=0).max()
        return np.log2(df + 1)

def cohort_scan(label, cfg):
    print(f"\n[Tx] {label}: {cfg['pos']} vs {cfg['neg']}")
    expr = load_expr(cfg["expr"], cfg["kind"])
    use_ensg = cfg["kind"] == "TPM_ENSG"
    rows = []
    scores = pd.read_csv(cfg["scores"]).set_index("Sample")
    samples = [s for s in expr.columns if s in scores.index]
    expr = expr[samples]
    resp = scores.loc[samples, "response_group"]
    keep = resp.isin([cfg["pos"], cfg["neg"]])
    expr = expr.loc[:, keep.values]; resp = resp[keep]
    print(f"  cohort n={len(resp)}: {cfg['pos']}={(resp==cfg['pos']).sum()}, "
          f"{cfg['neg']}={(resp==cfg['neg']).sum()}")
    for g in GENES:
        idx = ENSG[g] if use_ensg else g
        if idx not in expr.index:
            rows.append({"gene": g, "cohort": label, "in_data": False}); continue
        v = expr.loc[idx].astype(float).values
        a = v[resp.values == cfg["pos"]]
        b = v[resp.values == cfg["neg"]]
        if len(a) < 3 or len(b) < 3:
            rows.append({"gene": g, "cohort": label, "in_data": True,
                         "n_pos": len(a), "n_neg": len(b),
                         "mean_pos": float(a.mean()) if len(a) else np.nan,
                         "mean_neg": float(b.mean()) if len(b) else np.nan,
                         "log2FC_pos_minus_neg": np.nan, "mw_p": np.nan,
                         "auc_pos_vs_neg": np.nan}); continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        auc = u / (len(a) * len(b))
        rows.append({"gene": g, "cohort": label, "in_data": True,
                     "pos_label": cfg["pos"], "neg_label": cfg["neg"],
                     "n_pos": len(a), "n_neg": len(b),
                     "mean_pos": float(a.mean()), "mean_neg": float(b.mean()),
                     "log2FC_pos_minus_neg": float(a.mean() - b.mean()),
                     "mw_p": float(p), "auc_pos_vs_neg": float(auc)})
    return pd.DataFrame(rows)

tx_dfs = [cohort_scan(k, v) for k, v in COHORTS.items()]
tx_df = pd.concat(tx_dfs, ignore_index=True)
tx_df.to_csv(OUT / "validate5_treatment.csv", index=False)
print("\n=== Treatment cohorts ===")
print(tx_df[["gene","cohort","pos_label","mean_pos","mean_neg","log2FC_pos_minus_neg","mw_p","auc_pos_vs_neg"]].to_string(index=False))

print("\n" + "=" * 100)
print(" Combined per-gene summary ( = sig & favorable; favorable=ROI/NR > non/R)")
print("=" * 100)
wide = pd.DataFrame({"gene": GENES})
for label, _ in ST_DATASETS.items():
    sub = st_df[st_df["dataset"] == label].set_index("gene")
    wide[f"ST_{label}_Δ"]   = wide["gene"].map(sub["delta_roi_minus_non"])
    wide[f"ST_{label}_p"]   = wide["gene"].map(sub["mw_p"])
for label in COHORTS.keys():
    sub = tx_df[tx_df["cohort"] == label].set_index("gene")
    wide[f"Tx_{label}_log2FC"] = wide["gene"].map(sub["log2FC_pos_minus_neg"])
    wide[f"Tx_{label}_p"]      = wide["gene"].map(sub["mw_p"])

# composite count of "supportive" signals.
# ST: ROI > non-ROI direction (delta > 0) AND p < 0.05 → +1 each dataset (NR-surrogate enriched if direction expected to be NR)
# Tx: log2FC < 0 AND p < 0.05 means R < NR (favors candidate as NR-marker / resistance driver) → +1 each cohort
def support_score(r):
    s = 0
    for k in ST_DATASETS.keys():
        if pd.notna(r[f"ST_{k}_Δ"]) and r[f"ST_{k}_Δ"] > 0 and r[f"ST_{k}_p"] < 0.05: s += 1
    for k in COHORTS.keys():
        if pd.notna(r[f"Tx_{k}_log2FC"]) and r[f"Tx_{k}_log2FC"] < 0 and r[f"Tx_{k}_p"] < 0.05: s += 1
    return s
wide["support_score"] = wide.apply(support_score, axis=1)
wide = wide.sort_values("support_score", ascending=False)

# pretty-print
print(f"\n  {'gene':<10} {'ST_EMTAB_Δ':>12} {'ST_EMTAB_p':>12} {'ST_Okamura_Δ':>14} {'ST_Okamura_p':>14} "
      f"{'Tx_C1_FC':>10} {'Tx_C1_p':>10} {'Tx_C2_FC':>10} {'Tx_C2_p':>10} {'Tx_C3_FC':>10} {'Tx_C3_p':>10} {'score':>6}")
print("-" * 144)
def f(x, fmt=".3f"):
    if pd.isna(x): return "    —"
    return format(x, fmt)
for _, r in wide.iterrows():
    print(f"  {r['gene']:<10} "
          f"{f(r['ST_E-MTAB-13530_Δ']):>12} {f(r['ST_E-MTAB-13530_p'],'.1e'):>12} "
          f"{f(r['ST_Okamura_Δ']):>14} {f(r['ST_Okamura_p'],'.1e'):>14} "
          f"{f(r['Tx_GSE207422_log2FC']):>10} {f(r['Tx_GSE207422_p'],'.1e'):>10} "
          f"{f(r['Tx_GSE126044_log2FC']):>10} {f(r['Tx_GSE126044_p'],'.1e'):>10} "
          f"{f(r['Tx_GSE135222_log2FC']):>10} {f(r['Tx_GSE135222_p'],'.1e'):>10} "
          f"{int(r['support_score']):>6}")

wide.to_csv(OUT / "validate5_combined.csv", index=False)
print(f"\nWrote: {OUT/'validate5_st.csv'}, {OUT/'validate5_treatment.csv'}, {OUT/'validate5_combined.csv'}")
print("DONE.")
