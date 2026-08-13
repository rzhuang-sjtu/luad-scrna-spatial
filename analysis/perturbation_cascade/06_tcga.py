"""Fig 8D/E/F (top-3 candidates × T vs N), 8G (1×3 expr×effect), 8O (HSP90AB1 KM),
8P (SRSF9 KM) — restyled.

Per checklist: 8D=HSP90AB1 T vs N, 8E=SRSF9 T vs N, 8F=NDUFB2 T vs N (each its own panel).
8G is 1×3 grid for the same top-3.
8O = HSP90AB1 KM only (large), 8P = SRSF9 KM only (large).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import re

import sys
from pathlib import Path
# fig8_style lives with the figure it styles; make it importable when this
# script is run from its own directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "plotting" / "fig8"))
from fig8_style import setup, COL, save_panel, sig_stars, add_subtitle, style_axes
setup()

T = Path("${DATA_ROOT}/TCGA_LUAD_analysis")
DEPMAP = Path("${DATA_ROOT}/depmap/24Q2")
CAND_FILE = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn/8A_candidate_pool.csv")
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8D_tcga")
OUT.mkdir(parents=True, exist_ok=True)

cand = pd.read_csv(CAND_FILE)
ALL_CANDS = list(cand["Gene_name"])
TOP3 = ["HSP90AB1", "SRSF9", "NDUFB2"]

print("\nloading clinical+TPM...")
cln = pd.read_csv(T / "TCGA_LUAD_clinical.csv").rename(
    columns={"sample_barcode": "sample"})
tpm_rows = []
for chunk in pd.read_csv(T / "TCGA_LUAD_TPM_matrix.csv", chunksize=5000):
    chunk = chunk.rename(columns={chunk.columns[0]: "gene"})
    sub = chunk[chunk["gene"].isin(ALL_CANDS)]
    if len(sub): tpm_rows.append(sub)
tpm = pd.concat(tpm_rows, ignore_index=True).set_index("gene")
print(f"  TPM: {tpm.shape}")

sample_type_map = dict(zip(cln["sample"], cln["sample_type"]))
long = []
for g in tpm.index:
    for s in tpm.columns:
        st = sample_type_map.get(s, "Unknown")
        if st in ("Primary Tumor", "Solid Tissue Normal"):
            long.append({"gene": g, "sample": s,
                         "log2_TPM_p1": np.log2(tpm.loc[g, s] + 1),
                         "type": "Tumor" if st == "Primary Tumor" else "Normal"})
long = pd.DataFrame(long)
long.to_csv(OUT / "8D_data.csv", index=False)

stat_rows = []
for g in ALL_CANDS:
    sub = long[long["gene"] == g]
    t = sub.loc[sub["type"] == "Tumor", "log2_TPM_p1"]
    n = sub.loc[sub["type"] == "Normal", "log2_TPM_p1"]
    u, p = stats.mannwhitneyu(t, n, alternative="two-sided")
    stat_rows.append({"gene": g, "tumor_mean": t.mean(), "tumor_n": len(t),
                      "normal_mean": n.mean(), "normal_n": len(n),
                      "log2FC_T_minus_N": t.mean() - n.mean(), "wilcoxon_p": p})
stat_df = pd.DataFrame(stat_rows)
stat_df.to_csv(OUT / "8E_volcano_summary.csv", index=False)


def plot_tvn_single(g, out_stem):
    sub = long[long["gene"] == g]
    p = stat_df.loc[stat_df["gene"] == g, "wilcoxon_p"].iloc[0]
    fc = stat_df.loc[stat_df["gene"] == g, "log2FC_T_minus_N"].iloc[0]
    fig, ax = plt.subplots(figsize=(2.0, 2.4))
    sns.boxplot(data=sub, x="type", y="log2_TPM_p1", order=["Normal", "Tumor"],
                hue="type", hue_order=["Normal", "Tumor"], legend=False,
                palette={"Normal": COL["normal"], "Tumor": COL["tumor"]},
                showfliers=False, width=0.55, linewidth=0.5, ax=ax)
    sns.stripplot(data=sub, x="type", y="log2_TPM_p1", order=["Normal", "Tumor"],
                  color="black", size=0.8, alpha=0.35, ax=ax)
    add_subtitle(ax, f"{g}  log2FC={fc:+.2f}  {sig_stars(p)} p={p:.1e}")
    style_axes(ax, ylabel="log2(TPM+1)", xlabel="")
    fig.tight_layout()
    save_panel(fig, out_stem)
    plt.close(fig)


print("\n=== 8D/E/F (single-gene T vs N panels) ===")
for letter, g in zip("DEF", TOP3):
    plot_tvn_single(g, OUT / f"8{letter}_tvn_{g}")
    print(f"  wrote 8{letter}_tvn_{g}")

print("\n=== 8G (1×3 top-3) ===")
crispr_header = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", nrows=0).columns.tolist()
def name_of(c):
    m = re.match(r"^([A-Za-z0-9\-]+)\s*\(\d+\)$", c)
    return m.group(1) if m else c
sym2col = {name_of(c): c for c in crispr_header[1:]}
sel = [crispr_header[0]] + [sym2col[g] for g in TOP3 if g in sym2col]
crispr = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", usecols=sel).rename(
    columns={crispr_header[0]: "ModelID"})
crispr = crispr.rename(columns={sym2col[g]: g for g in TOP3 if g in sym2col})

model = pd.read_csv(DEPMAP / "Model.csv", low_memory=False)
luad_models = set(model.loc[model["OncotreeCode"] == "LUAD", "ModelID"])

expr_header = pd.read_csv(DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                          nrows=0).columns.tolist()
meta_cols = ["SequencingID", "ModelConditionID", "ModelID",
             "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
expr_id = expr_header[0]
expr_sym2col = {name_of(c): c for c in expr_header
                if c not in meta_cols and c != expr_id}
use = [expr_id, "ModelID", "IsDefaultEntryForModel"] + \
      [expr_sym2col[g] for g in TOP3 if g in expr_sym2col]
expr = pd.read_csv(DEPMAP / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", usecols=use)
expr = expr[expr["IsDefaultEntryForModel"] == "Yes"].rename(
    columns={expr_sym2col[g]: g for g in TOP3 if g in expr_sym2col})

luad_crispr = crispr[crispr["ModelID"].isin(luad_models)].set_index("ModelID")
luad_expr = expr[expr["ModelID"].isin(luad_models)].set_index("ModelID")
shared = luad_crispr.index.intersection(luad_expr.index)
print(f"  LUAD lines with both: {len(shared)}")

rows = []
for g in TOP3:
    if g not in sym2col or g not in expr_sym2col: continue
    for m in shared:
        rows.append({"gene": g, "ModelID": m,
                     "log2_TPM_p1": luad_expr.loc[m, g],
                     "gene_effect": luad_crispr.loc[m, g]})
sc = pd.DataFrame(rows).dropna()
sc.to_csv(OUT / "8G_data.csv", index=False)

fig, axes = plt.subplots(1, len(TOP3), figsize=(2.0 * len(TOP3), 2.0), squeeze=False)
for i, g in enumerate(TOP3):
    ax = axes[0, i]
    sub = sc[sc["gene"] == g]
    ax.scatter(sub["log2_TPM_p1"], sub["gene_effect"],
               s=12, alpha=0.75, color=COL["normal"], edgecolor="white", lw=0.3)
    if len(sub) >= 5:
        r, pp = stats.spearmanr(sub["log2_TPM_p1"], sub["gene_effect"])
        ax.text(0.04, 0.92, f"ρ={r:.2f}\np={pp:.1e}",
                transform=ax.transAxes, fontsize=6, va="top")
    ax.axhline(0, color="black", lw=0.3, ls=":")
    ax.axhline(-0.5, color=COL["ref_red"], lw=0.3, ls="--", alpha=0.5)
    style_axes(ax, xlabel="log2(TPM+1)", ylabel="CRISPR Gene Effect" if i == 0 else "")
    add_subtitle(ax, g)
fig.tight_layout()
save_panel(fig, OUT / "8G_scatter_top3")
plt.close(fig)

print("\n=== 8O / 8P KM ===")
tumor_samples = [s for s in cln.loc[cln["sample_type"] == "Primary Tumor", "sample"]
                 if s in tpm.columns]
surv = cln[cln["sample"].isin(tumor_samples)].copy()
surv["event"] = (surv["vital_status"] == "Dead").astype(int)
surv["time"] = np.where(surv["event"] == 1,
                        surv["days_to_death"], surv["days_to_last_follow_up"])
surv = surv.dropna(subset=["time"])
surv = surv[surv["time"] > 0].copy()


def plot_km_single(gene, out_stem):
    expr_g = tpm.loc[gene, surv["sample"].values].values
    s = surv.copy()
    s["expr"] = np.log2(expr_g + 1)
    med = s["expr"].median()
    s["group"] = np.where(s["expr"] >= med, "High", "Low")
    high = s[s["group"] == "High"]; low = s[s["group"] == "Low"]
    lr = logrank_test(high["time"], low["time"], high["event"], low["event"])
    p = lr.p_value
    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    KaplanMeierFitter().fit(high["time"], high["event"],
                             label=f"High (n={len(high)})").plot_survival_function(
        ax=ax, ci_show=False, color=COL["high"], lw=1.0)
    KaplanMeierFitter().fit(low["time"], low["event"],
                             label=f"Low (n={len(low)})").plot_survival_function(
        ax=ax, ci_show=False, color=COL["low"], lw=1.0)
    ax.set_xlim(0, s["time"].max())
    style_axes(ax, xlabel="Days", ylabel="Overall survival")
    ax.legend(fontsize=6, loc="upper right", frameon=False)
    add_subtitle(ax, f"{gene}  logrank p={p:.1e} {sig_stars(p)} · TCGA-LUAD n={len(s)} ({int(s['event'].sum())} events)")
    fig.tight_layout()
    save_panel(fig, out_stem)
    plt.close(fig)
    return p, len(high), len(low)


km_rows = []
for letter, g in zip(["O", "P"], ["HSP90AB1", "SRSF9"]):
    p, nh, nl = plot_km_single(g, OUT / f"8{letter}_km_{g}")
    km_rows.append({"gene": g, "n_high": nh, "n_low": nl, "logrank_p": p})
    print(f"  8{letter} {g}: p={p:.3g}")
pd.DataFrame(km_rows).to_csv(OUT / "8O_8P_km_data.csv", index=False)

# Cox forest (kept as supplementary data, not a main panel anymore)
cox_rows = []
for g in ALL_CANDS:
    if g not in tpm.index: continue
    expr_g = tpm.loc[g, surv["sample"].values].values
    df = pd.DataFrame({"time": surv["time"].values,
                       "event": surv["event"].values,
                       "expr": np.log2(expr_g + 1)})
    cph = CoxPHFitter()
    cph.fit(df, duration_col="time", event_col="event")
    s = cph.summary.loc["expr"]
    cox_rows.append({"gene": g, "HR": np.exp(s["coef"]),
                     "HR_low": np.exp(s["coef lower 95%"]),
                     "HR_high": np.exp(s["coef upper 95%"]),
                     "p": s["p"], "n": len(df), "events": df["event"].sum()})
pd.DataFrame(cox_rows).to_csv(OUT / "cox_forest_data.csv", index=False)

print("\nDONE.")
