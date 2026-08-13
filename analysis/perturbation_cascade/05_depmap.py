"""Fig 8B/C/G — DepMap 24Q2 validation, restyled to project standard."""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import sys
from pathlib import Path
# fig8_style lives with the figure it styles; make it importable when this
# script is run from its own directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "plotting" / "fig8"))
from fig8_style import setup, COL, save_panel, sig_stars, add_subtitle, style_axes
setup()

D = Path("${DATA_ROOT}/depmap/24Q2")
CAND_FILE = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn/8A_candidate_pool.csv")
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8B_depmap")
OUT.mkdir(parents=True, exist_ok=True)

cand = pd.read_csv(CAND_FILE)
CANDS = list(cand["Gene_name"])
print(f"candidates: {CANDS}")

model = pd.read_csv(D / "Model.csv", low_memory=False)
luad_models = set(model.loc[model["OncotreeCode"] == "LUAD", "ModelID"])
print(f"LUAD lines: {len(luad_models)}")

crispr_header = pd.read_csv(D / "CRISPRGeneEffect.csv", nrows=0).columns.tolist()
def name_of(c):
    m = re.match(r"^([A-Za-z0-9\-]+)\s*\(\d+\)$", c)
    return m.group(1) if m else c
sym2col = {name_of(c): c for c in crispr_header[1:]}
sel_cols = [crispr_header[0]] + [sym2col[g] for g in CANDS if g in sym2col]
crispr = pd.read_csv(D / "CRISPRGeneEffect.csv", usecols=sel_cols)
crispr = crispr.rename(columns={crispr_header[0]: "ModelID"})
crispr = crispr.rename(columns={sym2col[g]: g for g in CANDS if g in sym2col})
crispr["is_LUAD"] = crispr["ModelID"].isin(luad_models)

print("\n=== 8B ===")
long = crispr.melt(id_vars=["ModelID", "is_LUAD"],
                   value_vars=[g for g in CANDS if g in sym2col],
                   var_name="gene", value_name="gene_effect")
long["group"] = np.where(long["is_LUAD"], "LUAD", "non-LUAD")
long.to_csv(OUT / "8B_data.csv", index=False)

stat_rows = []
for g in CANDS:
    if g not in sym2col: continue
    a = long.loc[(long["gene"] == g) & long["is_LUAD"], "gene_effect"].dropna()
    b = long.loc[(long["gene"] == g) & ~long["is_LUAD"], "gene_effect"].dropna()
    u, p = (stats.mannwhitneyu(a, b, alternative="less") if len(a) >= 5 and len(b) >= 5
            else (np.nan, np.nan))
    stat_rows.append({"gene": g, "luad_mean": a.mean(), "luad_n": len(a),
                      "other_mean": b.mean(), "other_n": len(b),
                      "delta_luad_minus_other": a.mean() - b.mean(),
                      "mannwhitney_p_LUAD<other": p})
stat_df = pd.DataFrame(stat_rows).sort_values("luad_mean")
stat_df.to_csv(OUT / "8B_stats.csv", index=False)
print(stat_df.to_string(index=False))

order = stat_df["gene"].tolist()
fig, ax = plt.subplots(figsize=(5.0, 2.8))
sns.violinplot(data=long, x="gene", y="gene_effect", hue="group",
               order=order, hue_order=["LUAD", "non-LUAD"],
               split=True, inner="quartile", cut=0, linewidth=0.5,
               palette={"LUAD": COL["LUAD"], "non-LUAD": COL["other"]}, ax=ax)
ax.axhline(0, color="black", lw=0.4, ls=":")
ax.axhline(-0.5, color=COL["ref_red"], lw=0.4, ls="--", alpha=0.6,
           label="essentiality cutoff")
for i, g in enumerate(order):
    p = stat_df.loc[stat_df["gene"] == g, "mannwhitney_p_LUAD<other"].iloc[0]
    ax.text(i, ax.get_ylim()[1] * 0.96, sig_stars(p),
            ha="center", fontsize=6)
ax.set_xlabel("")
style_axes(ax, ylabel="CRISPR Gene Effect")
ax.legend(loc="lower right", fontsize=6)
add_subtitle(ax, f"LUAD n={int(stat_df['luad_n'].iloc[0])} vs non-LUAD n={int(stat_df['other_n'].iloc[0])} · DepMap 24Q2")
fig.tight_layout()
save_panel(fig, OUT / "8B_violin_LUAD_vs_other")
plt.close(fig)

print("\n=== 8C ===")
luad_only = long[long["is_LUAD"]].copy()
agg = luad_only.groupby("gene")["gene_effect"].agg(["mean", "std", "count"]).reset_index()
agg["sem"] = agg["std"] / np.sqrt(agg["count"])
agg = agg.sort_values("mean")
agg.to_csv(OUT / "8C_data.csv", index=False)
print(agg.to_string(index=False))

colors = [COL["tumor"] if m < -0.5 else (COL["Macro_SPP1"] if m < 0 else COL["other"])
          for m in agg["mean"]]
fig, ax = plt.subplots(figsize=(4.0, 2.6))
ax.bar(agg["gene"], agg["mean"], yerr=agg["sem"], capsize=2,
       color=colors, edgecolor="black", linewidth=0.4, error_kw={"lw": 0.5})
ax.axhline(0, color="black", lw=0.4)
ax.axhline(-0.5, color=COL["ref_red"], lw=0.4, ls="--", alpha=0.6,
           label="essentiality cutoff")
ax.legend(loc="lower right", fontsize=6)
style_axes(ax, ylabel="Mean CRISPR Gene Effect")
ax.tick_params(axis="x", rotation=0)
add_subtitle(ax, f"LUAD lines, n={int(agg['count'].iloc[0])} · DepMap 24Q2")
fig.tight_layout()
save_panel(fig, OUT / "8C_bar_luad_essentiality")
plt.close(fig)

print("\n=== 8G ===")
expr_header = pd.read_csv(D / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                          nrows=0).columns.tolist()
meta_cols = ["SequencingID", "ModelConditionID", "ModelID",
             "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
expr_id = expr_header[0]
expr_sym2col = {name_of(c): c for c in expr_header
                if c not in meta_cols and c != expr_id}
use = [expr_id, "ModelID", "IsDefaultEntryForModel"] + \
      [expr_sym2col[g] for g in CANDS if g in expr_sym2col]
expr = pd.read_csv(D / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", usecols=use)
expr = expr[expr["IsDefaultEntryForModel"] == "Yes"].copy()
expr = expr.rename(columns={expr_sym2col[g]: g for g in CANDS if g in expr_sym2col})

luad_crispr = crispr[crispr["is_LUAD"]].set_index("ModelID")
luad_expr = expr[expr["ModelID"].isin(luad_models)].set_index("ModelID")
shared = luad_crispr.index.intersection(luad_expr.index)
print(f"  LUAD lines with both: {len(shared)}")

rows = []
for g in CANDS:
    if g not in sym2col or g not in expr_sym2col: continue
    for m in shared:
        rows.append({"gene": g, "ModelID": m,
                     "log2_TPM_p1": luad_expr.loc[m, g],
                     "gene_effect": luad_crispr.loc[m, g]})
sc = pd.DataFrame(rows).dropna()
sc.to_csv(OUT / "8G_data.csv", index=False)

genes = [g for g in CANDS if g in sym2col and g in expr_sym2col]
n = len(genes); ncols = 4; nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(2.4 * ncols, 2.0 * nrows), squeeze=False)
for i, g in enumerate(genes):
    ax = axes[i // ncols, i % ncols]
    sub = sc[sc["gene"] == g]
    ax.scatter(sub["log2_TPM_p1"], sub["gene_effect"],
               s=10, alpha=0.7, color=COL["normal"],
               edgecolor="white", lw=0.3)
    if len(sub) >= 5:
        r, p = stats.spearmanr(sub["log2_TPM_p1"], sub["gene_effect"])
        ax.text(0.04, 0.92, f"ρ={r:.2f}\np={p:.1e}",
                transform=ax.transAxes, fontsize=6, va="top")
    ax.axhline(0, color="black", lw=0.3, ls=":")
    ax.axhline(-0.5, color=COL["ref_red"], lw=0.3, ls="--", alpha=0.5)
    style_axes(ax, xlabel="log2(TPM+1)", ylabel="Gene Effect")
    add_subtitle(ax, g)
for j in range(n, nrows * ncols):
    axes[j // ncols, j % ncols].axis("off")
fig.tight_layout()
save_panel(fig, OUT / "8G_scatter_expr_vs_effect")
plt.close(fig)

print("\nDONE.")
