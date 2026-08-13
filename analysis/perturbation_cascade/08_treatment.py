"""Fig 8M (cohort #1 boxplot 3 genes) + 8N (cohort #2 volcano) + S11D/E/F
(cohort #3 single-gene boxplot per gene).

Layout per checklist:
  8M : GSE207422 (neoadjuvant chemo-IO, MPR vs NMPR), 3 candidate genes side-by-side boxplot.
  8N : GSE126044 (anti-PD-1, R vs NR), genome-wide volcano plot, candidates highlighted.
  S11D / E / F : GSE135222 (anti-PD-1/L1, DCB R vs NR), one panel per top-3 gene boxplot.
"""
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

PD_ = Path("${WORK_ROOT}/luad_figures/fig_treatment")
CAND_FILE = Path("${PROJECT_ROOT}/results/fig8_plot_data/8A_venn/8A_candidate_pool.csv")
OUT = Path("${PROJECT_ROOT}/results/fig8_plot_data/8N_treatment")
OUT.mkdir(parents=True, exist_ok=True)

cand = pd.read_csv(CAND_FILE)
ALL_CANDS = list(cand["Gene_name"])
ENSG_OF = dict(zip(cand["Gene_name"], cand["Ensembl_ID"]))
TOP3 = ["HSP90AB1", "SRSF9", "NDUFB2"]

C1 = dict(name="GSE207422",
          expr="${DATA_ROOT}/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
          kind="log2TPM_symbol",
          scores=PD_ / "gse207422_mp_scores.csv",
          pos="MPR", neg="NMPR")
C2 = dict(name="GSE126044",
          expr="${DATA_ROOT}/GSE126044/GSE126044_counts.txt.gz",
          kind="counts_symbol",
          scores=PD_ / "gse126044_mp_scores.csv",
          pos="R", neg="NR")
C3 = dict(name="GSE135222",
          expr="${DATA_ROOT}/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
          kind="TPM_ENSG",
          scores=PD_ / "gse135222_mp_scores.csv",
          pos="R", neg="NR")


def load_expr_full(path, kind):
    """Whole-genome expression matrix, log2 scale, gene-symbol indexed."""
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
        # keep ENSG-indexed; we'll relabel candidates only
        return np.log2(df + 1)


def cohort_long_for(cohort_cfg, gene_list, ensg_map=False):
    """Return long DataFrame {sample,gene,expr,response} for the genes in gene_list.
    If ensg_map=True, gene_list are symbols and the matrix is ENSG-indexed."""
    expr = load_expr_full(cohort_cfg["expr"], cohort_cfg["kind"])
    sub_idx = [ENSG_OF[g] for g in gene_list] if ensg_map else gene_list
    sub_idx = [i for i in sub_idx if i in expr.index]
    expr_sub = expr.loc[sub_idx].copy()
    if ensg_map:
        e2s = {v: k for k, v in ENSG_OF.items()}
        expr_sub.index = [e2s[i] for i in expr_sub.index]
    scores = pd.read_csv(cohort_cfg["scores"]).set_index("Sample")
    samples = [s for s in expr_sub.columns if s in scores.index]
    expr_sub = expr_sub[samples]
    resp = scores.loc[samples, "response_group"]
    keep = resp.isin([cohort_cfg["pos"], cohort_cfg["neg"]])
    expr_sub = expr_sub.loc[:, keep.values]; resp = resp[keep]
    rows = []
    for g in expr_sub.index:
        for s in expr_sub.columns:
            rows.append({"sample": s, "gene": g, "expr": expr_sub.loc[g, s],
                         "response": resp[s]})
    return pd.DataFrame(rows), expr_sub, resp


print("\n=== 8M (GSE207422, MPR vs NMPR, 3 genes) ===")
long_c1, expr_c1, resp_c1 = cohort_long_for(C1, TOP3, ensg_map=False)
order_c1 = [C1["neg"], C1["pos"]]
fig, axes = plt.subplots(1, 3, figsize=(2.0 * 3, 2.4), squeeze=False)
for i, g in enumerate(TOP3):
    ax = axes[0, i]
    sub = long_c1[long_c1["gene"] == g]
    sns.boxplot(data=sub, x="response", y="expr", order=order_c1,
                hue="response", hue_order=order_c1, legend=False,
                palette={C1["pos"]: COL["R"], C1["neg"]: COL["NR"]},
                showfliers=False, width=0.55, linewidth=0.5, ax=ax)
    sns.stripplot(data=sub, x="response", y="expr", order=order_c1,
                  color="black", size=2, alpha=0.7, ax=ax)
    a = sub.loc[sub["response"] == C1["pos"], "expr"]
    b = sub.loc[sub["response"] == C1["neg"], "expr"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    auc = u / (len(a) * len(b))
    add_subtitle(ax, f"{g}  AUC={auc:.2f} p={p:.2g} {sig_stars(p)}")
    style_axes(ax, ylabel="log2(TPM+1)" if i == 0 else "", xlabel="")
fig.text(0.01, 1.00, f"GSE207422 — MPR (n={(resp_c1 == C1['pos']).sum()}) vs NMPR (n={(resp_c1 == C1['neg']).sum()})",
         fontsize=7, style="italic", va="top",
         transform=fig.transFigure)
fig.tight_layout()
save_panel(fig, OUT / "8M_GSE207422_3genes")
plt.close(fig)
long_c1.to_csv(OUT / "8M_long.csv", index=False)

print("\n=== 8N (GSE126044 R vs NR, volcano) ===")
expr_full_c2 = load_expr_full(C2["expr"], C2["kind"])
scores_c2 = pd.read_csv(C2["scores"]).set_index("Sample")
samples_c2 = [s for s in expr_full_c2.columns if s in scores_c2.index]
expr_full_c2 = expr_full_c2[samples_c2]
resp_c2 = scores_c2.loc[samples_c2, "response_group"]
keep = resp_c2.isin([C2["pos"], C2["neg"]])
expr_full_c2 = expr_full_c2.loc[:, keep.values]
resp_c2 = resp_c2[keep]
print(f"  n samples: {len(resp_c2)}; pos {C2['pos']} = {(resp_c2 == C2['pos']).sum()}; "
      f"neg {C2['neg']} = {(resp_c2 == C2['neg']).sum()}")

# filter low-expressed genes
mu = expr_full_c2.mean(axis=1)
expr_full_c2 = expr_full_c2.loc[mu > 0.5]
print(f"  genes after expr filter: {len(expr_full_c2)}")

pos_mask = (resp_c2 == C2["pos"]).values
neg_mask = (resp_c2 == C2["neg"]).values
mean_pos = expr_full_c2.loc[:, pos_mask].mean(axis=1)
mean_neg = expr_full_c2.loc[:, neg_mask].mean(axis=1)
log2fc = mean_pos - mean_neg   # already log2

# vectorized Wilcoxon would be slow for 18k genes — use t-test for volcano (standard practice)
from scipy.stats import ttest_ind
A = expr_full_c2.loc[:, pos_mask].values
B = expr_full_c2.loc[:, neg_mask].values
t_stat, pvals = ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
volc = pd.DataFrame({"gene": expr_full_c2.index, "log2FC": log2fc.values,
                     "p": pvals, "neg_log10_p": -np.log10(np.maximum(pvals, 1e-300))})
volc.to_csv(OUT / "8N_volcano_data.csv", index=False)

fig, ax = plt.subplots(figsize=(4.0, 3.2))
# background
sig_mask = (volc["p"] < 0.05)
ax.scatter(volc.loc[~sig_mask, "log2FC"], volc.loc[~sig_mask, "neg_log10_p"],
           s=2, color="lightgray", alpha=0.6, rasterized=True)
ax.scatter(volc.loc[sig_mask & (volc["log2FC"] > 0), "log2FC"],
           volc.loc[sig_mask & (volc["log2FC"] > 0), "neg_log10_p"],
           s=3, color=COL["NR"], alpha=0.7, rasterized=True,
           label=f"{C2['pos']} ↑ (p<0.05)")
ax.scatter(volc.loc[sig_mask & (volc["log2FC"] < 0), "log2FC"],
           volc.loc[sig_mask & (volc["log2FC"] < 0), "neg_log10_p"],
           s=3, color=COL["R"], alpha=0.7, rasterized=True,
           label=f"{C2['neg']} ↑ (p<0.05)")
# highlight candidates
for g in ALL_CANDS:
    sub = volc[volc["gene"] == g]
    if len(sub) == 0: continue
    x, y = sub["log2FC"].iloc[0], sub["neg_log10_p"].iloc[0]
    color = COL["high"] if g in TOP3 else "black"
    size = 30 if g in TOP3 else 15
    ax.scatter(x, y, s=size, color=color, edgecolor="black",
               linewidth=0.5, zorder=5)
    ax.annotate(g, (x, y), xytext=(4, 4), textcoords="offset points",
                fontsize=6, fontweight="bold" if g in TOP3 else "normal")
ax.axhline(-np.log10(0.05), color="black", lw=0.4, ls="--")
ax.axvline(0, color="black", lw=0.4, ls="--")
style_axes(ax, xlabel=f"log2 fold change ({C2['pos']} − {C2['neg']})",
           ylabel="-log10 p")
ax.legend(fontsize=6, frameon=False, loc="upper right")
add_subtitle(ax, f"{C2['name']} · R={(resp_c2 == C2['pos']).sum()} / NR={(resp_c2 == C2['neg']).sum()} · candidates highlighted")
fig.tight_layout()
save_panel(fig, OUT / "8N_GSE126044_volcano")
plt.close(fig)

print("\n=== S11D/E/F (GSE135222 R vs NR per gene) ===")
long_c3, expr_c3, resp_c3 = cohort_long_for(C3, TOP3, ensg_map=True)
order_c3 = [C3["neg"], C3["pos"]]
for letter, g in zip(["D", "E", "F"], TOP3):
    sub = long_c3[long_c3["gene"] == g]
    if len(sub) == 0: continue
    fig, ax = plt.subplots(figsize=(2.0, 2.4))
    sns.boxplot(data=sub, x="response", y="expr", order=order_c3,
                hue="response", hue_order=order_c3, legend=False,
                palette={C3["pos"]: COL["R"], C3["neg"]: COL["NR"]},
                showfliers=False, width=0.55, linewidth=0.5, ax=ax)
    sns.stripplot(data=sub, x="response", y="expr", order=order_c3,
                  color="black", size=2, alpha=0.7, ax=ax)
    a = sub.loc[sub["response"] == C3["pos"], "expr"]
    b = sub.loc[sub["response"] == C3["neg"], "expr"]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    auc = u / (len(a) * len(b))
    add_subtitle(ax, f"{g}  AUC={auc:.2f} p={p:.2g} {sig_stars(p)}")
    style_axes(ax, ylabel="log2(TPM+1)", xlabel="")
    fig.text(0.01, 0.98, f"GSE135222 · DCB {C3['pos']} (n={(resp_c3 == C3['pos']).sum()}) vs {C3['neg']} (n={(resp_c3 == C3['neg']).sum()})",
             fontsize=7, style="italic", va="top")
    fig.tight_layout()
    save_panel(fig, OUT / f"S11{letter}_GSE135222_{g}")
    plt.close(fig)
long_c3.to_csv(OUT / "S11_long.csv", index=False)

print("\nDONE.")
