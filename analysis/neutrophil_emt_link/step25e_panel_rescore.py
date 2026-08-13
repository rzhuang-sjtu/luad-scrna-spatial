"""Re-score literature TAN/NAN panels on own_raw (9698 genes) instead of HVG-3000.

Then refresh the panel heatmap + summary stub. Also flags markers MISSING from upstream
HVG-subset, so the user knows what's not available for downstream marker queries.
"""
import os, time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

t0 = time.time()
DATA = "${PROJECT_ROOT}/data/processed"
RES = "${PROJECT_ROOT}/results"
FIG = f"{RES}/figures"

# load own_raw (full 9698 genes) and own_annotated (HVG 3000 + obs)
print("[load] own_raw + own_annotated")
own = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad")
ann = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")

# transfer obs cols from annotated → raw
for col in ["scanvi_predicted", "scanvi_uncertainty", "leiden_0.6", "leiden_1.0"]:
    own.obs[col] = ann.obs.loc[own.obs.index, col].values
print(f"  own_raw: {own.shape}; merged obs cols: {[c for c in own.obs.columns if c.startswith(('scanvi','leiden'))]}")

# normalize + log
sc.pp.normalize_total(own, target_sum=1e4)
sc.pp.log1p(own)

# Comprehensive TAN/NAN literature panel
# TAN classification draws from Salcher 2022 + Wu 2023 PNAS + Zilionis 2019 Immunity
panels = {
    "TAN1_VEGFA":     ["VEGFA", "ICAM1", "CCRL2", "PROK2", "BCL2A1", "CCL3", "CCL4"],
    "TAN2_IFN":       ["IFIT1", "IFIT2", "IFIT3", "ISG15", "RSAD2", "MX1", "OAS1", "IFI6", "IFITM3"],
    "TAN3_aged":      ["CXCR4", "CXCL8", "IL1B", "PTGS2", "G0S2", "CCL3", "TNF"],
    "TAN4_NETs":      ["PADI4", "MPO", "ELANE", "S100A8", "S100A9", "DEFA3", "DEFA1B",
                       "AZU1", "PRTN3", "CTSG", "BPI", "LTF", "LCN2"],
    "NAN1_resident":  ["FCGR3B", "CXCR2", "CSF3R", "IFITM2", "FPR1", "FCGR2A"],
    "NAN2_home":      ["S100A12", "VCAN", "CD55", "PROK2", "MNDA", "FCGR3B"],
    "NAN3_prog":      ["LCN2", "MMP8", "MMP9", "BPI", "CAMP", "ARG1", "DEFA4", "RETN"],
    # core neutrophil identity (sanity)
    "Neutro_core":    ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "FPR1",
                       "MPO", "ELANE", "MNDA", "CXCR4"],
    # NETosis canonical
    "NETs_canonical": ["PADI4", "MPO", "ELANE", "H2AC1", "H2BC1", "AZU1", "PRTN3", "CTSG"],
    # Inflammation
    "Inflam":         ["TNF", "IL1B", "CXCL8", "IL6", "PTGS2", "NFKBIA", "ICAM1"],
    # Immunosuppressive
    "Immuno_supp":    ["ARG1", "CD274", "LGALS3", "VEGFA", "TGFB1"],
}

vn = set(own.var_names.astype(str))
report_rows = []
for name, genes in panels.items():
    have = [g for g in genes if g in vn]
    miss = [g for g in genes if g not in vn]
    report_rows.append({"panel": name, "n_genes": len(genes),
                        "n_present": len(have),
                        "missing": ",".join(miss),
                        "present": ",".join(have)})
    if have:
        sc.tl.score_genes(own, gene_list=have,
                          score_name=f"score_{name}",
                          ctrl_size=max(20, len(have)*5),
                          random_state=0)
report = pd.DataFrame(report_rows)
report.to_csv(f"{RES}/step25e_panel_gene_availability.csv", index=False)
print("[panel availability]")
print(report.to_string(index=False))

# panel score by scanvi_predicted + by leiden
score_cols = [c for c in own.obs.columns if c.startswith("score_")]
mean_scanvi = own.obs.groupby("scanvi_predicted", observed=True)[score_cols].mean()
mean_leiden = own.obs.groupby("leiden_0.6", observed=True)[score_cols].mean()

# z-score per panel (column) for cleaner heatmap
def zscore(df):
    return (df - df.mean()) / df.std().replace(0, 1)

mean_scanvi.to_csv(f"{RES}/step25e_panel_score_byscanvi.csv")
mean_leiden.to_csv(f"{RES}/step25e_panel_score_byleiden.csv")

# Big heatmap
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
sns.heatmap(zscore(mean_scanvi).T, ax=axes[0], cmap="RdBu_r", center=0,
            annot=mean_scanvi.T, fmt=".2f", annot_kws={"size": 7},
            cbar_kws={"label": "z-score (by panel)"})
axes[0].set_title("Panel score by scanvi_predicted (annot=raw mean)")
sns.heatmap(zscore(mean_leiden).T, ax=axes[1], cmap="RdBu_r", center=0,
            annot=mean_leiden.T, fmt=".2f", annot_kws={"size": 7},
            cbar_kws={"label": "z-score (by panel)"})
axes[1].set_title("Panel score by leiden_0.6 (annot=raw mean)")
fig.tight_layout()
fig.savefig(f"{FIG}/step25e_panel_heatmap_full.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25e_panel_heatmap_full.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# DE on full gene set per scanvi label, save
print("[de] full-gene rank_genes_groups by scanvi_predicted")
sc.tl.rank_genes_groups(own, "scanvi_predicted", method="wilcoxon", n_genes=80)
def rgg_to_df(adata, key):
    rg = adata.uns[key]
    out = []
    for grp in rg["names"].dtype.names:
        df = pd.DataFrame({
            "group": grp,
            "rank": np.arange(len(rg["names"][grp])),
            "gene": rg["names"][grp],
            "score": rg["scores"][grp],
            "logfc": rg["logfoldchanges"][grp],
            "pval": rg["pvals"][grp],
            "padj": rg["pvals_adj"][grp],
        })
        out.append(df)
    return pd.concat(out, ignore_index=True)
rgg_to_df(own, "rank_genes_groups").to_csv(
    f"{RES}/step25e_markers_scanvi_fullgene.csv", index=False)

print("\n[de] full-gene rank_genes_groups by leiden_0.6")
sc.tl.rank_genes_groups(own, "leiden_0.6", method="wilcoxon", n_genes=80)
rgg_to_df(own, "rank_genes_groups").to_csv(
    f"{RES}/step25e_markers_leiden_fullgene.csv", index=False)

# update annotated h5ad with the panel scores (transfer back)
print("[merge] writing panel scores back to annotated h5ad")
for c in score_cols:
    ann.obs[c] = own.obs.loc[ann.obs.index, c].values
ann.write_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad", compression="gzip")
print(f"  saved.")

# tag summary
print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
