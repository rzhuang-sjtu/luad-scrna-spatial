"""Downstream of step25c: figures, markers, and summary.

Reads:
  data/processed/luad_neutrophil_own_annotated.h5ad
  data/processed/luad_neutrophil_joint_scanvi.h5ad

Writes:
  results/figures/step25_neutrophil_umap_own.pdf       (own cells, by leiden + by scanvi_pred + by dataset)
  results/figures/step25_neutrophil_umap_joint.pdf     (joint Salcher+own, side-by-side)
  results/step25_neutrophil_markers_leiden.csv         (DE markers per leiden cluster)
  results/step25_neutrophil_markers_scanvi.csv         (DE markers per scanvi_predicted)
  results/step25_neutrophil_confusion.csv              (leiden × scanvi_predicted)
  results/step25_neutrophil_label_panel_score.csv      (literature TAN/NAN markers scored on own)
  results/step25_summary.md
"""
import os, time
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t0 = time.time()
DATA = "${PROJECT_ROOT}/data/processed"
RES = "${PROJECT_ROOT}/results"
FIG = f"{RES}/figures"
os.makedirs(FIG, exist_ok=True)

own = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")
joint = sc.read_h5ad(f"{DATA}/luad_neutrophil_joint_scanvi.h5ad")
print(f"own: {own.shape}; joint: {joint.shape}")
print(f"own.obs cols: {list(own.obs.columns)[:30]}")

# need lognorm layer for marker DE — own annotated h5ad only carries HVG counts.
# rebuild log1p of HVG counts in-place.
sc.pp.normalize_total(own, target_sum=1e4, layer=None)
sc.pp.log1p(own)
own.layers["lognorm"] = own.X.copy()

# 1. UMAPs — own cells
print("[fig] UMAP own")
sc.settings.set_figure_params(dpi=150, figsize=(5, 5), frameon=False)
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
sc.pl.umap(own, color="scanvi_predicted", ax=axes[0, 0], show=False, frameon=False,
           title="scANVI predicted (TAN/NAN)", legend_loc="on data", legend_fontsize=8)
sc.pl.umap(own, color="leiden_0.6", ax=axes[0, 1], show=False, frameon=False,
           title="Leiden res=0.6", legend_loc="on data", legend_fontsize=8)
sc.pl.umap(own, color="leiden_1.0", ax=axes[1, 0], show=False, frameon=False,
           title="Leiden res=1.0", legend_loc="on data", legend_fontsize=8)
sc.pl.umap(own, color="dataset", ax=axes[1, 1], show=False, frameon=False, title="Dataset")
fig.tight_layout()
fig.savefig(f"{FIG}/step25_neutrophil_umap_own.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25_neutrophil_umap_own.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 2. Joint UMAP — Salcher + own co-embedding
print("[fig] UMAP joint")
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
sc.pl.umap(joint, color="source", ax=axes[0], show=False, frameon=False, title="Source")
# True Salcher labels
sc.pl.umap(joint, color="transfer_label", ax=axes[1], show=False, frameon=False,
           title="True (Salcher) / Unknown (own)", legend_fontsize=7)
# scANVI predicted (Salcher = true, own = pred)
sc.pl.umap(joint, color="scanvi_predicted", ax=axes[2], show=False, frameon=False,
           title="scANVI predicted (all)", legend_fontsize=7)
fig.tight_layout()
fig.savefig(f"{FIG}/step25_neutrophil_umap_joint.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25_neutrophil_umap_joint.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 3. DE markers — per leiden + per scanvi predicted
print("[de] leiden_0.6 markers")
sc.tl.rank_genes_groups(own, "leiden_0.6", method="wilcoxon", n_genes=50, layer="lognorm", use_raw=False)
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
df_l = rgg_to_df(own, "rank_genes_groups")
df_l.to_csv(f"{RES}/step25_neutrophil_markers_leiden.csv", index=False)

print("[de] scanvi_predicted markers")
sc.tl.rank_genes_groups(own, "scanvi_predicted", method="wilcoxon", n_genes=50, layer="lognorm", use_raw=False)
df_s = rgg_to_df(own, "rank_genes_groups")
df_s.to_csv(f"{RES}/step25_neutrophil_markers_scanvi.csv", index=False)

# 4. Literature TAN/NAN marker panel score
#    (Salcher 2022 + Wu 2023 + Zilionis 2019)
print("[score] literature TAN/NAN panel")
panels = {
    "TAN-1_canonical":   ["VEGFA", "ICAM1", "CCRL2", "PROK2", "BCL2A1"],
    "TAN-2_IFN":         ["IFIT1", "IFIT2", "IFIT3", "ISG15", "RSAD2", "MX1", "OAS1"],
    "TAN-3_aged":        ["CXCR4", "CXCL8", "IL1B", "PTGS2", "G0S2"],
    "TAN-4_NETs":        ["PADI4", "MPO", "ELANE", "S100A8", "S100A9", "DEFA3", "DEFA1B"],
    "NAN-1_resident":    ["FCGR3B", "CXCR2", "CSF3R", "IFITM2"],
    "NAN-2_homeostatic": ["S100A12", "VCAN", "CD55"],
    "NAN-3_progenitor":  ["LCN2", "MMP8", "MMP9", "BPI"],
}
panel_present = {k: [g for g in v if g in own.var_names] for k, v in panels.items()}
for k, genes in panel_present.items():
    if genes:
        sc.tl.score_genes(own, gene_list=genes, score_name=f"score_{k}",
                          ctrl_size=min(50, max(10, len(genes)*5)), random_state=0)
score_cols = [c for c in own.obs.columns if c.startswith("score_")]
panel_df = own.obs.groupby("leiden_0.6")[score_cols].mean()
panel_df.to_csv(f"{RES}/step25_neutrophil_label_panel_score.csv")
panel_df_scanvi = own.obs.groupby("scanvi_predicted")[score_cols].mean()
panel_df_scanvi.to_csv(f"{RES}/step25_neutrophil_label_panel_score_byscanvi.csv")

# panel score heatmap
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
import seaborn as sns
sns.heatmap(panel_df.T, ax=axes[0], cmap="vlag", center=0, annot=True, fmt=".2f")
axes[0].set_title("Panel score by leiden_0.6")
sns.heatmap(panel_df_scanvi.T, ax=axes[1], cmap="vlag", center=0, annot=True, fmt=".2f")
axes[1].set_title("Panel score by scanvi_predicted")
fig.tight_layout()
fig.savefig(f"{FIG}/step25_neutrophil_panel_heatmap.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25_neutrophil_panel_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 5. Confusion matrix leiden × scanvi_predicted
print("[xtab] leiden × scanvi_predicted")
ct = pd.crosstab(own.obs["leiden_0.6"], own.obs["scanvi_predicted"])
ct.to_csv(f"{RES}/step25_neutrophil_confusion.csv")
ct_norm = ct.div(ct.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(ct_norm, ax=ax, cmap="Blues", annot=True, fmt=".2f")
ax.set_title("scANVI label fractions per leiden_0.6 cluster")
fig.tight_layout()
fig.savefig(f"{FIG}/step25_neutrophil_confusion.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25_neutrophil_confusion.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 6. Summary md
print("[md] summary")
ts = own.uns.get("transfer_summary", {})
md = []
md.append("# step25 — Neutrophil subclustering + TAN/NAN label transfer\n")
md.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
md.append("## Inputs")
md.append(f"- own neutrophils: `data/processed/luad_neutrophil_own_raw.h5ad` ({own.n_obs} cells)")
md.append(f"- Salcher reference: `${DATA_ROOT}/High-resolution/neutrophil_final.h5ad` (~19k neutrophils, 8 TAN/NAN labels)")
md.append("\n## Method")
md.append("- HGNC symbol harmonization (Salcher ENSG → symbol via `var['feature_name']`)")
md.append(f"- Shared genes: 9479; HVG (seurat_v3 on counts): {ts.get('n_hvg', 'N/A')}")
md.append("- Joint scVI (200 epochs early-stop) → scANVI fine-tune (25 epochs, n_samples_per_label=100)")
md.append("- batch_key = source-prefixed dataset (29 levels)")
md.append("- labels_key = `cell_type_neutro` for Salcher, `Unknown` for own")
md.append(f"\n## Outputs")
md.append(f"- own annotated: `data/processed/luad_neutrophil_own_annotated.h5ad`")
md.append(f"- joint: `data/processed/luad_neutrophil_joint_scanvi.h5ad`")
md.append("\n## Predicted label distribution (own 8.5k cells)")
md.append("```")
md.append(own.obs["scanvi_predicted"].value_counts().to_string())
md.append("```")
md.append(f"\nMean uncertainty: {ts.get('mean_uncertainty', float('nan')):.3f}")
md.append("\n## Leiden subclustering (independent)")
md.append(f"- res=0.6: {own.obs['leiden_0.6'].nunique()} clusters")
md.append(f"- res=1.0: {own.obs['leiden_1.0'].nunique()} clusters")
md.append("\n## Confusion (leiden_0.6 → scANVI top label)")
md.append("```")
top_per_l = ct.idxmax(axis=1)
agree = (ct.max(axis=1) / ct.sum(axis=1)).round(2)
md.append(pd.DataFrame({"top_label": top_per_l, "purity": agree}).to_string())
md.append("```")
md.append("\n## Files")
for f in [
    "results/figures/step25_neutrophil_umap_own.pdf",
    "results/figures/step25_neutrophil_umap_joint.pdf",
    "results/figures/step25_neutrophil_panel_heatmap.pdf",
    "results/figures/step25_neutrophil_confusion.pdf",
    "results/step25_neutrophil_markers_leiden.csv",
    "results/step25_neutrophil_markers_scanvi.csv",
    "results/step25_neutrophil_confusion.csv",
    "results/step25_neutrophil_label_panel_score.csv",
    "results/step25_neutrophil_label_panel_score_byscanvi.csv",
    "results/step25c_predictions.csv.gz",
]:
    md.append(f"- {f}")
with open(f"{RES}/step25_summary.md", "w") as fh:
    fh.write("\n".join(md) + "\n")

print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
