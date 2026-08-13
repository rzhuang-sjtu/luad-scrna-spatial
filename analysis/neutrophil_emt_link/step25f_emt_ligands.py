"""EMT ligand analysis: TAN/NAN × EMT secretome → setting up Figure 5D + tissue panel.

For each scANVI label (TAN-1..4, NAN-1..3, Neutrophils):
  1. mean / median expression per ligand
  2. dot plot (Wilcoxon DE significance + mean expression)
  3. fraction-expressing per ligand (CellChat-style cutoff)
  4. tissue distribution (TAN/NAN × tissue_type, stacked bar + crosstab)
  5. CSVs for downstream R plotting / Figure 5D-E

Logic chain we are setting up:
  TAN subtype (esp. TAN-4 NETs / TAN-1 VEGFA+) → EMT-promoting ligands (TGFB1, OSM, IL6, TNF, MMP9, SPP1)
    → MP3 (EMT/IFN program in malignant cells, shown in Fig 3) → invasion, metastasis, poor survival.

Reads:  data/processed/luad_neutrophil_own_annotated.h5ad
        data/processed/luad_neutrophil_own_raw.h5ad   (full 9698 genes for ligand expression)
Writes: results/step25f_emt_ligand_*.csv, figures/step25f_*.pdf/png
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

# load
print("[load] annotated (HVG-3000) + raw (9698)")
ann = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")
raw = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad")
# carry labels from annotated to raw
for col in ["scanvi_predicted", "scanvi_uncertainty", "leiden_0.6", "leiden_1.0"]:
    raw.obs[col] = ann.obs.loc[raw.obs.index, col].values
print(f"  raw: {raw.shape}; merged scanvi labels OK")

# normalize + log raw
sc.pp.normalize_total(raw, target_sum=1e4)
sc.pp.log1p(raw)

EMT_LIGANDS = {
    "TGFb_family": ["TGFB1", "TGFB2", "TGFB3"],
    "TNF_IL_axis": ["TNF", "IL6", "IL1B", "OSM", "IL1A"],
    "Chemokines":  ["CXCL8", "CXCL1", "CXCL2", "CXCL5", "CCL2", "CCL3", "CCL4", "CCL5"],
    "MMP_remodel": ["MMP9", "MMP2", "MMP8", "MMP25"],
    "Angio_EMT":   ["VEGFA", "VEGFB", "FN1", "SPP1", "SERPINE1", "PLAU", "PLAUR"],
    "Other_pro_EMT": ["HGF", "AREG", "EREG", "WNT5A", "PDGFB", "FGF2"],
}
ALL_EMT = sorted({g for v in EMT_LIGANDS.values() for g in v})
print(f"\n[panel] {len(ALL_EMT)} unique EMT ligands queried")

vn = set(raw.var_names.astype(str))
avail = {k: [g for g in v if g in vn] for k, v in EMT_LIGANDS.items()}
miss  = {k: [g for g in v if g not in vn] for k, v in EMT_LIGANDS.items()}
flat_avail = sorted({g for v in avail.values() for g in v})
flat_miss  = sorted({g for v in miss.values() for g in v})

# availability report
avail_rows = []
for k, all_g in EMT_LIGANDS.items():
    avail_rows.append({"family": k, "n_total": len(all_g),
                       "n_present": len(avail[k]),
                       "present": ",".join(avail[k]),
                       "missing": ",".join(miss[k])})
pd.DataFrame(avail_rows).to_csv(f"{RES}/step25f_ligand_availability.csv", index=False)
print(f"  available: {len(flat_avail)} / {len(ALL_EMT)}")
print(f"  missing:   {flat_miss}")

print("\n[1] mean / fraction-expressing per scANVI label")

# pick a stable order for scANVI labels
SCANVI_ORDER = ["TAN-1", "TAN-2", "TAN-3", "TAN-4",
                "NAN-1", "NAN-2", "NAN-3", "Neutrophils"]
present_order = [x for x in SCANVI_ORDER if x in raw.obs["scanvi_predicted"].unique()]
raw.obs["scanvi_predicted"] = pd.Categorical(raw.obs["scanvi_predicted"],
                                              categories=present_order, ordered=True)

X = raw[:, flat_avail].X
if hasattr(X, "toarray"):
    X = X.toarray()
mean_df = pd.DataFrame(X, index=raw.obs.index, columns=flat_avail)
mean_df["scanvi_predicted"] = raw.obs["scanvi_predicted"].values
mean_by_scanvi = mean_df.groupby("scanvi_predicted", observed=True)[flat_avail].mean()

# fraction-expressing (binary: x > 0)
binX = (X > 0).astype(np.float32)
frac_df = pd.DataFrame(binX, index=raw.obs.index, columns=flat_avail)
frac_df["scanvi_predicted"] = raw.obs["scanvi_predicted"].values
frac_by_scanvi = frac_df.groupby("scanvi_predicted", observed=True)[flat_avail].mean()

mean_by_scanvi.to_csv(f"{RES}/step25f_ligand_mean_byscanvi.csv")
frac_by_scanvi.to_csv(f"{RES}/step25f_ligand_frac_byscanvi.csv")

# also by leiden
mean_by_leiden = mean_df.assign(leiden=raw.obs["leiden_0.6"].values).groupby("leiden", observed=True)[flat_avail].mean()
frac_by_leiden = frac_df.assign(leiden=raw.obs["leiden_0.6"].values).groupby("leiden", observed=True)[flat_avail].mean()
mean_by_leiden.to_csv(f"{RES}/step25f_ligand_mean_byleiden.csv")
frac_by_leiden.to_csv(f"{RES}/step25f_ligand_frac_byleiden.csv")

print(f"  mean_by_scanvi shape: {mean_by_scanvi.shape}")
print("\n  mean expression (top-3 by row):")
top3 = mean_by_scanvi.apply(lambda r: ",".join(r.nlargest(3).index), axis=1)
print(top3.to_string())

print("\n[2] dot plot")
# group ligands in family order on x-axis
gene_order = []
for k in ["TGFb_family", "TNF_IL_axis", "Chemokines", "MMP_remodel", "Angio_EMT", "Other_pro_EMT"]:
    gene_order += [g for g in EMT_LIGANDS[k] if g in vn]

sc.settings.set_figure_params(dpi=150, frameon=False, fontsize=9)
ax = sc.pl.dotplot(
    raw, var_names=gene_order, groupby="scanvi_predicted",
    standard_scale="var", dendrogram=False, swap_axes=False,
    categories_order=present_order, return_fig=True,
    cmap="Reds", figsize=(0.32 * len(gene_order), 0.45 * len(present_order)),
)
ax.savefig(f"{FIG}/step25f_emt_dotplot_scanvi.pdf", bbox_inches="tight")
ax.savefig(f"{FIG}/step25f_emt_dotplot_scanvi.png", dpi=150, bbox_inches="tight")
plt.close()

# also leiden version
ax2 = sc.pl.dotplot(
    raw, var_names=gene_order, groupby="leiden_0.6",
    standard_scale="var", dendrogram=True, swap_axes=False,
    return_fig=True, cmap="Reds",
    figsize=(0.32 * len(gene_order), 0.45 * raw.obs["leiden_0.6"].nunique()),
)
ax2.savefig(f"{FIG}/step25f_emt_dotplot_leiden.pdf", bbox_inches="tight")
ax2.savefig(f"{FIG}/step25f_emt_dotplot_leiden.png", dpi=150, bbox_inches="tight")
plt.close()

print("\n[3] DE: scanvi_predicted")
sc.tl.rank_genes_groups(raw, "scanvi_predicted", method="wilcoxon", n_genes=200, pts=True)

def rgg_long(adata, key="rank_genes_groups"):
    rg = adata.uns[key]
    out = []
    for grp in rg["names"].dtype.names:
        df = pd.DataFrame({
            "group": grp,
            "rank": np.arange(len(rg["names"][grp])),
            "gene": rg["names"][grp],
            "logfc": rg["logfoldchanges"][grp],
            "pval": rg["pvals"][grp],
            "padj": rg["pvals_adj"][grp],
        })
        out.append(df)
    return pd.concat(out, ignore_index=True)
de = rgg_long(raw)
de_emt = de[de["gene"].isin(flat_avail)].copy()
de_emt = de_emt.sort_values(["group", "padj"]).reset_index(drop=True)
de_emt.to_csv(f"{RES}/step25f_de_emt_byscanvi.csv", index=False)

# significant up-regulated EMT ligands per scANVI label (padj<0.05, logfc>0.5)
sig = de_emt[(de_emt["padj"] < 0.05) & (de_emt["logfc"] > 0.5)].copy()
sig_pivot = sig.pivot_table(index="gene", columns="group", values="logfc", fill_value=0)
sig_pivot.to_csv(f"{RES}/step25f_emt_logfc_matrix.csv")
print("  significant up EMT ligands per scANVI label (padj<0.05, lfc>0.5):")
print(sig.groupby("group")["gene"].apply(list).to_string())

# heatmap of logfc
fig, ax = plt.subplots(figsize=(6, 0.35 * len(sig_pivot.index) + 1.5))
sns.heatmap(sig_pivot.reindex(columns=present_order, fill_value=0),
            cmap="RdBu_r", center=0, annot=True, fmt=".2f",
            cbar_kws={"label": "logFC vs rest"}, ax=ax)
ax.set_title("EMT ligand logFC (significant up, padj<0.05, lfc>0.5)\nby scANVI label")
fig.tight_layout()
fig.savefig(f"{FIG}/step25f_emt_logfc_heatmap.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25f_emt_logfc_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n[4] tissue distribution")
tissue_xtab = pd.crosstab(raw.obs["scanvi_predicted"], raw.obs["tissue_type"])
tissue_xtab_col = pd.crosstab(raw.obs["scanvi_predicted"], raw.obs["tissue_type"], normalize="columns")
tissue_xtab_row = pd.crosstab(raw.obs["scanvi_predicted"], raw.obs["tissue_type"], normalize="index")
tissue_xtab.to_csv(f"{RES}/step25f_tissue_xtab_raw.csv")
tissue_xtab_col.to_csv(f"{RES}/step25f_tissue_xtab_colnorm.csv")
tissue_xtab_row.to_csv(f"{RES}/step25f_tissue_xtab_rownorm.csv")
print("  scANVI × tissue (column-normalized — fraction of each tissue that is each label):")
print(tissue_xtab_col.round(3).to_string())

# stacked bar: per-tissue scANVI composition (column norm)
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
tissue_xtab_col.T.plot(kind="bar", stacked=True, ax=axes[0],
                       colormap="tab10", width=0.85, edgecolor="white", linewidth=0.5)
axes[0].set_title("scANVI composition per tissue (column-normalized)")
axes[0].set_ylabel("fraction")
axes[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=30, ha="right")

tissue_xtab_row.plot(kind="bar", stacked=True, ax=axes[1],
                    colormap="tab20", width=0.85, edgecolor="white", linewidth=0.5)
axes[1].set_title("Tissue composition per scANVI label (row-normalized)")
axes[1].set_ylabel("fraction")
axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=30, ha="right")
fig.tight_layout()
fig.savefig(f"{FIG}/step25f_tissue_distribution.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25f_tissue_distribution.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n[5] tumor-vs-normal ligand contrast")
TUMOR_TT = ["Primary_Tumor", "Brain_Metastasis", "LN_Metastasis", "Distant_Metastasis", "Pleural_Effusion"]
NORMAL_TT = ["Normal_Lung", "Adjacent_Normal", "Normal_LN"]
raw.obs["tissue_group"] = raw.obs["tissue_type"].astype(str).map(
    lambda x: "tumor" if x in TUMOR_TT else ("normal" if x in NORMAL_TT else "other"))
print(f"  tissue_group counts: {raw.obs['tissue_group'].value_counts().to_dict()}")

# per scANVI × tissue_group mean
mat = mean_df.assign(tissue_group=raw.obs["tissue_group"].values).groupby(
    ["scanvi_predicted", "tissue_group"], observed=True)[flat_avail].mean()
mat.to_csv(f"{RES}/step25f_ligand_mean_byscanvi_tissuegroup.csv")

# Tumor vs Normal logFC of mean ligand expression per scANVI label
tumor = mat.xs("tumor", level=1)
normal = mat.xs("normal", level=1) if "normal" in raw.obs["tissue_group"].unique() else None
if normal is not None:
    delta = (tumor + 1e-3).div(normal + 1e-3).apply(np.log2)
    delta = delta.reindex(present_order)
    delta.to_csv(f"{RES}/step25f_ligand_log2_tumor_vs_normal.csv")
    fig, ax = plt.subplots(figsize=(0.4 * len(flat_avail) + 1.5, 0.45 * len(present_order) + 1.5))
    sns.heatmap(delta, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "log2(tumor+ε / normal+ε)"}, ax=ax)
    ax.set_title("EMT ligand log2 fold (tumor vs normal) per TAN/NAN label")
    fig.tight_layout()
    fig.savefig(f"{FIG}/step25f_emt_tumor_vs_normal.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/step25f_emt_tumor_vs_normal.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

print("\n[merge] adding EMT ligand expression as obs cols (for cross-cell-type plots later)")
# also save log-normalized EMT ligand per cell as obsm  for fast retrieval
ann.obsm["emt_ligand_lognorm"] = mean_df.loc[ann.obs.index, flat_avail].values
ann.uns["emt_ligand_genes"] = flat_avail
ann.write_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad", compression="gzip")

print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
