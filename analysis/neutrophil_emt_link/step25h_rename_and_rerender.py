"""step25h: apply functional rename to TAN/NAN labels everywhere + re-render Figure 5 panels.

Rename map (per user, 2026-04-27):
  TAN-1 → Neu_Inflammatory      (high CXCL8/IL1B/PLAU)
  TAN-2 → Neu_IFN_response      (Salcher original IFN type; sparse in our data)
  TAN-3 → Neu_Angiogenic        (high VEGFA/CXCL2)
  TAN-4 → Neu_Metastatic        (enriched in LN/Brain/Pleural)
  NAN-1 → Neu_ECM_remodeling    (MMP9 specialist)
  NAN-2 → Neu_OSM_priming       (high OSM/PLAUR)
  NAN-3 → Neu_OSM_low           (OSM weaker version)
  Neutrophils → Neu_unclassified (low-quality / transitional)

Operations:
  1. Add `neu_subtype` (categorical, ordered) to annotated h5ad and own_raw h5ad.
     Keep `scanvi_predicted` for archeology.
  2. Re-render selected PDFs from existing CSVs with new names:
       - step26_tan_mp_correlation_matrix → step26b_*
       - step26_scatter_tan_mp3 → step26b_scatter_neu_mp_axes
       - step25f_tissue_distribution → step25f2_*
       - step27 LIANA focus heatmap + dotplots → step27b_*
  3. Apply rename inside CSVs (write *_renamed.csv with both old and new label cols
     for downstream R-side plotting consistency).
"""
import os, time
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

t0 = time.time()
DATA = "${PROJECT_ROOT}/data/processed"
RES = "${PROJECT_ROOT}/results"
FIG = f"{RES}/figures"

RENAME = {
    "TAN-1": "Neu_Inflammatory",
    "TAN-2": "Neu_IFN_response",
    "TAN-3": "Neu_Angiogenic",
    "TAN-4": "Neu_Metastatic",
    "NAN-1": "Neu_ECM_remodeling",
    "NAN-2": "Neu_OSM_priming",
    "NAN-3": "Neu_OSM_low",
    "Neutrophils": "Neu_unclassified",
}
ORDER_NEW = ["Neu_Inflammatory", "Neu_Angiogenic", "Neu_Metastatic", "Neu_IFN_response",
             "Neu_OSM_priming", "Neu_OSM_low", "Neu_ECM_remodeling", "Neu_unclassified"]
# colors: keep TAN-side in reds/oranges, NAN-side in blues, unclassified gray
PALETTE = {
    "Neu_Inflammatory":   "#d62728",  # crimson
    "Neu_Angiogenic":     "#ff7f0e",  # orange
    "Neu_Metastatic":     "#9467bd",  # purple
    "Neu_IFN_response":   "#bcbd22",  # olive
    "Neu_OSM_priming":    "#1f77b4",  # blue
    "Neu_OSM_low":        "#17becf",  # cyan
    "Neu_ECM_remodeling": "#2ca02c",  # green
    "Neu_unclassified":   "#7f7f7f",  # gray
}

def add_neu_subtype(adata):
    if "scanvi_predicted" not in adata.obs.columns:
        return False
    new = adata.obs["scanvi_predicted"].astype(str).map(RENAME)
    adata.obs["neu_subtype"] = pd.Categorical(new.values, categories=ORDER_NEW, ordered=True)
    return True

print("[1] update annotated h5ad")
ann = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")
add_neu_subtype(ann)
print(f"  neu_subtype counts:\n{ann.obs['neu_subtype'].value_counts().to_string()}")
ann.write_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad", compression="gzip")

print("\n[2] update own_raw h5ad")
raw = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad")
raw.obs["scanvi_predicted"] = ann.obs.loc[raw.obs.index, "scanvi_predicted"].values
add_neu_subtype(raw)
raw.write_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad", compression="gzip")

# write rename map for reference
pd.DataFrame([{"scanvi_label": k, "neu_subtype": v} for k, v in RENAME.items()]).to_csv(
    f"{RES}/step25h_rename_map.csv", index=False)

print("\n[3] re-render step26 correlation matrix")
combined = pd.read_csv(f"{RES}/step26_patient_combined.csv", index_col=0)
# rename frac_* columns
combined.columns = [c.replace("frac_TAN-1", f"frac_{RENAME['TAN-1']}")
                     .replace("frac_TAN-2", f"frac_{RENAME['TAN-2']}")
                     .replace("frac_TAN-3", f"frac_{RENAME['TAN-3']}")
                     .replace("frac_TAN-4", f"frac_{RENAME['TAN-4']}")
                     .replace("frac_NAN-1", f"frac_{RENAME['NAN-1']}")
                     .replace("frac_NAN-2", f"frac_{RENAME['NAN-2']}")
                     .replace("frac_NAN-3", f"frac_{RENAME['NAN-3']}")
                     .replace("frac_Neutrophils", f"frac_{RENAME['Neutrophils']}")
                    for c in combined.columns]
combined.to_csv(f"{RES}/step26b_patient_combined_renamed.csv")

from scipy.stats import spearmanr
mp_cols = ["MP1", "MP2", "MP3", "MP4", "MP5"]
frac_cols_new = [f"frac_{n}" for n in ORDER_NEW]
cor_rows = []
for fc in frac_cols_new:
    if fc not in combined.columns: continue
    nm = fc.replace("frac_", "")
    for mc in mp_cols:
        rho, p = spearmanr(combined[fc], combined[mc])
        cor_rows.append({"neu_subtype": nm, "MP": mc, "spearman_rho": rho, "p": p})
cor_df = pd.DataFrame(cor_rows)
cor_df.to_csv(f"{RES}/step26b_neu_mp_spearman_renamed.csv", index=False)

rho_mat = cor_df.pivot_table(index="neu_subtype", columns="MP",
                              values="spearman_rho").reindex(ORDER_NEW)
p_mat = cor_df.pivot_table(index="neu_subtype", columns="MP",
                            values="p").reindex(ORDER_NEW)
def annot_star(rho, p):
    s = f"{rho:.2f}"
    if p < 0.001: s += "***"
    elif p < 0.01: s += "**"
    elif p < 0.05: s += "*"
    return s
annot = pd.DataFrame(np.vectorize(annot_star)(rho_mat.values, p_mat.values),
                     index=rho_mat.index, columns=rho_mat.columns)
fig, ax = plt.subplots(figsize=(5.5, 5.5))
sns.heatmap(rho_mat, cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
            annot=annot, fmt="", cbar_kws={"label": "Spearman ρ"}, ax=ax,
            annot_kws={"size": 9})
ax.set_title(f"Neu subtype × malignant MP score\n(per-patient, n={len(combined)} samples)")
fig.tight_layout()
fig.savefig(f"{FIG}/step26b_neu_mp_correlation_matrix.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step26b_neu_mp_correlation_matrix.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# scatter — both axes
print("\n[3b] re-render step26 scatters with functional names")
fig, axes = plt.subplots(2, 4, figsize=(20, 9))
for ax, lbl in zip(axes.ravel(), ORDER_NEW):
    fc = f"frac_{lbl}"
    if fc not in combined.columns: continue
    x = combined[fc].values
    # primary axis: MP1 (immune/EMT) for inflammatory/angiogenic; MP2 (cycling) for metastatic
    primary_mp = "MP1" if lbl in ["Neu_Inflammatory", "Neu_Angiogenic", "Neu_OSM_priming",
                                  "Neu_OSM_low", "Neu_ECM_remodeling"] else "MP2"
    y = combined[primary_mp].values
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=combined["n_neutrophils"]/10, alpha=0.7,
               c=combined["n_malignant"], cmap="viridis", edgecolor="black", linewidth=0.4)
    if len(x) > 5:
        z = np.polyfit(x, y, 1)
        xr = np.linspace(x.min(), x.max(), 50)
        ax.plot(xr, np.polyval(z, xr), color=PALETTE[lbl], lw=2.0, alpha=0.9)
    ax.set_xlabel(f"frac of {lbl}")
    ax.set_ylabel(f"mean malignant {primary_mp}")
    ax.set_title(f"{lbl} ↔ {primary_mp}\nρ={rho:.2f}, p={p:.2g}", color=PALETTE[lbl])
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Per-patient Neu-subtype composition vs malignant MP score (primary axis)",
             y=1.02, fontsize=13)
fig.tight_layout()
fig.savefig(f"{FIG}/step26b_scatter_neu_primary_axes.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step26b_scatter_neu_primary_axes.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n[4] re-render step27 focus heatmap")
focus = pd.read_csv(f"{RES}/step27_liana_focus_pathways.csv")
print(f"  rows: {len(focus)}; cols: {list(focus.columns)[:10]}")
NEU_RENAME = {f"Neu_{k}": f"Neu_{v.replace('Neu_', '')}".replace("Neu_Neu_", "Neu_") for k, v in RENAME.items()}
# replace in source/target — note source has "Neu_TAN-1" → want "Neu_Inflammatory"
def remap(s):
    s = str(s)
    if s.startswith("Neu_"):
        old = s[4:]
        if old in RENAME: return RENAME[old]
    return s
focus["source_renamed"] = focus["source"].apply(remap)
focus["target_renamed"] = focus["target"]   # malignant labels Mal_MP1..4 unchanged
focus.to_csv(f"{RES}/step27b_liana_focus_pathways_renamed.csv", index=False)

# focus heatmap with new names
lig_col = "ligand_complex" if "ligand_complex" in focus.columns else "ligand"
rec_col = "receptor_complex" if "receptor_complex" in focus.columns else "receptor"
key_score = "magnitude_rank" if "magnitude_rank" in focus.columns else "lr_means"
focus["lr"] = focus[lig_col].astype(str) + "→" + focus[rec_col].astype(str)
focus["sender_target"] = focus["source_renamed"] + " → " + focus["target_renamed"].str.replace("Mal_", "")
pivot = focus.pivot_table(index="lr", columns="sender_target", values=key_score, aggfunc="min")

# order columns: cluster by sender (Neu_Inflammatory first, etc.) × MP1-4
sender_order = [n for n in ORDER_NEW if n != "Neu_unclassified"]
target_order = ["MP1", "MP2", "MP3", "MP4"]
ord_cols = [f"{s} → {t}" for s in sender_order for t in target_order if f"{s} → {t}" in pivot.columns]
pivot = pivot[ord_cols]

fig, ax = plt.subplots(figsize=(0.50 * pivot.shape[1] + 4, 0.32 * pivot.shape[0] + 2))
sns.heatmap(pivot, cmap="viridis_r", ax=ax,
            cbar_kws={"label": f"{key_score} (lower = stronger)"})
ax.set_title("Neu subtype → Malignant MP — focus EMT pathway LR strength")
ax.tick_params(axis='x', labelrotation=45)
plt.setp(ax.get_xticklabels(), ha="right")
fig.tight_layout()
fig.savefig(f"{FIG}/step27b_liana_focus_heatmap_renamed.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step27b_liana_focus_heatmap_renamed.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# also rename the per-sender top-30 CSVs
print("\n[4b] rename per-sender top-30 CSVs")
import glob
for f in glob.glob(f"{RES}/step27_liana_top30_*_to_malignant.csv"):
    df = pd.read_csv(f)
    # find old label from filename
    base = os.path.basename(f).replace("step27_liana_top30_", "").replace("_to_malignant.csv", "")
    # base is e.g. TAN1 / TAN4 / NAN2 etc. Need to add the dash back.
    base_dash = base[:3] + "-" + base[3:]
    new_label = RENAME.get(base_dash, base)
    df["source_renamed"] = df["source"].apply(remap)
    df["target_renamed"] = df["target"]
    df.to_csv(f"{RES}/step27b_liana_top30_{new_label}_to_malignant.csv", index=False)
    print(f"  {base_dash} → {new_label}: {f} → step27b_*")

print("\n[5] re-render tissue distribution stacked bar")
tissue_col = pd.read_csv(f"{RES}/step25f_tissue_xtab_colnorm.csv", index_col=0)
tissue_row = pd.read_csv(f"{RES}/step25f_tissue_xtab_rownorm.csv", index_col=0)
# reindex with renamed
tissue_col.index = [RENAME.get(x, x) for x in tissue_col.index]
tissue_row.index = [RENAME.get(x, x) for x in tissue_row.index]
tissue_col = tissue_col.reindex(ORDER_NEW)
tissue_row = tissue_row.reindex(ORDER_NEW)
tissue_col.to_csv(f"{RES}/step25f2_tissue_xtab_colnorm_renamed.csv")
tissue_row.to_csv(f"{RES}/step25f2_tissue_xtab_rownorm_renamed.csv")

fig, axes = plt.subplots(1, 2, figsize=(15, 6))
colors_list = [PALETTE[x] for x in tissue_col.index]
tissue_col.T.plot(kind="bar", stacked=True, ax=axes[0],
                  color=colors_list, width=0.85, edgecolor="white", linewidth=0.5)
axes[0].set_title("Neu subtype composition per tissue")
axes[0].set_ylabel("fraction")
axes[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=30, ha="right")

tissue_row.plot(kind="bar", stacked=True, ax=axes[1],
                colormap="tab20", width=0.85, edgecolor="white", linewidth=0.5)
axes[1].set_title("Tissue composition per Neu subtype")
axes[1].set_ylabel("fraction")
axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=30, ha="right")
fig.tight_layout()
fig.savefig(f"{FIG}/step25f2_tissue_distribution_renamed.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25f2_tissue_distribution_renamed.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("\n[6] re-render EMT dot plot")
ann.obs["neu_subtype"] = pd.Categorical(
    ann.obs["scanvi_predicted"].astype(str).map(RENAME),
    categories=ORDER_NEW, ordered=True)
sc.pp.normalize_total(ann, target_sum=1e4) if "lognorm" not in ann.layers else None
if "lognorm" not in ann.layers and (ann.X.max() > 50):
    sc.pp.normalize_total(ann, target_sum=1e4); sc.pp.log1p(ann)

# load gene order from step25f
EMT_GENES = ['TGFB1','TNF','IL6','IL1B','OSM','IL1A','CXCL8','CXCL1','CXCL2','CCL2',
             'CCL3','CCL4','CCL5','MMP9','VEGFA','VEGFB','FN1','SPP1','SERPINE1','PLAU','PLAUR',
             'AREG','EREG','WNT5A','PDGFB']
gene_order = [g for g in EMT_GENES if g in ann.var_names]
ax = sc.pl.dotplot(
    ann, var_names=gene_order, groupby="neu_subtype",
    standard_scale="var", dendrogram=False,
    categories_order=ORDER_NEW, return_fig=True,
    cmap="Reds", figsize=(0.32*len(gene_order), 0.45*len(ORDER_NEW)),
)
ax.savefig(f"{FIG}/step25f2_emt_dotplot_renamed.pdf", bbox_inches="tight")
ax.savefig(f"{FIG}/step25f2_emt_dotplot_renamed.png", dpi=150, bbox_inches="tight")
plt.close()

print("\n[7] re-render UMAP with functional colors")
sc.settings.set_figure_params(dpi=150, frameon=False)
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
sc.pl.umap(ann, color="neu_subtype", ax=axes[0], show=False, frameon=False,
           palette=PALETTE, title="Neu subtype (functional)",
           legend_loc="right margin", legend_fontsize=8)
sc.pl.umap(ann, color="leiden_0.6", ax=axes[1], show=False, frameon=False, title="Leiden res=0.6")
sc.pl.umap(ann, color="dataset", ax=axes[2], show=False, frameon=False, title="Dataset")
fig.tight_layout()
fig.savefig(f"{FIG}/step25h_umap_neu_functional.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step25h_umap_neu_functional.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
