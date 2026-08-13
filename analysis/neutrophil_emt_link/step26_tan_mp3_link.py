"""step26: TAN/NAN × MP correlation per patient (Path A core).

Closes the causal loop: Fig 3 (MP3=EMT/IFN, drives poor survival) ↔ Fig 5 (TAN/NAN secretome).

Inputs:
  data/processed/luad_neutrophil_own_annotated.h5ad   (8549 neutrophils × scANVI labels)
  data/processed/luad_malignant_scored.h5ad           (45791 malignant cells × MP1-5 scores)

Outputs:
  results/step26_patient_tan_fractions.csv             (per patient × TAN/NAN fraction)
  results/step26_patient_mp_means.csv                  (per patient × mean MP score)
  results/step26_tan_mp_spearman.csv                   (cor matrix + p)
  results/step26_tan_signatures.csv                    (top-50 markers per TAN/NAN; reused in step28)
  figures/step26_tan_mp_correlation_matrix.pdf         (heatmap)
  figures/step26_scatter_tan_mp3.pdf                   (TAN-1/3/4 vs MP3 scatter triple)
  figures/step26_scatter_nan1_mmp9_axis.pdf            (NAN-1 vs MP3, MP1, etc)
  figures/step26_canonical_marker_dotplot.pdf          (Fig 5B candidate)
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
from scipy.stats import spearmanr, pearsonr

t0 = time.time()
DATA = "${PROJECT_ROOT}/data/processed"
RES = "${PROJECT_ROOT}/results"
FIG = f"{RES}/figures"
os.makedirs(FIG, exist_ok=True)

# ----- load -----
print("[load] neutrophil + malignant")
neu = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_annotated.h5ad")
mal = sc.read_h5ad(f"{DATA}/luad_malignant_scored.h5ad", backed="r")
print(f"  neu: {neu.shape}; mal: {mal.shape}")

# canonicalize sample_id sets
print(f"  neu n_samples: {neu.obs['sample_id'].nunique()}; mal n_samples: {mal.obs['sample_id'].nunique()}")
shared = sorted(set(neu.obs["sample_id"].astype(str)) & set(mal.obs["sample_id"].astype(str)))
print(f"  shared samples (have both neutrophil + malignant): {len(shared)}")

# minimum cells per sample (avoid noise)
MIN_NEU = 10
MIN_MAL = 20

# ----- per-patient TAN/NAN fractions -----
print(f"\n[1] per-patient TAN/NAN fractions (min {MIN_NEU} neutrophils + {MIN_MAL} malignant)")
neu_by = neu.obs.groupby("sample_id", observed=True)
SCANVI_LABELS = ["TAN-1", "TAN-2", "TAN-3", "TAN-4", "NAN-1", "NAN-2", "NAN-3", "Neutrophils"]
rows = []
for sid, sub in neu_by:
    if len(sub) < MIN_NEU or sid not in shared:
        continue
    fracs = sub["scanvi_predicted"].value_counts(normalize=True).to_dict()
    rows.append({"sample_id": sid, "n_neutrophils": len(sub),
                 **{f"frac_{lbl}": fracs.get(lbl, 0.0) for lbl in SCANVI_LABELS}})
tan_frac = pd.DataFrame(rows).set_index("sample_id")
tan_frac.to_csv(f"{RES}/step26_patient_tan_fractions.csv")
print(f"  shape: {tan_frac.shape}")

# ----- per-patient MP score means (over malignant cells) -----
print("\n[2] per-patient MP score means (over malignant cells)")
mal_obs = mal.obs[mal.obs["sample_id"].astype(str).isin(set(tan_frac.index))][
    ["sample_id", "MP1_score", "MP2_score", "MP3_score", "MP4_score", "MP5_score"]
].copy()
mp_mean = mal_obs.groupby("sample_id", observed=True).agg(
    n_malignant=("MP1_score", "size"),
    MP1=("MP1_score", "mean"), MP2=("MP2_score", "mean"),
    MP3=("MP3_score", "mean"), MP4=("MP4_score", "mean"), MP5=("MP5_score", "mean"),
)
mp_mean = mp_mean[mp_mean["n_malignant"] >= MIN_MAL]
mp_mean.to_csv(f"{RES}/step26_patient_mp_means.csv")
print(f"  shape: {mp_mean.shape}")

# ----- merge -----
combined = tan_frac.join(mp_mean, how="inner")
print(f"\n[3] merged per-patient table: {combined.shape}")
print(combined.head().to_string())
combined.to_csv(f"{RES}/step26_patient_combined.csv")

# ----- Spearman correlation matrix -----
print("\n[4] Spearman correlation TAN/NAN frac × MP")
frac_cols = [f"frac_{l}" for l in SCANVI_LABELS]
mp_cols = ["MP1", "MP2", "MP3", "MP4", "MP5"]
cor_rows = []
for fc in frac_cols:
    for mc in mp_cols:
        x = combined[fc].values
        y = combined[mc].values
        rho, p = spearmanr(x, y)
        cor_rows.append({"tan_label": fc.replace("frac_", ""), "MP": mc,
                         "spearman_rho": rho, "p": p, "n": len(combined)})
cor_df = pd.DataFrame(cor_rows)
# adj p (BH within MP3 column for primary hypothesis, then full)
from scipy.stats import false_discovery_control
try:
    cor_df["padj_full"] = false_discovery_control(cor_df["p"].values, method="bh")
except Exception:
    # fallback for older scipy
    from statsmodels.stats.multitest import multipletests
    cor_df["padj_full"] = multipletests(cor_df["p"], method="fdr_bh")[1]
cor_df.to_csv(f"{RES}/step26_tan_mp_spearman.csv", index=False)
print(cor_df.pivot_table(index="tan_label", columns="MP", values="spearman_rho").round(3).to_string())
print()
print("Top |rho| pairs (sorted):")
print(cor_df.reindex(cor_df["spearman_rho"].abs().sort_values(ascending=False).index).head(15).to_string(index=False))

# ----- heatmap -----
print("\n[5] heatmap TAN/NAN × MP")
rho_mat = cor_df.pivot_table(index="tan_label", columns="MP", values="spearman_rho").reindex(SCANVI_LABELS)
p_mat = cor_df.pivot_table(index="tan_label", columns="MP", values="p").reindex(SCANVI_LABELS)
# annotate with rho with star for sig
def annot(rho, p):
    s = f"{rho:.2f}"
    if p < 0.001: s += "***"
    elif p < 0.01: s += "**"
    elif p < 0.05: s += "*"
    return s
annot_mat = pd.DataFrame(
    np.vectorize(annot)(rho_mat.values, p_mat.values),
    index=rho_mat.index, columns=rho_mat.columns)
fig, ax = plt.subplots(figsize=(5.5, 5.0))
sns.heatmap(rho_mat, cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
            annot=annot_mat, fmt="", cbar_kws={"label": "Spearman ρ"}, ax=ax,
            annot_kws={"size": 9})
ax.set_title(f"TAN/NAN composition × malignant MP score\n(per-patient, n={len(combined)} samples)")
ax.set_xlabel(""); ax.set_ylabel("scANVI label")
fig.tight_layout()
fig.savefig(f"{FIG}/step26_tan_mp_correlation_matrix.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step26_tan_mp_correlation_matrix.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ----- scatter: each TAN/NAN vs MP3 -----
print("\n[6] scatter: TAN-1/3/4, NAN-1 × MP3")
fig, axes = plt.subplots(2, 4, figsize=(20, 9))
for ax, lbl in zip(axes.ravel(), SCANVI_LABELS):
    x = combined[f"frac_{lbl}"].values
    y = combined["MP3"].values
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=combined["n_neutrophils"] / 10, alpha=0.7,
               c=combined["n_malignant"], cmap="viridis", edgecolor="black", linewidth=0.4)
    if len(x) > 5:
        z = np.polyfit(x, y, 1)
        xr = np.linspace(x.min(), x.max(), 50)
        ax.plot(xr, np.polyval(z, xr), color="red", lw=1.5, alpha=0.8)
    ax.set_xlabel(f"fraction of {lbl} (per patient)")
    ax.set_ylabel("mean malignant MP3 score")
    ax.set_title(f"{lbl} ↔ MP3   ρ={rho:.2f}  p={p:.2g}")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Per-patient TAN/NAN composition vs malignant MP3 (EMT/IFN program)", y=1.02, fontsize=13)
fig.tight_layout()
fig.savefig(f"{FIG}/step26_scatter_tan_mp3.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step26_scatter_tan_mp3.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ----- NAN-1 (MMP9 specialist) vs all MPs -----
print("\n[7] scatter: NAN-1 (MMP9 axis) × MP1-5")
fig, axes = plt.subplots(1, 5, figsize=(22, 4.5))
for ax, mp in zip(axes, mp_cols):
    x = combined["frac_NAN-1"].values
    y = combined[mp].values
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=combined["n_neutrophils"] / 10, alpha=0.7,
               c=combined["n_malignant"], cmap="viridis", edgecolor="black", linewidth=0.4)
    if len(x) > 5:
        z = np.polyfit(x, y, 1)
        xr = np.linspace(x.min(), x.max(), 50)
        ax.plot(xr, np.polyval(z, xr), color="red", lw=1.5, alpha=0.8)
    ax.set_xlabel("frac NAN-1")
    ax.set_ylabel(f"mean malignant {mp}")
    ax.set_title(f"NAN-1 ↔ {mp}   ρ={rho:.2f}  p={p:.2g}")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("NAN-1 (MMP9 specialist) per-patient × all malignant MPs", y=1.02, fontsize=12)
fig.tight_layout()
fig.savefig(f"{FIG}/step26_scatter_nan1_axes.pdf", bbox_inches="tight")
fig.savefig(f"{FIG}/step26_scatter_nan1_axes.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ----- Bonus: Fig 5B canonical marker dot plot -----
print("\n[8] Fig 5B candidate — canonical TAN/NAN marker dot plot (own data, full 9698 genes)")
raw = sc.read_h5ad(f"{DATA}/luad_neutrophil_own_raw.h5ad")
raw.obs["scanvi_predicted"] = neu.obs.loc[raw.obs.index, "scanvi_predicted"].values
sc.pp.normalize_total(raw, target_sum=1e4); sc.pp.log1p(raw)
SCANVI_ORDER = [l for l in SCANVI_LABELS if l in raw.obs["scanvi_predicted"].unique()]
raw.obs["scanvi_predicted"] = pd.Categorical(raw.obs["scanvi_predicted"],
                                              categories=SCANVI_ORDER, ordered=True)

# Canonical markers (Salcher 2022 + Wu 2023 + Zilionis 2019), filtered to what's available
CANON = {
    "TAN-1": ["VEGFA", "ICAM1", "BCL2A1", "CCRL2"],
    "TAN-2": ["IFIT1", "IFIT2", "ISG15", "MX1", "RSAD2", "OAS1", "IFI6"],
    "TAN-3": ["CXCR4", "CXCL8", "PTGS2", "G0S2"],
    "TAN-4": ["S100A8", "S100A9", "LCN2"],
    "NAN-1": ["CSF3R", "FPR1", "FCGR2A", "IFITM2"],
    "NAN-2": ["VCAN", "CD55", "MNDA"],
    "NAN-3": ["LCN2", "MMP9", "RETN"],
    "core":  ["CXCL8", "IL1B", "MNDA", "CXCR4"],  # neutrophil identity
}
gene_order = []
for k in ["TAN-1", "TAN-2", "TAN-3", "TAN-4", "NAN-1", "NAN-2", "NAN-3", "core"]:
    for g in CANON[k]:
        if g in raw.var_names and g not in gene_order:
            gene_order.append(g)

ax = sc.pl.dotplot(
    raw, var_names=gene_order, groupby="scanvi_predicted",
    standard_scale="var", dendrogram=False,
    categories_order=SCANVI_ORDER, return_fig=True,
    cmap="Reds", figsize=(0.32 * len(gene_order), 0.45 * len(SCANVI_ORDER)),
)
ax.savefig(f"{FIG}/step26_canonical_marker_dotplot.pdf", bbox_inches="tight")
ax.savefig(f"{FIG}/step26_canonical_marker_dotplot.png", dpi=150, bbox_inches="tight")
plt.close()

# ----- TAN/NAN signatures (top-50) for downstream step28 -----
print("\n[9] export TAN/NAN signatures for TCGA ssGSEA (step28)")
m = pd.read_csv(f"{RES}/step25e_markers_scanvi_fullgene.csv")
sig_rows = []
for grp, sub in m.groupby("group"):
    sub = sub[(sub["padj"] < 0.05) & (sub["logfc"] > 0.5)].sort_values("padj").head(50)
    for r in sub.itertuples():
        sig_rows.append({"group": grp, "gene": r.gene, "logfc": r.logfc, "padj": r.padj})
sig = pd.DataFrame(sig_rows)
sig.to_csv(f"{RES}/step26_tan_signatures.csv", index=False)
print(sig.groupby("group").size().to_string())

print(f"\nelapsed: {time.time()-t0:.1f}s")
print("DONE.")
