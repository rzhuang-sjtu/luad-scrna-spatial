"""step32: prepare CellChat input h5ads (all + Normal/Tumor/Metastasis).

Build a unified `cellchat_label` per cell:
  Priority:
    Mal_MP1/2/3/4   from luad_malignant_scored.h5ad (45,791 cells, all aneuploid)
    Neu_*           from luad_neutrophil_own_annotated.h5ad (8,549 cells, drop Neu_unclassified)
    Macro_SPP1/...  from fig4 panel_major_type_metadata.csv.gz myeloid_subtype_refined
                    (also includes cDC1/2/LAMP3, pDC, Mono_nonclassical)
    Fibroblast / T_NK / B / Plasma / Endothelial / Mast  from celltype_coarse

Subsample within each tissue group:
  Normal      <= 30000   (Normal_Lung + Adjacent_Normal + Normal_LN)
  Tumor       <= 30000   (Primary_Tumor + Precancerous)
  Metastasis  <= 30000   (LN_Met + Brain_Met + Distant_Met + Pleural_Effusion)
  All         <= 50000

Outputs (~/luad/data/processed/cellchat_input_*.h5ad):
  cellchat_input_all.h5ad
  cellchat_input_Normal.h5ad
  cellchat_input_Tumor.h5ad
  cellchat_input_Metastasis.h5ad
"""
import os, time, gc
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

t0 = time.time()
DATA = Path("${PROJECT_ROOT}/data/processed")

print("[1] load merged_annotated (raw counts)")
adata = sc.read_h5ad(DATA / "luad_merged_annotated.h5ad")
print(f"  {adata.shape}; X integer? {np.allclose(adata.X[:200,:200].toarray() if hasattr(adata.X[:200,:200],'toarray') else adata.X[:200,:200], np.round(adata.X[:200,:200].toarray() if hasattr(adata.X[:200,:200],'toarray') else adata.X[:200,:200]))}")

# Save raw counts to layer, then log-normalize (CellChat expects normalized)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(f"  log-normalized; X range: [{adata.X[:200,:200].min():.3f}, {adata.X[:200,:200].max():.3f}]")

print("\n[2] load auxiliary labels")
mal = ad.read_h5ad(DATA / "luad_malignant_scored.h5ad", backed="r")
neu = ad.read_h5ad(DATA / "luad_neutrophil_own_annotated.h5ad", backed="r")
fig4 = pd.read_csv("${WORK_ROOT}/luad_figures/fig4/panel_major_type_metadata.csv.gz")
print(f"  malignant: {mal.shape}; neutrophil: {neu.shape}; fig4 myeloid: {fig4.shape}")

mal_label = pd.Series(index=mal.obs.index.astype(str),
                       data=("Mal_" + mal.obs["dominant_MP"].astype(str)).values)
mal_label = mal_label[~mal_label.str.endswith("nan")]
print(f"  Mal labels: {mal_label.value_counts().to_dict()}")

neu_label = pd.Series(index=neu.obs.index.astype(str),
                       data=neu.obs["neu_subtype"].astype(str).values)
# Drop Neu_unclassified — too noisy for cell-cell comm
neu_label = neu_label[neu_label != "Neu_unclassified"]
print(f"  Neu labels (excl. unclassified): {neu_label.value_counts().to_dict()}")

mye_label = fig4.set_index("barcode")["myeloid_subtype_refined"]
mye_label = mye_label.dropna()
# drop "Neutrophil" from myeloid (handled by neu_label which has finer subtypes)
mye_label = mye_label[mye_label != "Neutrophil"]
print(f"  Myeloid (refined, excl. Neutrophil): {mye_label.value_counts().to_dict()}")

print("\n[3] compose cellchat_label per cell")
n = adata.n_obs
labels = pd.Series(index=adata.obs.index.astype(str), data=pd.NA, dtype="object")

mal_set = set(mal_label.index)
neu_set = set(neu_label.index)
mye_set = set(mye_label.index)
print(f"  set sizes: mal={len(mal_set)}, neu={len(neu_set)}, mye={len(mye_set)}")

# vectorize composition
idx = labels.index.astype(str)
# Step 1: malignant
mask_mal = idx.isin(mal_set)
labels[mask_mal] = mal_label.reindex(idx[mask_mal]).values
# Step 2: neutrophil (subset of myeloid)
mask_neu = (~mask_mal) & idx.isin(neu_set)
labels[mask_neu] = neu_label.reindex(idx[mask_neu]).values
# Step 3: other myeloid (refined subtypes)
mask_mye = (~mask_mal) & (~mask_neu) & idx.isin(mye_set)
labels[mask_mye] = mye_label.reindex(idx[mask_mye]).values
# Step 4: fall back to celltype_coarse for non-malignant non-myeloid
fallback_types = {"Endothelial","Fibroblast","T_NK","B","Plasma","Mast","Epithelial"}
mask_other = labels.isna()
cc = adata.obs["celltype_coarse"].astype(str).values
fb_assign = np.where(np.isin(cc, list(fallback_types)), cc, np.nan)
labels[mask_other] = fb_assign[mask_other]
# Map non-malignant Epithelial → Epithelial_Normal so CellChat doesn't mix them with Mal_*
labels[labels == "Epithelial"] = "Epithelial_Normal"

n_label = labels.notna().sum()
print(f"  cells with label: {n_label} / {n}  ({100*n_label/n:.1f}%)")
adata.obs["cellchat_label"] = labels.values
print("  cellchat_label distribution (top 25):")
print(adata.obs["cellchat_label"].value_counts().head(25).to_string())

adata = adata[adata.obs["cellchat_label"].notna()].copy()
print(f"\n  after drop unlabeled: {adata.shape}")
gc.collect()

print("\n[4] tissue group assignment")
tissue_groups = {
    "Normal":     ["Normal_Lung", "Adjacent_Normal", "Normal_LN"],
    "Tumor":      ["Primary_Tumor", "Precancerous"],
    "Metastasis": ["LN_Metastasis", "Brain_Metastasis",
                   "Distant_Metastasis", "Pleural_Effusion"],
}
adata.obs["tissue_group"] = "other"
for g, ts in tissue_groups.items():
    adata.obs.loc[adata.obs["tissue_type"].isin(ts), "tissue_group"] = g
print(adata.obs["tissue_group"].value_counts().to_string())

def stratified_subsample(adata, max_n, group_col="cellchat_label", seed=42):
    rng = np.random.default_rng(seed)
    if adata.n_obs <= max_n: return adata.copy()
    counts = adata.obs[group_col].value_counts()
    # Allocate proportionally with min cap of 100
    target_per = (max_n * counts / counts.sum()).round().astype(int)
    target_per = np.maximum(target_per, np.minimum(counts, 100))
    idx_list = []
    for grp, n_target in target_per.items():
        ix = np.where(adata.obs[group_col].values == grp)[0]
        if len(ix) > n_target:
            ix = rng.choice(ix, n_target, replace=False)
        idx_list.append(ix)
    keep_idx = np.sort(np.concatenate(idx_list))
    return adata[keep_idx].copy()

# minimize obs/var for export
SLIM_OBS = ["cellchat_label", "tissue_group", "tissue_type", "dataset",
            "patient_id", "sample_id"]

print("\n[5] subsample + save")
# all
all_sub = stratified_subsample(adata, 50000)
all_sub.obs = all_sub.obs[[c for c in SLIM_OBS if c in all_sub.obs.columns]]
out = DATA / "cellchat_input_all.h5ad"
all_sub.write_h5ad(out, compression="gzip")
print(f"  all → {all_sub.shape} → {out} ({os.path.getsize(out)/1e6:.0f} MB)")
print(f"  all label dist: {all_sub.obs['cellchat_label'].value_counts().to_dict()}")
del all_sub; gc.collect()

# per tissue
for g in ["Normal", "Tumor", "Metastasis"]:
    sub = adata[adata.obs["tissue_group"] == g]
    print(f"\n  {g} pre-sub: {sub.shape}")
    if sub.n_obs == 0:
        continue
    sub = stratified_subsample(sub, 30000)
    sub.obs = sub.obs[[c for c in SLIM_OBS if c in sub.obs.columns]]
    out = DATA / f"cellchat_input_{g}.h5ad"
    sub.write_h5ad(out, compression="gzip")
    print(f"  {g} → {sub.shape} → {out} ({os.path.getsize(out)/1e6:.0f} MB)")
    print(f"    label dist: {sub.obs['cellchat_label'].value_counts().head(20).to_dict()}")
    del sub; gc.collect()

print(f"\nelapsed: {(time.time()-t0)/60:.1f} min")
print("DONE.")
