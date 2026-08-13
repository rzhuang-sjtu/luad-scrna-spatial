"""scANVI label transfer: Salcher TAN/NAN → own 8.5k Neutrophils.

Pipeline:
  1. Load Salcher neutrophil_final.h5ad, re-key var_names to HGNC symbols.
  2. Load own neutrophil subset (luad_neutrophil_own_raw.h5ad).
  3. Restrict to shared HGNC symbols, concatenate.
  4. Build batch_key = 'batch' (Salcher dataset ∪ own dataset, prefixed),
     labels_key = 'transfer_label' ('cell_type_neutro' for Salcher, 'Unknown' for own).
  5. HVG=3000 with batch_key on Salcher cells only (avoid HVG bias from own).
  6. Train scVI (max 200 epochs, early stop), then scANVI (max 25 epochs).
  7. Predict labels + soft probabilities on own cells.
  8. Independent leiden subclustering on own cells (scVI latent restricted to own).
  9. UMAP on own-only latent.
 10. Save: data/processed/luad_neutrophil_own_annotated.h5ad
          + results/step25_neutrophil_transfer_summary.csv
          + obsm copies.

Heavy step — runs ~20-40 min on RTX 3080.
"""
import os, sys, time, gc, json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
import torch

t0 = time.time()
LOG = "${PROJECT_ROOT}/results/step25c_scanvi.log"
RNG = 0
scvi.settings.seed = RNG
np.random.seed(RNG)

SAL = "${DATA_ROOT}/High-resolution/neutrophil_final.h5ad"
OWN = "${PROJECT_ROOT}/data/processed/luad_neutrophil_own_raw.h5ad"
OUT_DIR = "${PROJECT_ROOT}/data/processed"
RES_DIR = "${PROJECT_ROOT}/results"

print(f"[env] scvi={scvi.__version__} scanpy={sc.__version__} torch={torch.__version__} cuda={torch.cuda.is_available()}")

# --- 1. Load + harmonize Salcher ---
print(f"\n[load] Salcher: {SAL}")
sal = sc.read_h5ad(SAL)
print(f"  {sal.shape}")

# raw counts in layers['count'] (singular, per inspection)
if "count" in sal.layers:
    sal.X = sal.layers["count"].copy()
    print("  Salcher .X ← layers['count']")
elif "counts" in sal.layers:
    sal.X = sal.layers["counts"].copy()
xs = sal.X[:200, :200].toarray() if hasattr(sal.X, "toarray") else sal.X[:200, :200]
print(f"  Salcher X integer-ish: {np.allclose(xs, np.round(xs))}  range=[{xs.min():.1f},{xs.max():.1f}]")

# var_names → HGNC symbols
sal.var["ensembl_id"] = sal.var_names.astype(str)
sal.var_names = sal.var["feature_name"].astype(str).values
sal.var_names_make_unique()
print(f"  Salcher var renamed to HGNC; n_var={sal.n_vars}, unique={sal.var_names.is_unique}")

# Restrict Salcher labels to its own column
sal.obs["transfer_label"] = sal.obs["cell_type_neutro"].astype(str)
sal.obs["source"] = "salcher"
sal.obs["batch"] = "salcher_" + sal.obs["dataset"].astype(str)
sal_keep = ["transfer_label", "source", "batch", "dataset", "study", "platform",
            "origin", "tissue", "disease", "sample"]
sal_obs_clean = sal.obs[sal_keep].copy()

# --- 2. Load own ---
print(f"\n[load] own: {OWN}")
own = sc.read_h5ad(OWN)
print(f"  {own.shape}")
own.obs["transfer_label"] = "Unknown"
own.obs["source"] = "own"
own.obs["batch"] = "own_" + own.obs["dataset"].astype(str)
own_keep = ["transfer_label", "source", "batch", "dataset", "patient_id", "sample_id",
            "tissue_type", "tissue_stage", "stage", "chemotherapy"]
own_obs_clean = own.obs[own_keep].copy()

# --- 3. Shared genes + concat ---
shared = sal.var_names.intersection(own.var_names)
print(f"\n[shared genes] {len(shared)}")

# subset both to shared genes (preserve order from own)
shared_ordered = [g for g in own.var_names if g in shared]
sal_sub = sal[:, shared_ordered].copy()
own_sub = own[:, shared_ordered].copy()
del sal, own; gc.collect()

# slim Salcher obs to shared cols only
sal_sub.obs = sal_obs_clean
own_sub.obs = own_obs_clean
# pad to a common obs schema
all_cols = sorted(set(sal_sub.obs.columns) | set(own_sub.obs.columns))
for c in all_cols:
    if c not in sal_sub.obs.columns: sal_sub.obs[c] = pd.NA
    if c not in own_sub.obs.columns: own_sub.obs[c] = pd.NA
sal_sub.obs = sal_sub.obs[all_cols]
own_sub.obs = own_sub.obs[all_cols]

print(f"  Salcher sub: {sal_sub.shape}")
print(f"  own sub:     {own_sub.shape}")

concat = ad.concat([sal_sub, own_sub], axis=0, join="outer", merge="same",
                   index_unique=None, label="source_marker")
print(f"  concat: {concat.shape}")
print(f"  concat.obs.transfer_label counts:\n{concat.obs['transfer_label'].value_counts()}")
print(f"  concat.obs.batch n_levels: {concat.obs['batch'].nunique()}")

# raw counts must live in layers['counts'] for scvi
concat.layers["counts"] = concat.X.copy()

# --- 4. HVG (combined data; no batch_key — seurat_v3+batch hits LOESS singularity) ---
print("\n[HVG] selecting 3000 genes (seurat_v3 on counts, no batch_key)...")
# cast to int32 in counts layer to silence "non-integer" warning
import scipy.sparse as sp
if sp.issparse(concat.layers["counts"]):
    concat.layers["counts"].data = np.rint(concat.layers["counts"].data).astype(np.int32)
else:
    concat.layers["counts"] = np.rint(concat.layers["counts"]).astype(np.int32)
sc.pp.highly_variable_genes(
    concat, n_top_genes=3000, flavor="seurat_v3",
    layer="counts", subset=True,
)
print(f"  after HVG: {concat.shape}")

# --- 5. scVI ---
print("\n[scVI] setup + train")
scvi.model.SCVI.setup_anndata(
    concat, layer="counts", batch_key="batch",
    categorical_covariate_keys=None,
)
scvi_model = scvi.model.SCVI(concat, n_layers=2, n_latent=30, gene_likelihood="nb")
scvi_model.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=15,
    batch_size=512,
    accelerator="gpu" if torch.cuda.is_available() else "cpu",
)
print(f"  scVI trained, history len={len(scvi_model.history.get('elbo_train', []))}")

concat.obsm["X_scVI"] = scvi_model.get_latent_representation()

# --- 6. scANVI ---
print("\n[scANVI] from_scvi_model + train")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, adata=concat,
    labels_key="transfer_label", unlabeled_category="Unknown",
)
scanvi_model.train(
    max_epochs=25,
    n_samples_per_label=100,
    accelerator="gpu" if torch.cuda.is_available() else "cpu",
    batch_size=512,
)
concat.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

# --- 7. Predict labels on own cells ---
own_mask = (concat.obs["source"] == "own").values
preds = scanvi_model.predict(concat[own_mask], soft=False)
probs = scanvi_model.predict(concat[own_mask], soft=True)
concat.obs["scanvi_predicted"] = pd.Series(index=concat.obs.index, dtype=object)
concat.obs.loc[own_mask, "scanvi_predicted"] = np.asarray(preds)
# salcher cells keep their true label
concat.obs.loc[~own_mask, "scanvi_predicted"] = concat.obs.loc[~own_mask, "transfer_label"].values
# uncertainty = 1 - max prob
prob_arr = np.asarray(probs)
concat.obs["scanvi_uncertainty"] = np.nan
concat.obs.loc[own_mask, "scanvi_uncertainty"] = 1.0 - prob_arr.max(axis=1)

print("\n[predict] own-cell label dist:")
print(concat.obs.loc[own_mask, "scanvi_predicted"].value_counts())
print(f"\n  mean uncertainty: {concat.obs.loc[own_mask, 'scanvi_uncertainty'].mean():.3f}")
print(f"  median uncertainty: {concat.obs.loc[own_mask, 'scanvi_uncertainty'].median():.3f}")

# --- 8. Subset to own cells, independent leiden + UMAP on scANVI latent ---
print("\n[own-only] neighbors + leiden + UMAP")
own_ad = concat[own_mask].copy()
sc.pp.neighbors(own_ad, use_rep="X_scANVI", n_neighbors=15)
sc.tl.leiden(own_ad, resolution=0.6, key_added="leiden_0.6", random_state=RNG)
sc.tl.leiden(own_ad, resolution=1.0, key_added="leiden_1.0", random_state=RNG)
sc.tl.umap(own_ad, min_dist=0.3, random_state=RNG)
print(f"  leiden_0.6 clusters: {own_ad.obs['leiden_0.6'].nunique()}")
print(f"  leiden_1.0 clusters: {own_ad.obs['leiden_1.0'].nunique()}")

# also compute UMAP for the joint anndata (Salcher+own) so we can render co-embedding fig
print("[joint] UMAP on combined scANVI latent")
sc.pp.neighbors(concat, use_rep="X_scANVI", n_neighbors=15)
sc.tl.umap(concat, min_dist=0.3, random_state=RNG)

# --- 9. Save ---
own_ad.uns["transfer_summary"] = {
    "salcher_n": int((~own_mask).sum()),
    "own_n": int(own_mask.sum()),
    "n_hvg": int(concat.n_vars),
    "scvi_epochs": int(len(scvi_model.history.get("elbo_train", []))),
    "scanvi_unlabeled_category": "Unknown",
    "n_predicted_labels": int(own_ad.obs["scanvi_predicted"].nunique()),
    "mean_uncertainty": float(own_ad.obs["scanvi_uncertainty"].mean()),
}

OUT_OWN = f"{OUT_DIR}/luad_neutrophil_own_annotated.h5ad"
OUT_JOINT = f"{OUT_DIR}/luad_neutrophil_joint_scanvi.h5ad"

print(f"\n[write] own-only annotated → {OUT_OWN}")
own_ad.write_h5ad(OUT_OWN, compression="gzip")
print(f"  size: {os.path.getsize(OUT_OWN)/1e6:.1f} MB")

print(f"[write] joint (Salcher+own) → {OUT_JOINT}")
# slim down concat before writing (no need for full counts again)
concat_slim = concat.copy()
concat_slim.layers.clear()
concat_slim.write_h5ad(OUT_JOINT, compression="gzip")
print(f"  size: {os.path.getsize(OUT_JOINT)/1e6:.1f} MB")

# also save the scvi/scanvi models for reuse
scvi_model.save("${PROJECT_ROOT}/data/processed/scvi_neutrophil_model", overwrite=True)
scanvi_model.save("${PROJECT_ROOT}/data/processed/scanvi_neutrophil_model", overwrite=True)

# --- 10. Predictions table ---
pred_table = own_ad.obs[["scanvi_predicted", "scanvi_uncertainty",
                          "leiden_0.6", "leiden_1.0",
                          "dataset", "patient_id", "sample_id", "tissue_type"]].copy()
pred_table.to_csv(f"{RES_DIR}/step25c_predictions.csv.gz", compression="gzip")
print(f"\n[write] predictions table → {RES_DIR}/step25c_predictions.csv.gz")

print(f"\nelapsed total: {(time.time()-t0)/60:.1f} min")
