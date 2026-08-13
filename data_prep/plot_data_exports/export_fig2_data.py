"""
Export data for Figure 2 (A-F) and Figure S2 (cNMF QC).

Fig 2:
  2A: GEP hierarchical clustering heatmap (Spearman correlation matrix)
  2B: UMAP × 4 MP score heatmaps
  2C: UMAP colored by dominant MP
  2D: MP distribution across tissue types
  2E: Dot plot of MP marker genes
  2F: Hallmark GSEA scores across MPs (heatmap)

Fig S2:
  S2A: K selection error curves
  S2B: K distribution histogram + stability
  S2C: GEP proportion across samples
  S2D: Patient mixing entropy (GEP level)
  S2E: GEP linkage matrix + dendrogram order
  S2F: GEP vs MP entropy comparison

Outputs -> ${WORK_ROOT}/luad_figures/fig2/ and fig_s2/
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"

import glob
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

BASE = Path("${PROJECT_ROOT}")
FIG2_OUT = Path("${WORK_ROOT}/luad_figures/fig2")
FIGS2_OUT = Path("${WORK_ROOT}/luad_figures/fig_s2")
FIG2_OUT.mkdir(parents=True, exist_ok=True)
FIGS2_OUT.mkdir(parents=True, exist_ok=True)


def log(msg):
    print(f"[export] {msg}", flush=True)


# 1. Load existing results
log("Loading existing results...")

gep_mp = pd.read_csv(BASE / "results/step6_gep_mp_assignment.csv")
mp_sigs = pd.read_csv(BASE / "results/step6_mp_signatures_top100.csv")
mp_entropy = pd.read_csv(BASE / "results/step6_mp_patient_mixing.csv")
gsea = pd.read_csv(BASE / "results/step7_mp_hallmark_gsea.csv")
cell_scores = pd.read_csv(BASE / "results/step7_mp_cell_scores.csv")
ds_comp = pd.read_csv(BASE / "results/step7_mp_dataset_composition.csv")
k_summary = pd.read_csv(BASE / "results/step5d_k_selection_summary.csv")

log(f"  gep_mp: {gep_mp.shape}, cell_scores: {cell_scores.shape}")
log(f"  gsea: {gsea.shape}, k_summary: {k_summary.shape}")
log(f"  mp_sigs: {mp_sigs.shape} (long format: {mp_sigs.columns.tolist()})")

# 2A / S2E: GEP Spearman correlation + linkage
log("Loading GEP correlation matrix (precomputed by step6) ...")
corr_path_pre = BASE / "results/step6_gep_correlation.csv"
if corr_path_pre.exists():
    corr_df = pd.read_csv(corr_path_pre, index_col=0)
    log(f"  GEP corr (precomputed): {corr_df.shape}")
else:
    log("  step6_gep_correlation.csv missing — recomputing from pool ...")
    pool = pd.read_csv(BASE / "data/cnmf_output/gep_pool_zscore.csv", index_col=0)
    from scipy.stats import rankdata
    ranks = np.empty_like(pool.values, dtype=np.float32)
    for j in range(pool.shape[1]):
        ranks[:, j] = rankdata(pool.values[:, j]).astype(np.float32)
    corr = np.corrcoef(ranks, rowvar=False).astype(np.float32)
    corr_df = pd.DataFrame(corr, index=pool.columns, columns=pool.columns)

# Linkage
dist = (1.0 - corr_df.values).astype(np.float64)
dist = (dist + dist.T) / 2
np.fill_diagonal(dist, 0.0)
dist = np.clip(dist, 0.0, 2.0)
condensed = squareform(dist, checks=False)
Z = linkage(condensed, method="average")
leaf_order = leaves_list(Z)
gep_ids = corr_df.index.to_numpy()

# Save
corr_df.to_csv(FIG2_OUT / "gep_spearman_corr.csv", float_format="%.4f")
np.savetxt(FIG2_OUT / "gep_linkage.csv", Z, delimiter=",",
           header="idx1,idx2,distance,n_members", comments="")
pd.DataFrame({"gep": gep_ids[leaf_order], "order": range(len(leaf_order))}).to_csv(
    FIG2_OUT / "gep_dendrogram_order.csv", index=False
)
log(f"  Saved gep_spearman_corr.csv ({corr_df.shape})")
log(f"  Saved gep_linkage.csv + gep_dendrogram_order.csv")

# GEP -> MP annotation (colour bar)
gep_mp.to_csv(FIG2_OUT / "gep_mp_annotation.csv", index=False)
log(f"  Saved gep_mp_annotation.csv ({len(gep_mp)} GEPs)")

# 2B/2C: Malignant UMAP + MP scores
log("Loading malignant h5ad for UMAP + scores...")

adata_mal = ad.read_h5ad(BASE / "data/processed/luad_malignant_scored.h5ad")
log(f"  Malignant adata: {adata_mal.shape}")
log(f"  obsm keys: {list(adata_mal.obsm.keys())}")

# tissue_type from full adata (load in backed mode, extract, close)
log("Loading copykat h5ad for tissue_type (backed) ...")
adata_full = ad.read_h5ad(BASE / "data/processed/luad_copykat.h5ad", backed="r")
tissue_map = adata_full.obs[["tissue_type"]].copy()
adata_full.file.close()

umap_df = pd.DataFrame(
    adata_mal.obsm["X_umap_mal"],
    index=adata_mal.obs_names,
    columns=["UMAP1", "UMAP2"],
)

# MP score columns (strict pattern, avoid celltype_marker_score)
import re
mp_score_cols = sorted(
    [c for c in adata_mal.obs.columns if re.fullmatch(r"MP\d+_score", c)],
    key=lambda c: int(c.replace("MP", "").replace("_score", "")),
)
log(f"  MP score columns: {mp_score_cols}")

for c in mp_score_cols + ["dominant_MP", "dominant_MP_score"]:
    if c in adata_mal.obs.columns:
        umap_df[c] = adata_mal.obs[c].values

for c in ["dataset", "patient_id", "sample_id"]:
    if c in adata_mal.obs.columns:
        umap_df[c] = adata_mal.obs[c].values

umap_df = umap_df.join(tissue_map, how="left")

umap_df.to_csv(FIG2_OUT / "malignant_umap_metadata.csv.gz", compression="gzip")
log(f"  Saved malignant_umap_metadata.csv.gz ({umap_df.shape})")

# 2D: MP proportion by tissue / dataset / sample
log("Computing MP proportion tables ...")

if "tissue_type" in umap_df.columns:
    ct = pd.crosstab(umap_df["tissue_type"], umap_df["dominant_MP"], normalize="index")
    ct_cnt = pd.crosstab(umap_df["tissue_type"], umap_df["dominant_MP"])
    ct.to_csv(FIG2_OUT / "mp_proportion_by_tissue.csv")
    ct_cnt.to_csv(FIG2_OUT / "mp_count_by_tissue.csv")
    log(f"  Saved mp_proportion_by_tissue.csv + counts ({ct.shape})")

ct2 = pd.crosstab(umap_df["dataset"], umap_df["dominant_MP"], normalize="index")
ct2.to_csv(FIG2_OUT / "mp_proportion_by_dataset.csv")
log(f"  Saved mp_proportion_by_dataset.csv")

if "sample_id" in umap_df.columns:
    ct3 = pd.crosstab(umap_df["sample_id"], umap_df["dominant_MP"], normalize="index")
    ct3_cnt = pd.crosstab(umap_df["sample_id"], umap_df["dominant_MP"])
    ct3.to_csv(FIG2_OUT / "mp_proportion_by_sample.csv")
    ct3_cnt.to_csv(FIG2_OUT / "mp_count_by_sample.csv")
    log(f"  Saved mp_proportion_by_sample.csv + counts ({ct3.shape})")

# 2E: Dot plot data (top-10 marker genes per MP)
log("Computing dot plot data for MP markers (top 10 per MP) ...")

TOP_N = 10
marker_genes = {}
for mp in sorted(mp_sigs["MP"].unique(), key=lambda m: int(m.replace("MP", ""))):
    sub = mp_sigs[mp_sigs["MP"] == mp].sort_values("rank")
    marker_genes[mp] = sub["gene"].head(TOP_N).tolist()

# Unique, preserve order
all_markers = []
for genes in marker_genes.values():
    for g in genes:
        if g not in all_markers:
            all_markers.append(g)
log(f"  {len(all_markers)} unique marker genes across MPs")

# adata_mal.X is already log-normalised (step 7 sets it before saving)
gene_mask = adata_mal.var_names.isin(all_markers)
expr = adata_mal[:, gene_mask]
X = expr.X
if sp.issparse(X):
    X_dense = X.toarray()
else:
    X_dense = np.asarray(X)
gene_names_found = expr.var_names.tolist()
missing = sorted(set(all_markers) - set(gene_names_found))
if missing:
    log(f"  WARN: {len(missing)} markers not found in var_names: {missing[:5]}")

dominant = adata_mal.obs["dominant_MP"].astype(str).values
rows = []
for mp_label in sorted(set(dominant), key=lambda m: int(m.replace("MP", ""))):
    mask = dominant == mp_label
    X_sub = X_dense[mask]
    n_cells = int(mask.sum())
    for i, gene in enumerate(gene_names_found):
        col = X_sub[:, i]
        rows.append({
            "MP": mp_label,
            "gene": gene,
            "frac_expressing": float((col > 0).mean()),
            "mean_expression": float(col.mean()),
            "mean_expr_pos": float(col[col > 0].mean()) if (col > 0).any() else 0.0,
            "n_cells": n_cells,
        })

dotplot_df = pd.DataFrame(rows)

gene_group = {}
for mp, genes in marker_genes.items():
    for g in genes:
        gene_group.setdefault(g, mp)
dotplot_df["gene_group"] = dotplot_df["gene"].map(gene_group)

dotplot_df.to_csv(FIG2_OUT / "mp_dotplot_markers.csv", index=False)
log(f"  Saved mp_dotplot_markers.csv ({dotplot_df.shape})")

gene_order = pd.DataFrame(
    [{"gene": g, "MP": mp, "rank": i}
     for mp, genes in marker_genes.items()
     for i, g in enumerate(genes, start=1)]
)
gene_order.to_csv(FIG2_OUT / "mp_marker_gene_order.csv", index=False)
log(f"  Saved mp_marker_gene_order.csv")

# 2F: Hallmark GSEA heatmap data
log("Preparing Hallmark GSEA heatmap ...")

gsea.to_csv(FIG2_OUT / "hallmark_gsea_full.csv", index=False)

if {"MP", "Term", "NES"}.issubset(gsea.columns):
    nes_pivot = gsea.pivot_table(index="Term", columns="MP",
                                   values="NES", aggfunc="first")
    nes_pivot.to_csv(FIG2_OUT / "hallmark_nes_heatmap.csv")
    log(f"  Saved hallmark_nes_heatmap.csv ({nes_pivot.shape})")
    if "FDR q-val" in gsea.columns:
        fdr_pivot = gsea.pivot_table(index="Term", columns="MP",
                                      values="FDR q-val", aggfunc="first")
        fdr_pivot.to_csv(FIG2_OUT / "hallmark_fdr_heatmap.csv")
        log(f"  Saved hallmark_fdr_heatmap.csv")
else:
    log(f"  WARN: unexpected GSEA columns {gsea.columns.tolist()}")

# MP signature top-100 (raw)
mp_sigs.to_csv(FIG2_OUT / "mp_signatures_top100.csv", index=False)
log(f"  Saved mp_signatures_top100.csv")

# Fig S2: cNMF QC
log("")
log("=== Figure S2: cNMF QC ===")

# S2B: K summary
k_summary.to_csv(FIGS2_OUT / "k_selection_summary.csv", index=False)
log(f"  Saved k_selection_summary.csv ({k_summary.shape})")

# S2A: per-patient error curves
log("Extracting K selection error curves...")
cnmf_dir = BASE / "data/cnmf_input"
error_rows = []
patients_with_curves = 0

for patient_dir in sorted(cnmf_dir.glob("*/cnmf_out")):
    patient_key = patient_dir.parent.name
    stats_files = list(patient_dir.glob(f"{patient_key}/{patient_key}.k_selection_stats.df.npz"))
    if not stats_files:
        continue
    try:
        with np.load(stats_files[0], allow_pickle=True) as z:
            stats_df = pd.DataFrame(**z)
        for _, r in stats_df.iterrows():
            row = {"patient": patient_key}
            row.update(r.to_dict())
            error_rows.append(row)
        patients_with_curves += 1
    except Exception as e:
        log(f"  WARN: failed {stats_files[0]}: {e}")

if error_rows:
    err_df = pd.DataFrame(error_rows)
    err_df.to_csv(FIGS2_OUT / "k_selection_error_curves.csv", index=False)
    log(f"  Saved k_selection_error_curves.csv "
        f"({patients_with_curves} patients, {len(err_df)} rows)")
else:
    log("  WARN: no error curves extracted")

# S2D: GEP patient mixing entropy (MP-level already exists — rename for clarity)
mp_entropy.to_csv(FIGS2_OUT / "mp_patient_mixing.csv", index=False)
log(f"  Saved mp_patient_mixing.csv ({mp_entropy.shape})")

# GEP-level entropy: each GEP belongs to 1 patient by construction, so the
# "mixing" concept only makes sense at MP level. Instead, we save per-MP
# how many patients contribute + what fraction of each MP's GEPs per patient.
log("Computing per-MP patient distribution table (GEP counts per patient)...")
pat_counts = (
    gep_mp.groupby(["MP", "patient_key"]).size()
    .reset_index(name="n_GEPs")
)
pat_counts["total_in_MP"] = pat_counts.groupby("MP")["n_GEPs"].transform("sum")
pat_counts["frac_of_MP"] = pat_counts["n_GEPs"] / pat_counts["total_in_MP"]
pat_counts.to_csv(FIGS2_OUT / "gep_per_mp_per_patient.csv", index=False)
log(f"  Saved gep_per_mp_per_patient.csv ({pat_counts.shape})")

# S2C: GEP proportion across samples — MP-level per-sample already saved.
# Additional: per-patient MP GEP count
pat_mp_ct = gep_mp.pivot_table(index="patient_key", columns="MP",
                                  values="gep_id", aggfunc="count",
                                  fill_value=0)
pat_mp_ct.to_csv(FIGS2_OUT / "gep_count_per_patient_per_mp.csv")
log(f"  Saved gep_count_per_patient_per_mp.csv ({pat_mp_ct.shape})")

# S2E: linkage + dendrogram + corr matrix (copy from fig2)
for f in ["gep_linkage.csv", "gep_dendrogram_order.csv", "gep_spearman_corr.csv"]:
    shutil.copy(FIG2_OUT / f, FIGS2_OUT / f)
log("  Copied linkage/corr to fig_s2/")

# S2F: GEP vs MP entropy comparison
# GEP-level: each GEP is a single patient -> entropy = 0.
# But for display, we characterise per-MP how concentrated its GEPs are.
# Build a 2-row table: (MP, n_patients, entropy, H_norm) + per-GEP "entropy=0" marker.
if {"MP", "entropy", "H_norm", "n_GEPs"}.issubset(mp_entropy.columns):
    rows = []
    for _, r in mp_entropy.iterrows():
        rows.append({
            "level": "MP",
            "name": r["MP"],
            "entropy": float(r["entropy"]),
            "H_norm": float(r["H_norm"]),
            "n_items": int(r["n_GEPs"]),
            "n_patients": int(r["n_patients"]),
        })
    # For GEPs, entropy is 0 by construction (single-patient); add summary:
    rows.append({
        "level": "GEP",
        "name": "all",
        "entropy": 0.0,
        "H_norm": 0.0,
        "n_items": int(mp_entropy["n_GEPs"].sum()),
        "n_patients": 1,
    })
    pd.DataFrame(rows).to_csv(FIGS2_OUT / "gep_vs_mp_entropy.csv", index=False)
    log("  Saved gep_vs_mp_entropy.csv")

ds_comp.to_csv(FIGS2_OUT / "mp_dataset_composition.csv", index=False)
log(f"  Saved mp_dataset_composition.csv")

# Summary
log("")
log("=== Export complete ===")
log(f"Fig 2 files in: {FIG2_OUT}")
for f in sorted(FIG2_OUT.iterdir()):
    log(f"  {f.name}  ({f.stat().st_size/1024:.0f} KB)")

log(f"")
log(f"Fig S2 files in: {FIGS2_OUT}")
for f in sorted(FIGS2_OUT.iterdir()):
    log(f"  {f.name}  ({f.stat().st_size/1024:.0f} KB)")
