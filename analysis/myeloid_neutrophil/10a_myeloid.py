"""Step 10a: LUAD Myeloid subtyping.

  1. Subset Myeloid from luad_copykat.h5ad (cell_type_coarse=='Myeloid')
  2. HVG(3000) on counts layer, drop MT-/RPS/RPL/MRPS/MRPL; scale
  3. PCA(50) → Harmony(dataset) via harmonypy
  4. rapids-singlecell GPU: neighbors(30) → UMAP → Leiden(1.0 & 1.5)
  5. Score each cluster against subtype panels (Macro-* / Mono / DC / Neut)
  6. Patient-level subtype % ↔ mean MP3 score Spearman
  7. M1/M2 polarization scoring per macrophage subtype
  8. Export h5ad + Fig 4 tables.

Outputs:
  - ~/luad/data/processed/luad_myeloid.h5ad
  - ~/luad/results/step10_myeloid_summary.md
  - ~/luad/results/step10_myeloid_markers.csv
  - ~/luad/results/step10_myeloid_mp3_correlation.csv
  - ${WORK_ROOT}/luad_figures/fig4/{myeloid_umap_metadata.csv.gz,
        myeloid_dotplot_markers.csv, myeloid_proportion_by_tissue.csv,
        myeloid_mp3_association.csv, myeloid_m1m2_scores.csv}
"""

from __future__ import annotations
import os, gc, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

IN = Path.home() / "luad/data/processed/luad_copykat.h5ad"
OUT_H5 = Path.home() / "luad/data/processed/luad_myeloid.h5ad"
MP_SCORES = Path.home() / "luad/results/step7_mp_cell_scores.csv"
RES_DIR = Path.home() / "luad/results"
FIGDIR = Path("${WORK_ROOT}/luad_figures/fig4")
FIGDIR.mkdir(parents=True, exist_ok=True)

SUBTYPE_MARKERS = {
    "Macro_general": ["CD68", "CSF1R", "APOE", "C1QA", "C1QB", "C1QC", "LYZ"],
    "Macro_C1QC":    ["C1QA", "C1QB", "C1QC", "APOE"],
    "Macro_SPP1":    ["SPP1", "MMP9", "VEGFA", "TREM2"],
    "Macro_FCN1":    ["S100A8", "S100A9", "VCAN", "CD14", "LYZ"],
    "Macro_FOLR2":   ["FOLR2", "MRC1", "CD163", "SELENOP"],
    "Macro_MARCO":   ["MARCO", "MSR1", "FABP4", "PPARG"],
    "Macro_prolif":  ["TOP2A", "STMN1", "HMGB2"],
    "Mono_classical":     ["CD14", "S100A8", "S100A9", "VCAN"],
    "Mono_nonclassical":  ["FCGR3A", "CDKN1C", "MTSS1"],
    "cDC1":  ["CLEC9A", "BATF3", "IRF8", "WDFY4"],
    "cDC2":  ["CD1C", "FCER1A", "CLEC10A"],
    "cDC_LAMP3": ["LAMP3", "CCR7", "FSCN1", "CCL22"],
    "pDC":   ["TCF4", "IL3RA", "JCHAIN", "LILRA4", "GZMB"],
    "Neutrophil": ["CSF3R", "S100A8", "S100A9"],  # limited coverage post-HVG
}

M1_MARKERS = ["TNF", "IL1B", "CXCL10", "CXCL9", "CXCL11"]   # from var + marker_counts
M2_MARKERS = ["CD163", "MRC1", "MSR1", "TGFB1"]
TISSUE_ORDER = ["Precancerous", "Adjacent_Normal", "Normal_Lung", "Normal_LN",
                "Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
                "Distant_Metastasis", "Pleural_Effusion"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  full shape={a.shape}")

    my = a[a.obs["celltype_coarse"] == "Myeloid"].copy()
    log(f"  Myeloid subset: {my.shape}")
    del a; gc.collect()

    # X is counts? use layers['counts'] explicitly for normalize
    my.X = my.layers["counts"].astype(np.float32).copy()

    # Remove MT/RPS/RPL/MRPS/MRPL before HVG
    bad_prefix = ("MT-", "MTRNR", "RPS", "RPL", "MRPS", "MRPL", "MT.")
    keep_gene = ~my.var_names.str.startswith(bad_prefix)
    log(f"  keep {int(keep_gene.sum())}/{my.n_vars} genes after dropping MT/RPS/RPL")
    my = my[:, keep_gene].copy()

    log("normalize_total + log1p ...")
    sc.pp.normalize_total(my, target_sum=1e4)
    sc.pp.log1p(my)

    log("HVG 3000 (per-dataset)")
    sc.pp.highly_variable_genes(my, n_top_genes=3000, batch_key="dataset",
                                flavor="seurat")
    my.raw = my.copy()     # keep full expression in raw for marker scoring later
    my_hvg = my[:, my.var["highly_variable"]].copy()
    log(f"  HVG subset: {my_hvg.shape}")

    log("scale + PCA(50)")
    sc.pp.scale(my_hvg, max_value=10)
    sc.tl.pca(my_hvg, n_comps=50, random_state=0)

    log("Harmony(dataset) via harmonypy ...")
    import harmonypy
    ho = harmonypy.run_harmony(my_hvg.obsm["X_pca"],
                                my_hvg.obs[["dataset"]], "dataset",
                                max_iter_harmony=15, verbose=True)
    n_cells = my_hvg.n_obs
    try:
        Z = ho.Z_corr.detach().cpu().numpy()  # torch backend
    except AttributeError:
        Z = np.asarray(ho.Z_corr)
    log(f"  harmony Z shape raw: {Z.shape}")
    if Z.shape[0] != n_cells:
        Z = Z.T
    log(f"  harmony Z shape after align: {Z.shape}")
    my_hvg.obsm["X_pca_harmony"] = Z.astype("float32")

    log("GPU neighbors/UMAP/Leiden via rapids-singlecell ...")
    import rapids_singlecell as rsc
    rsc.pp.neighbors(my_hvg, use_rep="X_pca_harmony", n_neighbors=30, random_state=0)
    rsc.tl.umap(my_hvg, random_state=0)
    rsc.tl.leiden(my_hvg, resolution=1.0, random_state=0, key_added="leiden_1.0")
    rsc.tl.leiden(my_hvg, resolution=1.5, random_state=0, key_added="leiden_1.5")

    # Transfer results back to full-gene my
    my.obsm["X_pca"]         = my_hvg.obsm["X_pca"]
    my.obsm["X_pca_harmony"] = my_hvg.obsm["X_pca_harmony"]
    my.obsm["X_umap"]        = my_hvg.obsm["X_umap"]
    my.obs["leiden_1.0"]     = my_hvg.obs["leiden_1.0"].values
    my.obs["leiden_1.5"]     = my_hvg.obs["leiden_1.5"].values

    del my_hvg; gc.collect()

    log("scoring subtype marker panels on full gene set ...")
    # Use marker_counts obsm to supplement missing genes (CLEC9A, CD1C, etc.)
    mc_df = my.obsm.get("marker_counts", None)
    # For genes not in var_names but in marker_counts, score manually via z-score
    for subtype, genes in SUBTYPE_MARKERS.items():
        present_var = [g for g in genes if g in my.var_names]
        present_mc = [g for g in genes if (g not in my.var_names) and
                      (mc_df is not None) and (g in mc_df.columns)]
        if present_var:
            sc.tl.score_genes(my, gene_list=present_var, score_name=f"sc_{subtype}",
                              random_state=0, use_raw=False)
        else:
            my.obs[f"sc_{subtype}"] = 0.0
        if present_mc:
            # add mean z-score of supplementary markers
            sub = mc_df[present_mc]
            z = (sub - sub.mean(axis=0)) / sub.std(axis=0).replace(0, 1)
            my.obs[f"sc_{subtype}"] = my.obs[f"sc_{subtype}"].values + z.mean(axis=1).values
        log(f"  {subtype}: var={len(present_var)}/{len(genes)} mc_extra={len(present_mc)}")

    # Per-cluster mean score → vote label
    log("cluster → subtype voting (leiden_1.0)")
    score_cols = [c for c in my.obs.columns if c.startswith("sc_")]
    clust = my.obs.groupby("leiden_1.0", observed=True)[score_cols].mean()
    best = clust.idxmax(axis=1).str.replace("sc_", "")
    second = clust.apply(lambda r: r.nlargest(2).index[-1], axis=1).str.replace("sc_", "")
    margin = clust.apply(lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[-1], axis=1)
    vote = pd.DataFrame({
        "cluster": clust.index,
        "best_subtype": best.values,
        "second_subtype": second.values,
        "margin": margin.values,
        "n_cells": my.obs.groupby("leiden_1.0", observed=True).size().reindex(clust.index).values,
    })
    vote.to_csv(RES_DIR / "step10_myeloid_cluster_votes.csv", index=False)
    log(vote.to_string(index=False))

    my.obs["myeloid_subtype"] = my.obs["leiden_1.0"].map(
        dict(zip(vote["cluster"], vote["best_subtype"]))
    ).astype("category")

    # Apply Neutrophil rescue: if cluster has high CSF3R AND low CD68/APOE/C1QA → Neutrophil
    csf3r = my.raw[:, "CSF3R"].X.toarray().ravel() if "CSF3R" in my.raw.var_names else None
    cd68 = my.raw[:, "CD68"].X.toarray().ravel() if "CD68" in my.raw.var_names else None
    c1qa = my.raw[:, "C1QA"].X.toarray().ravel() if "C1QA" in my.raw.var_names else None
    if csf3r is not None and cd68 is not None and c1qa is not None:
        neut_rule_rows = []
        for cl, g in my.obs.groupby("leiden_1.0", observed=True).groups.items():
            idx = my.obs.index.get_indexer(g)
            neut_rule_rows.append({
                "cluster": cl, "n": len(idx),
                "mean_CSF3R": float(csf3r[idx].mean()),
                "mean_CD68":  float(cd68[idx].mean()),
                "mean_C1QA":  float(c1qa[idx].mean()),
            })
        neut_rule = pd.DataFrame(neut_rule_rows)
        neut_rule["neutrophil_likely"] = (
            (neut_rule["mean_CSF3R"] > neut_rule["mean_CSF3R"].median() * 1.5) &
            (neut_rule["mean_CD68"] < neut_rule["mean_CD68"].median() * 0.5) &
            (neut_rule["mean_C1QA"] < neut_rule["mean_C1QA"].median() * 0.5)
        )
        neut_rule.to_csv(RES_DIR / "step10_myeloid_neutrophil_rule.csv", index=False)
        rescued = neut_rule[neut_rule["neutrophil_likely"]]["cluster"].tolist()
        if rescued:
            log(f"  Neutrophil rescue applied to clusters: {rescued}")
            mask = my.obs["leiden_1.0"].isin(rescued)
            my.obs.loc[mask, "myeloid_subtype"] = "Neutrophil"
            # refresh category
            my.obs["myeloid_subtype"] = my.obs["myeloid_subtype"].astype("category")

    log("rank_genes_groups per myeloid_subtype")
    sc.tl.rank_genes_groups(my, groupby="myeloid_subtype",
                             method="wilcoxon", use_raw=True,
                             n_genes=50, pts=True)
    # Export top DEGs
    rg = my.uns["rank_genes_groups"]
    deg_rows = []
    for grp in rg["names"].dtype.names:
        for rank_i in range(min(30, len(rg["names"][grp]))):
            deg_rows.append({
                "subtype": grp,
                "rank": rank_i + 1,
                "gene": rg["names"][grp][rank_i],
                "score": float(rg["scores"][grp][rank_i]),
                "logFC": float(rg["logfoldchanges"][grp][rank_i]),
                "pval_adj": float(rg["pvals_adj"][grp][rank_i]),
            })
    pd.DataFrame(deg_rows).to_csv(RES_DIR / "step10_myeloid_markers.csv", index=False)

    log("proportion tables")
    prop_tis = (pd.crosstab(my.obs["tissue_type"], my.obs["myeloid_subtype"],
                             normalize="index")
                .reindex(index=[t for t in TISSUE_ORDER if t in
                                my.obs["tissue_type"].unique()]))
    prop_tis.to_csv(FIGDIR / "myeloid_proportion_by_tissue.csv")

    prop_ds = pd.crosstab(my.obs["dataset"], my.obs["myeloid_subtype"],
                           normalize="index")
    prop_ds.to_csv(FIGDIR / "myeloid_proportion_by_dataset.csv")

    # Dotplot marker table (subtype × gene → mean log1p + pct expressing)
    dot_genes = sorted(set(g for gs in SUBTYPE_MARKERS.values() for g in gs))
    dot_genes = [g for g in dot_genes if g in my.raw.var_names]
    raw_sub = my.raw[:, dot_genes].X
    dot_rows = []
    for st in my.obs["myeloid_subtype"].unique():
        mask = (my.obs["myeloid_subtype"] == st).values
        if mask.sum() == 0: continue
        vals = raw_sub[mask]
        mean_expr = np.asarray(vals.mean(axis=0)).ravel()
        pct_expr = np.asarray((vals > 0).mean(axis=0)).ravel()
        for i, g in enumerate(dot_genes):
            dot_rows.append({"subtype": st, "gene": g,
                             "mean_log1p": mean_expr[i],
                             "pct_expressing": pct_expr[i]})
    pd.DataFrame(dot_rows).to_csv(FIGDIR / "myeloid_dotplot_markers.csv", index=False)

    log("M1/M2 polarization scoring")
    for name, genes in [("M1", M1_MARKERS), ("M2", M2_MARKERS)]:
        present_var = [g for g in genes if g in my.var_names]
        present_mc = [g for g in genes if (g not in my.var_names) and
                      (mc_df is not None) and (g in mc_df.columns)]
        if present_var:
            sc.tl.score_genes(my, gene_list=present_var, score_name=f"{name}_score",
                              random_state=0, use_raw=False)
        else:
            my.obs[f"{name}_score"] = 0.0
        if present_mc:
            sub = mc_df[present_mc]
            z = (sub - sub.mean(axis=0)) / sub.std(axis=0).replace(0, 1)
            my.obs[f"{name}_score"] = my.obs[f"{name}_score"].values + z.mean(axis=1).values

    my.obs["M1_M2_ratio"] = my.obs["M1_score"] - my.obs["M2_score"]
    m1m2 = (my.obs.groupby("myeloid_subtype", observed=True)[["M1_score","M2_score","M1_M2_ratio"]]
            .agg(["mean","median","std","count"]))
    m1m2.columns = [f"{a}_{b}" for a, b in m1m2.columns]
    m1m2.to_csv(FIGDIR / "myeloid_m1m2_scores.csv")
    log(m1m2.round(3).to_string())

    log("patient-level Myeloid subtype % vs mean MP3")
    mp = pd.read_csv(MP_SCORES, index_col=0)
    mp_patient = mp.groupby(mp["patient_id"] if "patient_id" in mp.columns
                             else mp.index.map(lambda x: x.split("_")[0]))[["MP1_score","MP2_score","MP3_score","MP4_score"]].mean()
    # Above pattern may not work for all barcodes; use the malignant h5ad for patient_id
    # Simpler: re-load malignant obs for patient_id
    mal = sc.read_h5ad(Path.home()/"luad/data/processed/luad_malignant_scored.h5ad",
                       backed="r")
    mal_df = mal.obs[["patient_id","MP1_score","MP2_score","MP3_score","MP4_score"]]
    mp_patient = mal_df.groupby("patient_id", observed=True).mean()

    my_patient = (pd.crosstab(my.obs["patient_id"], my.obs["myeloid_subtype"],
                               normalize="index"))
    joint = my_patient.join(mp_patient, how="inner")
    log(f"  joint patient table: {joint.shape}")

    from scipy.stats import spearmanr
    corr_rows = []
    for subtype in my_patient.columns:
        for mpx in ["MP1_score","MP2_score","MP3_score","MP4_score"]:
            sub = joint[[subtype, mpx]].dropna()
            if len(sub) < 5: continue
            rho, p = spearmanr(sub[subtype], sub[mpx])
            corr_rows.append({"subtype": subtype, "MP": mpx.replace("_score",""),
                              "n_patients": len(sub),
                              "spearman_rho": rho, "p": p})
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(RES_DIR / "step10_myeloid_mp3_correlation.csv", index=False)
    corr.to_csv(FIGDIR / "myeloid_mp3_association.csv", index=False)
    log("MP3-focused associations:")
    log(corr[corr["MP"] == "MP3"].sort_values("spearman_rho", ascending=False).to_string(index=False))

    log("exporting UMAP metadata")
    umap_df = pd.DataFrame({
        "barcode": my.obs.index,
        "dataset": my.obs["dataset"].astype(str).values,
        "patient_id": my.obs["patient_id"].astype(str).values,
        "tissue_type": my.obs["tissue_type"].astype(str).values,
        "leiden_1.0": my.obs["leiden_1.0"].astype(str).values,
        "leiden_1.5": my.obs["leiden_1.5"].astype(str).values,
        "myeloid_subtype": my.obs["myeloid_subtype"].astype(str).values,
        "UMAP1": my.obsm["X_umap"][:, 0],
        "UMAP2": my.obsm["X_umap"][:, 1],
        "M1_score": my.obs["M1_score"].values,
        "M2_score": my.obs["M2_score"].values,
    })
    umap_df.to_csv(FIGDIR / "myeloid_umap_metadata.csv.gz",
                    index=False, compression="gzip")

    # Save h5ad — drop raw to keep file size manageable
    log(f"saving {OUT_H5}")
    my_out = my.copy()
    del my_out.raw  # we keep full-gene X; raw is redundant
    my_out.write(OUT_H5)

    # Summary md
    log("writing summary md")
    with open(RES_DIR / "step10_myeloid_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 10a — LUAD Myeloid subtyping\n\n")
        f.write(f"- Total Myeloid cells: {my.n_obs}\n")
        f.write(f"- Clusters (Leiden 1.0): {my.obs['leiden_1.0'].nunique()}\n")
        f.write(f"- Distinct subtypes assigned: {my.obs['myeloid_subtype'].nunique()}\n\n")
        f.write("## Subtype cell counts\n\n")
        f.write(my.obs["myeloid_subtype"].value_counts().to_frame("n_cells").to_markdown() + "\n\n")
        f.write("## Cluster → subtype vote\n\n")
        f.write(vote.round(3).to_markdown(index=False) + "\n\n")
        f.write("## M1/M2 polarization per subtype\n\n")
        f.write(m1m2.round(3).to_markdown() + "\n\n")
        f.write("## Top-20 MP3 ↔ Myeloid subtype % Spearman (patient-level)\n\n")
        top = (corr[corr["MP"] == "MP3"]
               .sort_values("spearman_rho", ascending=False).head(20))
        f.write(top.round(4).to_markdown(index=False) + "\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
