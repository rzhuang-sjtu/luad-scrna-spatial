"""Step 11: LUAD T_NK subclustering (423k cells).

Pipeline follows Step 10 Myeloid framework with adjustments:
  - Vote at Leiden 1.5 directly (finer res → higher chance of resolving minor subtypes).
  - Stringent Treg vs CD8-Exhausted disambiguation rescue rules.
  - Auto-detect ambiguous clusters (margin < 0.05) and flag them.
  - Marker panels pared to what survives after 9881-HVG filtering;
    each panel supplemented from obsm['marker_counts'] where possible.

Outputs:
  - ~/luad/data/processed/luad_tnk.h5ad
  - ~/luad/results/step11_tnk_summary.md
  - ~/luad/results/step11_tnk_cluster_votes.csv
  - ~/luad/results/step11_tnk_markers.csv
  - ~/luad/results/step11_tnk_mp3_correlation.csv
  - ~/luad/results/step11_marker_matrix.csv
  - ~/luad/results/step11_diagnosis.md
  - ${WORK_ROOT}/luad_figures/fig_tnk/:
      tnk_umap_metadata.csv.gz, tnk_dotplot_markers.csv,
      tnk_proportion_by_tissue.csv, tnk_mp_association.csv,
      tnk_exhausted_treg_profile.csv
"""

from __future__ import annotations
import os, gc, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

IN = Path.home() / "luad/data/processed/luad_copykat.h5ad"
OUT_H5 = Path.home() / "luad/data/processed/luad_tnk.h5ad"
MAL = Path.home() / "luad/data/processed/luad_malignant_scored.h5ad"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig_tnk")
FIG.mkdir(parents=True, exist_ok=True)

# Marker panels (only use genes likely present in var or obsm['marker_counts'])
PANELS = {
    "CD8_lineage":       ["CD8A", "CD8B"],
    "CD4_lineage":       ["CD4", "IL7R"],
    "T_general":         ["CD3D", "CD3E", "CD3G"],
    "Naive_like":        ["SELL", "TCF7"],                     # weaker without CCR7
    "Effector_CTL":      ["GZMB", "GZMA", "PRF1", "NKG7", "GNLY", "CCL5"],
    "Exhausted_CD8":     ["PDCD1", "HAVCR2", "LAG3", "TIGIT"],  # + marker_counts has HAVCR2, LAG3
    "TRM":               ["ITGAE", "CD69", "RUNX3"],
    "Treg":              ["FOXP3", "IL2RA", "CTLA4"],          # + IKZF2, TNFRSF4
    "NK":                ["FCGR3A", "KLRD1", "NKG7", "GNLY"],
    "NKT":               ["ZBTB16"],
    "gdT":               ["TRGC2"],
    "Prolif":            ["TOP2A", "STMN1", "HMGB2"],
}

# Genes used for rescue-rule expression lookups
RESCUE_MARKERS = ["CD8A", "CD8B", "CD4", "FOXP3", "IL2RA", "CTLA4",
                  "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOP2A", "MKI67",
                  "TRGC2", "ZBTB16", "FCGR3A", "NKG7", "GNLY", "ITGAE",
                  "CD69", "RUNX3", "SELL", "TCF7", "GZMB", "GZMK", "KLRB1"]

LEIDEN_RES = 1.5    # vote at finer res first


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def score_panels(a, panels, mc_df=None):
    """Score each marker panel on adata.obs as sc_<panel>. Returns dict of (n_var, n_mc)."""
    info = {}
    for panel_name, genes in panels.items():
        pv = [g for g in genes if g in a.var_names]
        pm = [g for g in genes if (g not in a.var_names) and
              (mc_df is not None) and (g in mc_df.columns)]
        if pv:
            sc.tl.score_genes(a, gene_list=pv, score_name=f"sc_{panel_name}",
                              random_state=0, use_raw=False)
        else:
            a.obs[f"sc_{panel_name}"] = 0.0
        if pm:
            sub = mc_df[pm]
            z = (sub - sub.mean(axis=0)) / sub.std(axis=0).replace(0, 1)
            a.obs[f"sc_{panel_name}"] = a.obs[f"sc_{panel_name}"].values + z.mean(axis=1).values
        info[panel_name] = (len(pv), len(pm))
    return info


def vote_at_resolution(a, res_col, score_cols):
    """Mean score per cluster → best/2nd subtype + margin."""
    clust = a.obs.groupby(res_col, observed=True)[score_cols].mean()
    best = clust.idxmax(axis=1).str.replace("sc_", "")
    second = clust.apply(lambda r: r.nlargest(2).index[-1], axis=1).str.replace("sc_", "")
    margin = clust.apply(lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[-1], axis=1)
    vote = pd.DataFrame({
        "cluster": clust.index.astype(str),
        "best_panel": best.values,
        "second_panel": second.values,
        "margin": margin.values,
        "n_cells": a.obs.groupby(res_col, observed=True).size().reindex(clust.index).values,
    })
    return vote, clust


def main() -> None:
    t0 = time.time()
    log(f"loading {IN}")
    a = sc.read_h5ad(IN)
    log(f"  full shape={a.shape}")

    tnk = a[a.obs["celltype_coarse"] == "T_NK"].copy()
    log(f"  T_NK subset: {tnk.shape}")
    del a; gc.collect()

    tnk.X = tnk.layers["counts"].astype(np.float32).copy()

    bad_prefix = ("MT-", "MTRNR", "RPS", "RPL", "MRPS", "MRPL", "MT.")
    keep = ~tnk.var_names.str.startswith(bad_prefix)
    log(f"  keep {int(keep.sum())}/{tnk.n_vars} genes after dropping MT/RPS/RPL")
    tnk = tnk[:, keep].copy()

    log("normalize_total + log1p")
    # Reset any pre-existing log1p uns flag so scanpy applies fresh log1p
    if "log1p" in tnk.uns: del tnk.uns["log1p"]
    sc.pp.normalize_total(tnk, target_sum=1e4)
    sc.pp.log1p(tnk)

    log("HVG 3000 (per-dataset)")
    sc.pp.highly_variable_genes(tnk, n_top_genes=3000, batch_key="dataset",
                                flavor="seurat")
    tnk_hvg = tnk[:, tnk.var["highly_variable"]].copy()
    log(f"  HVG subset: {tnk_hvg.shape}")

    log("scale + PCA(50)")
    sc.pp.scale(tnk_hvg, max_value=10)
    sc.tl.pca(tnk_hvg, n_comps=50, random_state=0)

    log("Harmony(dataset) via harmonypy")
    import harmonypy
    ho = harmonypy.run_harmony(tnk_hvg.obsm["X_pca"],
                                tnk_hvg.obs[["dataset"]], "dataset",
                                max_iter_harmony=15, verbose=True)
    n_cells = tnk_hvg.n_obs
    try:
        Z = ho.Z_corr.detach().cpu().numpy()
    except AttributeError:
        Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != n_cells:
        Z = Z.T
    tnk_hvg.obsm["X_pca_harmony"] = Z.astype("float32")
    log(f"  Harmony done, Z shape: {Z.shape}")

    log("GPU neighbors/UMAP/Leiden")
    import rapids_singlecell as rsc
    rsc.pp.neighbors(tnk_hvg, use_rep="X_pca_harmony", n_neighbors=30,
                     random_state=0)
    rsc.tl.umap(tnk_hvg, random_state=0)
    rsc.tl.leiden(tnk_hvg, resolution=1.0, random_state=0, key_added="leiden_1.0")
    rsc.tl.leiden(tnk_hvg, resolution=1.5, random_state=0, key_added="leiden_1.5")

    # Transfer back
    for k in ["X_pca", "X_pca_harmony", "X_umap"]:
        tnk.obsm[k] = tnk_hvg.obsm[k]
    tnk.obs["leiden_1.0"] = tnk_hvg.obs["leiden_1.0"].values
    tnk.obs["leiden_1.5"] = tnk_hvg.obs["leiden_1.5"].values
    del tnk_hvg; gc.collect()

    log("scoring subtype panels")
    mc_df = tnk.obsm.get("marker_counts", None)
    info = score_panels(tnk, PANELS, mc_df=mc_df)
    for k, (v, m) in info.items():
        log(f"  {k:18s} var={v} + mc={m}")

    # Build cluster × rescue-marker mean matrix for downstream rules
    log("computing per-cluster marker matrix")
    present_rescue = [g for g in RESCUE_MARKERS if g in tnk.var_names]
    log(f"  rescue markers present in var: {len(present_rescue)}/{len(RESCUE_MARKERS)}")
    X_rescue = tnk[:, present_rescue].X
    df_rescue = pd.DataFrame(
        X_rescue.toarray() if hasattr(X_rescue, "toarray") else np.asarray(X_rescue),
        index=tnk.obs.index, columns=present_rescue,
    )
    df_rescue["cluster_1.5"] = tnk.obs["leiden_1.5"].astype(str).values
    marker_mean = df_rescue.groupby("cluster_1.5").mean()
    marker_mean.to_csv(RES / "step11_marker_matrix.csv")

    log("vote at Leiden 1.5")
    score_cols = [c for c in tnk.obs.columns if c.startswith("sc_")]
    vote, clust = vote_at_resolution(tnk, "leiden_1.5", score_cols)
    # attach marker-hint columns for interpretation
    hint = ["CD8A", "CD8B", "CD4", "FOXP3", "PDCD1", "HAVCR2", "LAG3", "TOP2A",
            "FCGR3A", "NKG7", "TRGC2", "ITGAE", "GZMB"]
    hint = [h for h in hint if h in marker_mean.columns]
    vote = vote.merge(marker_mean[hint].round(3), left_on="cluster",
                      right_index=True, how="left")
    vote.to_csv(RES / "step11_tnk_cluster_votes.csv", index=False)

    # Simple mapping from panel-name best → canonical subtype
    PANEL2STYPE = {
        "CD8_lineage":     "CD8_generic",
        "CD4_lineage":     "CD4_generic",
        "T_general":       "CD4_generic",
        "Naive_like":      "T_Naive_like",
        "Effector_CTL":    "CD8_Effector",
        "Exhausted_CD8":   "CD8_Exhausted",
        "TRM":             "CD8_TRM",
        "Treg":            "Treg",
        "NK":              "NK",
        "NKT":             "NKT",
        "gdT":             "gdT",
        "Prolif":          "T_Proliferating",
    }
    initial = vote.set_index("cluster")["best_panel"].map(PANEL2STYPE).fillna("Unknown")
    mapping = initial.to_dict()

    med = marker_mean.median()
    FOXP3_TH = max(0.6, marker_mean["FOXP3"].quantile(0.80)) if "FOXP3" in marker_mean else 0.6
    IL2RA_TH = marker_mean["IL2RA"].quantile(0.70) if "IL2RA" in marker_mean else 0.3
    CD8_TH   = max(1.0, marker_mean["CD8A"].quantile(0.50)) if "CD8A" in marker_mean else 1.0
    CD4_TH   = max(0.8, marker_mean["CD4"].quantile(0.50)) if "CD4" in marker_mean else 0.8
    EXH_TH   = marker_mean["PDCD1"].quantile(0.75) if "PDCD1" in marker_mean else 0.3
    TRGC2_TH = marker_mean["TRGC2"].quantile(0.90) if "TRGC2" in marker_mean else 0.5
    TOP2A_TH = marker_mean["TOP2A"].quantile(0.90) if "TOP2A" in marker_mean else 0.5

    log(f"  thresholds: FOXP3={FOXP3_TH:.2f} IL2RA={IL2RA_TH:.2f} "
        f"CD8={CD8_TH:.2f} CD4={CD4_TH:.2f} PDCD1={EXH_TH:.2f} "
        f"TRGC2={TRGC2_TH:.2f} TOP2A={TOP2A_TH:.2f}")

    # Rule priority:
    #   1. Proliferating (if TOP2A top in cluster AND TOP2A > threshold)
    #   2. Treg: FOXP3 >= FOXP3_TH AND IL2RA >= IL2RA_TH AND CD4 >= CD4_TH
    #   3. CD8_Exhausted: PDCD1+ AND CD8A/B dominant AND FOXP3 < FOXP3_TH
    #   4. gdT: TRGC2 >> median (top decile)
    #   5. NK: CD3-low AND NKG7/GNLY high (approximated via best panel NK)
    #   6. CD8 vs CD4 lineage refinement for "generic" labels

    decisions = {}
    for cl in vote["cluster"]:
        row = marker_mean.loc[cl]
        decision = mapping[cl]

        # 1. Proliferating
        if (row.get("TOP2A", 0) >= TOP2A_TH and
            decision not in ("NK", "NKT")):
            decision = "T_Proliferating"

        # 2. Treg (strict co-expression FOXP3 + IL2RA + CD4)
        if (row.get("FOXP3", 0) >= FOXP3_TH and
            row.get("IL2RA", 0) >= IL2RA_TH and
            row.get("CD4", 0)   >= CD4_TH):
            decision = "Treg"

        # 3. CD8_Exhausted — PDCD1 high, CD8+, not Treg
        elif (row.get("PDCD1", 0) >= EXH_TH and
              max(row.get("CD8A", 0), row.get("CD8B", 0)) >= CD8_TH and
              row.get("FOXP3", 0) < FOXP3_TH):
            decision = "CD8_Exhausted"

        # 4. gdT
        elif row.get("TRGC2", 0) >= TRGC2_TH:
            decision = "gdT"

        # 5. CD8 vs CD4 disambiguation for generic
        if decision in ("CD8_generic", "CD4_generic", "T_Naive_like"):
            cd8max = max(row.get("CD8A", 0), row.get("CD8B", 0))
            cd4exp = row.get("CD4", 0)
            if cd8max > cd4exp + 0.2:   # CD8 dominant
                if decision == "T_Naive_like":
                    decision = "CD8_Naive_CM"
                else:
                    decision = "CD8_generic"
            elif cd4exp > cd8max + 0.2:  # CD4 dominant
                if decision == "T_Naive_like":
                    decision = "CD4_Naive_CM"
                else:
                    decision = "CD4_generic"
            # else keep generic label (ambiguous lineage)

        decisions[cl] = decision

    tnk.obs["tnk_subtype"] = (
        tnk.obs["leiden_1.5"].astype(str).map(decisions).astype("category")
    )
    stype_counts = tnk.obs["tnk_subtype"].value_counts()
    log("tnk_subtype counts:")
    log(stype_counts.to_string())

    ambig = vote[vote["margin"] < 0.05].copy()
    if len(ambig):
        log(f"Ambiguous clusters (margin<0.05): {ambig['cluster'].tolist()}")
    ambig.to_csv(RES / "step11_tnk_ambiguous_clusters.csv", index=False)

    log("rank_genes_groups per tnk_subtype (n_genes=30)")
    sc.tl.rank_genes_groups(tnk, groupby="tnk_subtype",
                             method="wilcoxon", use_raw=False,
                             n_genes=30, pts=True)
    rg = tnk.uns["rank_genes_groups"]
    deg_rows = []
    for grp in rg["names"].dtype.names:
        for i in range(len(rg["names"][grp])):
            deg_rows.append({
                "subtype": grp, "rank": i+1,
                "gene": rg["names"][grp][i],
                "score": float(rg["scores"][grp][i]),
                "logFC": float(rg["logfoldchanges"][grp][i]),
                "pval_adj": float(rg["pvals_adj"][grp][i]),
            })
    pd.DataFrame(deg_rows).to_csv(RES / "step11_tnk_markers.csv", index=False)

    log("proportion tables")
    TISSUE_ORDER = ["Precancerous", "Adjacent_Normal", "Normal_Lung", "Normal_LN",
                    "Primary_Tumor", "LN_Metastasis", "Brain_Metastasis",
                    "Distant_Metastasis", "Pleural_Effusion"]
    prop_tis = (pd.crosstab(tnk.obs["tissue_type"], tnk.obs["tnk_subtype"],
                             normalize="index")
                .reindex(index=[t for t in TISSUE_ORDER if t in
                                tnk.obs["tissue_type"].unique()]))
    prop_tis.to_csv(FIG / "tnk_proportion_by_tissue.csv")

    prop_ds = pd.crosstab(tnk.obs["dataset"], tnk.obs["tnk_subtype"],
                           normalize="index")
    prop_ds.to_csv(FIG / "tnk_proportion_by_dataset.csv")

    # Dotplot
    all_markers = sorted(set(g for gs in PANELS.values() for g in gs))
    all_markers = [g for g in all_markers if g in tnk.var_names]
    X_all = tnk[:, all_markers].X
    X_all = X_all.toarray() if hasattr(X_all, "toarray") else np.asarray(X_all)
    dot_rows = []
    for st in tnk.obs["tnk_subtype"].cat.categories:
        mask = (tnk.obs["tnk_subtype"] == st).values
        if not mask.any(): continue
        sub = X_all[mask]
        mean_expr = sub.mean(axis=0)
        pct_expr = (sub > 0).mean(axis=0)
        for i, g in enumerate(all_markers):
            dot_rows.append({"subtype": st, "gene": g,
                             "mean_log1p": float(mean_expr[i]),
                             "pct_expressing": float(pct_expr[i])})
    pd.DataFrame(dot_rows).to_csv(FIG / "tnk_dotplot_markers.csv", index=False)

    # For each of {Treg, CD8_Exhausted}, export mean expression of checkpoint
    # markers + lineage markers across all rescue_markers
    checkpoint_genes = ["FOXP3", "IL2RA", "CTLA4", "PDCD1", "HAVCR2", "LAG3",
                        "TIGIT", "CD4", "CD8A", "CD8B"]
    cg = [g for g in checkpoint_genes if g in tnk.var_names]
    profile_rows = []
    for st in tnk.obs["tnk_subtype"].cat.categories:
        mask = (tnk.obs["tnk_subtype"] == st).values
        if not mask.any(): continue
        X_cg = tnk[mask, cg].X
        X_cg = X_cg.toarray() if hasattr(X_cg, "toarray") else np.asarray(X_cg)
        m = X_cg.mean(axis=0)
        pct = (X_cg > 0).mean(axis=0)
        for i, g in enumerate(cg):
            profile_rows.append({"subtype": st, "gene": g,
                                 "mean_log1p": float(m[i]),
                                 "pct_expressing": float(pct[i])})
    pd.DataFrame(profile_rows).to_csv(FIG / "tnk_exhausted_treg_profile.csv",
                                       index=False)

    log("patient-level subtype% × MP score Spearman")
    mal = sc.read_h5ad(MAL, backed="r")
    mp_pat = (mal.obs[["patient_id", "MP1_score", "MP2_score",
                        "MP3_score", "MP4_score"]]
                 .groupby("patient_id", observed=True).mean())
    tnk_pat = pd.crosstab(tnk.obs["patient_id"], tnk.obs["tnk_subtype"],
                           normalize="index")
    joint = tnk_pat.join(mp_pat, how="inner")
    log(f"  joint patient table: {joint.shape}")
    rows = []
    for st in tnk_pat.columns:
        for mpx in ["MP1_score", "MP2_score", "MP3_score", "MP4_score"]:
            sub = joint[[st, mpx]].dropna()
            if len(sub) < 5: continue
            rho, p = spearmanr(sub[st], sub[mpx])
            rows.append({"subtype": st, "MP": mpx.replace("_score", ""),
                         "n_patients": len(sub),
                         "spearman_rho": rho, "p": p})
    mp_corr = pd.DataFrame(rows)
    mp_corr.to_csv(RES / "step11_tnk_mp3_correlation.csv", index=False)
    mp_corr.to_csv(FIG / "tnk_mp_association.csv", index=False)

    log("exporting UMAP metadata")
    umap_df = pd.DataFrame({
        "barcode": tnk.obs.index,
        "dataset": tnk.obs["dataset"].astype(str).values,
        "patient_id": tnk.obs["patient_id"].astype(str).values,
        "tissue_type": tnk.obs["tissue_type"].astype(str).values,
        "leiden_1.0": tnk.obs["leiden_1.0"].astype(str).values,
        "leiden_1.5": tnk.obs["leiden_1.5"].astype(str).values,
        "tnk_subtype": tnk.obs["tnk_subtype"].astype(str).values,
        "UMAP1": tnk.obsm["X_umap"][:, 0],
        "UMAP2": tnk.obsm["X_umap"][:, 1],
    })
    umap_df.to_csv(FIG / "tnk_umap_metadata.csv.gz", index=False,
                    compression="gzip")

    # Save h5ad (minus raw)
    log(f"saving {OUT_H5}")
    tnk_out = tnk.copy()
    if tnk_out.raw is not None:
        del tnk_out.raw
    tnk_out.write(OUT_H5)

    log("writing summary md")
    with open(RES / "step11_tnk_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 11 — LUAD T_NK subclustering\n\n")
        f.write(f"- Total T_NK cells: {tnk.n_obs}\n")
        f.write(f"- Leiden 1.5 clusters: {tnk.obs['leiden_1.5'].nunique()}\n")
        f.write(f"- Assigned subtypes: {tnk.obs['tnk_subtype'].nunique()}\n")
        f.write(f"- Ambiguous clusters (margin<0.05): {len(ambig)}\n\n")
        f.write("## Subtype cell counts\n\n")
        f.write(stype_counts.to_frame("n_cells").to_markdown() + "\n\n")
        f.write("## Cluster → panel-best vote (Leiden 1.5) with marker hints\n\n")
        f.write(vote.round(3).to_markdown(index=False) + "\n\n")
        f.write("## Top-15 MP ↔ subtype % Spearman (patient-level)\n\n")
        top_any = (mp_corr.sort_values("spearman_rho",
                                        ascending=False, key=abs).head(20))
        f.write(top_any.round(4).to_markdown(index=False) + "\n\n")
        f.write("## Thresholds used in rescue rules\n\n")
        f.write(f"- FOXP3 ≥ {FOXP3_TH:.2f}, IL2RA ≥ {IL2RA_TH:.2f} (Treg)\n")
        f.write(f"- CD8 ≥ {CD8_TH:.2f} (CD8 lineage); CD4 ≥ {CD4_TH:.2f}\n")
        f.write(f"- PDCD1 ≥ {EXH_TH:.2f} (exhaustion)\n")
        f.write(f"- TRGC2 ≥ {TRGC2_TH:.2f} (gdT)\n")
        f.write(f"- TOP2A ≥ {TOP2A_TH:.2f} (proliferation)\n")

    with open(RES / "step11_diagnosis.md", "w", encoding="utf-8") as f:
        f.write("# Step 11 diagnosis — marker availability and ambiguous clusters\n\n")
        f.write("## Marker-panel availability post-HVG\n\n")
        info_df = pd.DataFrame(info).T
        info_df.columns = ["n_var", "n_marker_counts"]
        f.write(info_df.to_markdown() + "\n\n")
        f.write("## Ambiguous clusters (margin<0.05)\n\n")
        if len(ambig):
            f.write(ambig.round(3).to_markdown(index=False) + "\n\n")
        else:
            f.write("None.\n\n")
        f.write("## Treg/Exhausted separation check\n\n")
        ex = tnk.obs[tnk.obs["tnk_subtype"].isin(["Treg", "CD8_Exhausted"])]
        if len(ex):
            f.write(f"- Treg cells: {int((ex['tnk_subtype']=='Treg').sum())}\n")
            f.write(f"- CD8_Exhausted cells: {int((ex['tnk_subtype']=='CD8_Exhausted').sum())}\n")
            # cross with checkpoint co-expression in profile table above
            f.write("- Checkpoint marker expression per subtype: see "
                    "fig_tnk/tnk_exhausted_treg_profile.csv\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
