"""Step 10f: Compute Fig 4 14-panel auxiliary data.

Outputs (fig_treatment-style → fig4/):
  panel_F_antitumor_genes.csv         13-subtype × anti-tumor gene mean+pct
  panel_GM_subset_markers.csv         per-Macro-subset functional marker (long form)
  panel_major_type_metadata.csv.gz    metadata with myeloid_major_type column added
  panel_N_spp1_vs_c1qc_deg.csv        DEG ranked list
  panel_N_spp1_vs_c1qc_gsea.csv       Hallmark GSEA result

Usage:
  python analysis/myeloid_neutrophil/10f_panel_data.py
"""
from __future__ import annotations
import os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

H5 = Path.home() / "luad/data/processed/luad_myeloid.h5ad"
LBL = Path.home() / "luad/results/step10c_obs_labels.csv.gz"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig4")

# ---- Panel F: anti-tumor / pro-immune gene panel ----
# Gene list taken from the HCC study cited in the manuscript.
ANTI_TUMOR = [
    # MHC-II antigen presentation
    "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", "HLA-DPB1", "HLA-DRA", "HLA-DMA",
    "HLA-DMB", "HLA-DPA1", "HLA-DRB1", "CD74",
    # Pro-inflammatory cytokines
    "IL1A", "IL1B", "IL6", "TNF",
    # Th1 polarization
    "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
    # IRF / TF
    "IRF1", "IRF8", "STAT1",
    # M1 chemokines
    "CSF2", "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4",
]

# ---- Panel G-M: per macro subset functional markers ----
SUBSET_MARKERS = {
    "Macro_MARCO":  ["MARCO", "FABP4", "PPARG", "MCEMP1", "CD163", "VSIG4",
                      "MSR1", "C1QA", "APOE"],                    # alveolar mac
    "Macro_C1QC":   ["C1QA", "C1QB", "C1QC", "APOE", "HLA-DRA",
                      "HLA-DPB1", "AXL", "GPR34"],                  # antigen-presenting
    "Macro_FCN1":   ["FCN1", "VCAN", "S100A8", "S100A9", "IL1B",
                      "EREG", "FCGR3A", "CD14", "IRF1", "IRF8"],    # mono-derived inflammatory
    "Macro_FOLR2":  ["FOLR2", "MRC1", "CD163", "SELENOP", "LYVE1",
                      "MAF", "STAB1"],                              # tissue-resident
    "Macro_SPP1":   ["SPP1", "TREM2", "VEGFA", "MMP9", "ADAM8",
                      "FN1", "P4HA1", "LDHA", "ENO1"],              # hypoxic / pro-tumor
    "Macro_general": ["APOE", "C1QA", "CD68", "CSF1R", "CTSB",
                      "CTSD", "CTSL"],                              # baseline
    "Macro_prolif": ["MKI67", "TOP2A", "STMN1", "CDK1", "CCNB1"],  # cell cycle
}

# Major type collapse rule
MAJOR_TYPE_MAP = {
    "Macro_FCN1": "Macrophage", "Macro_C1QC": "Macrophage",
    "Macro_FOLR2": "Macrophage", "Macro_MARCO": "Macrophage",
    "Macro_SPP1": "Macrophage", "Macro_general": "Macrophage",
    "Macro_prolif": "Macrophage",
    "Mono_nonclassical": "Mono_NC",
    "Neutrophil": "Neutrophil",
    "cDC1": "cDC1", "cDC2": "cDC2", "cDC_LAMP3": "cDC_LAMP3",
    "pDC": "pDC",
}

GENE_SETS_GSEA = ["MSigDB_Hallmark_2020"]


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main():
    t0 = time.time()
    FIG.mkdir(parents=True, exist_ok=True)

    log(f"load {H5}")
    a = sc.read_h5ad(H5)
    log(f"  shape={a.shape}")

    # Apply refined labels
    log(f"apply refined labels {LBL}")
    lbl = pd.read_csv(LBL, index_col=0)
    a.obs["myeloid_subtype_refined"] = (
        lbl["myeloid_subtype_refined"].reindex(a.obs.index).astype("category")
    )
    a.obs["myeloid_major_type"] = (
        a.obs["myeloid_subtype_refined"].astype(str).map(MAJOR_TYPE_MAP)
        .fillna(a.obs["myeloid_subtype_refined"].astype(str)).astype("category")
    )
    log("  major-type counts:")
    log(a.obs["myeloid_major_type"].value_counts().to_string())

    # ---- Save updated metadata (Panel A/B + downstream) ----
    umap = a.obsm["X_umap"]
    md = pd.DataFrame({
        "barcode": a.obs.index,
        "dataset": a.obs["dataset"].astype(str).values,
        "patient_id": a.obs["patient_id"].astype(str).values,
        "tissue_type": a.obs["tissue_type"].astype(str).values,
        "myeloid_subtype_refined": a.obs["myeloid_subtype_refined"].astype(str).values,
        "myeloid_major_type": a.obs["myeloid_major_type"].astype(str).values,
        "UMAP1": umap[:, 0], "UMAP2": umap[:, 1],
        "M1_score": a.obs["M1_score"].values,
        "M2_score": a.obs["M2_score"].values,
    })
    md.to_csv(FIG / "panel_major_type_metadata.csv.gz",
              index=False, compression="gzip")
    log(f"  panel_major_type_metadata.csv.gz: {len(md)} rows")

    # ---- Panel F: anti-tumor gene matrix ----
    log("Panel F: anti-tumor gene mean+pct per subtype")
    present = [g for g in ANTI_TUMOR if g in a.var_names]
    log(f"  {len(present)}/{len(ANTI_TUMOR)} genes present")
    X = a[:, present].X
    Xa = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    df_expr = pd.DataFrame(Xa, index=a.obs.index, columns=present)
    df_expr["sub"] = a.obs["myeloid_subtype_refined"].astype(str).values
    rows = []
    for st, sub in df_expr.groupby("sub", sort=False):
        if st in ("nan", ""):
            continue
        mat = sub[present].values
        for i, g in enumerate(present):
            rows.append({
                "subtype": st, "gene": g,
                "mean_log1p": float(mat[:, i].mean()),
                "pct_expressing": float((mat[:, i] > 0).mean()),
            })
    pd.DataFrame(rows).to_csv(FIG / "panel_F_antitumor_genes.csv", index=False)
    log(f"  panel_F: {len(rows)} rows, {len(present)} genes × "
        f"{df_expr['sub'].nunique()} subtypes")

    # ---- Panel G-M: per-subset functional markers ----
    log("Panel G-M: functional markers per Macro subset")
    all_subset_genes = sorted(set(g for v in SUBSET_MARKERS.values() for g in v))
    present_sub = [g for g in all_subset_genes if g in a.var_names]
    log(f"  pooled marker genes: {len(present_sub)}/{len(all_subset_genes)} present")
    Xs = a[:, present_sub].X
    Xsa = Xs.toarray() if hasattr(Xs, "toarray") else np.asarray(Xs)
    df_e = pd.DataFrame(Xsa, index=a.obs.index, columns=present_sub)
    df_e["sub"] = a.obs["myeloid_subtype_refined"].astype(str).values
    rows = []
    for st, sub in df_e.groupby("sub", sort=False):
        if st in ("nan", ""):
            continue
        mat = sub[present_sub].values
        for i, g in enumerate(present_sub):
            rows.append({
                "subtype": st, "gene": g,
                "mean_log1p": float(mat[:, i].mean()),
                "pct_expressing": float((mat[:, i] > 0).mean()),
            })
    # also tag which panel each gene "belongs" to (its origin macro subset)
    panel_assign = {}
    for ms, glist in SUBSET_MARKERS.items():
        for g in glist:
            panel_assign.setdefault(g, ms)  # first-claim wins
    df_gm = pd.DataFrame(rows)
    df_gm["panel_origin"] = df_gm["gene"].map(panel_assign)
    df_gm.to_csv(FIG / "panel_GM_subset_markers.csv", index=False)
    log(f"  panel_GM: {len(df_gm)} rows")

    # ---- Panel N: SPP1 vs C1QC DEG + GSEA Hallmark ----
    log("Panel N: SPP1 vs C1QC rank_genes_groups + Hallmark GSEA")
    if {"Macro_SPP1", "Macro_C1QC"}.issubset(set(a.obs["myeloid_subtype_refined"].cat.categories)):
        sc.tl.rank_genes_groups(
            a, groupby="myeloid_subtype_refined",
            groups=["Macro_SPP1"], reference="Macro_C1QC",
            method="wilcoxon", n_genes=5000, pts=True,
            key_added="SPP1_vs_C1QC",
        )
        deg = sc.get.rank_genes_groups_df(a, group="Macro_SPP1", key="SPP1_vs_C1QC")
        deg = deg.rename(columns={
            "names": "gene", "logfoldchanges": "logFC",
            "pvals": "pval", "pvals_adj": "pval_adj", "scores": "score",
        })
        deg = deg[["gene", "score", "logFC", "pval", "pval_adj"]]
        deg.to_csv(FIG / "panel_N_spp1_vs_c1qc_deg.csv", index=False)
        log(f"  DEG rows: {len(deg)}")

        # Pre-ranked GSEA
        try:
            import gseapy as gp
            ranked = deg.dropna(subset=["score"]).sort_values("score", ascending=False)
            ranked = ranked[["gene", "score"]].copy()
            ranked.columns = [0, 1]
            log("  running prerank GSEA Hallmark")
            res = gp.prerank(rnk=ranked, gene_sets=GENE_SETS_GSEA[0],
                             outdir=None, threads=4, min_size=10, max_size=500,
                             permutation_num=1000, seed=42, verbose=False)
            r = res.res2d
            r.to_csv(FIG / "panel_N_spp1_vs_c1qc_gsea.csv", index=False)
            log(f"  GSEA terms: {len(r)}; top 5 by NES:")
            top5 = r.sort_values("NES", ascending=False).head(5)[["Term", "NES", "FDR q-val"]]
            log(top5.round(3).to_string(index=False))
        except Exception as e:
            log(f"  GSEA failed: {e}")
    else:
        log("  WARN: missing Macro_SPP1 or Macro_C1QC; skip Panel N")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
