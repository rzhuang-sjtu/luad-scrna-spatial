"""Step 12: LIANA cell-cell communication analysis of LUAD tumor ecosystem.

Objective: show MP3-high malignant cells recruit/interact with
Neutrophil + Macro_SPP1 (or other TAMs) via SPP1/TGFB1/CSF1/etc. axes,
and separately compare chemokine ligand expression across MP1-4.

Pipeline:
  1. Assemble labeled AnnData: malignant cells by dominant_MP (MP1-4),
     myeloid by myeloid_subtype_refined, T_NK by tnk_subtype, others
     by celltype_coarse.
  2. Downsample each label to max CAP_PER_LABEL cells.
  3. Inject PDCD1 from obsm['marker_counts'] into X to rescue PD-1/PD-L1.
  4. Run LIANA rank_aggregate (5 methods) on consensus LR resource.
  5. Extract MP3→Neutrophil / MP3→Macro_SPP1 LR subsets; compare
     MP1/MP2/MP3/MP4 → Neutrophil strengths.
  6. Separate analysis: mean chemokine ligand expression per MP label
     (CXCL1/2/5/8, CCL2/5, CSF1/2/3, TGFB1, SPP1).
"""

from __future__ import annotations
import os, gc, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

IN_FULL = Path.home() / "luad/data/processed/luad_copykat.h5ad"
IN_MY   = Path.home() / "luad/data/processed/luad_myeloid.h5ad"
IN_TNK  = Path.home() / "luad/data/processed/luad_tnk.h5ad"
MP_SCORES = Path.home() / "luad/results/step7_mp_cell_scores.csv"
RES = Path.home() / "luad/results"
FIG = Path("${WORK_ROOT}/luad_figures/fig6")
FIG.mkdir(parents=True, exist_ok=True)

CAP_PER_LABEL = 1500
CHEMOKINE_LIGANDS = ["CXCL1", "CXCL2", "CXCL8", "CXCL9", "CXCL10", "CXCL11",
                    "CCL2", "CCL3", "CCL4", "CCL5",
                    "CSF1", "TGFB1", "TGFB3", "SPP1", "VEGFA", "MMP9"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    t0 = time.time()

    # --- 1. Build cell_label on each source file; collect barcodes ---
    log("loading MP scores (malignant assignments)")
    mp = pd.read_csv(MP_SCORES, index_col=0)
    mp_labels = mp["dominant_MP"].astype(str).map(lambda x: f"Malignant_{x}")
    mp_labels_dict = mp_labels.to_dict()      # barcode → "Malignant_MP1"/etc
    log(f"  malignant labeled: {len(mp_labels_dict)}")

    log("loading Myeloid labels")
    my = sc.read_h5ad(IN_MY, backed="r")
    # Prefer refined subtype if present
    col = "myeloid_subtype_refined" if "myeloid_subtype_refined" in my.obs.columns \
          else "myeloid_subtype"
    my_labels = my.obs[col].astype(str).to_dict()
    log(f"  myeloid labeled: {len(my_labels)} using {col}")
    # drop unresolved/unknown if any
    keep_my = {k: v for k, v in my_labels.items()
               if v not in ("Unknown", "T_unresolved")}
    log(f"  after filtering unresolved: {len(keep_my)}")

    log("loading T_NK labels")
    tnk = sc.read_h5ad(IN_TNK, backed="r")
    tnk_labels = tnk.obs["tnk_subtype"].astype(str).to_dict()
    keep_tnk = {k: v for k, v in tnk_labels.items()
                if v not in ("T_unresolved", "T_cytotoxic_ambiguous", "CD8_other")}
    log(f"  T_NK labeled: {len(tnk_labels)}, kept: {len(keep_tnk)}")

    # --- 2. Load main h5ad in memory and build labels ---
    log("loading full h5ad into memory")
    a = sc.read_h5ad(IN_FULL)
    log(f"  shape={a.shape}")

    labels = pd.Series("_drop_", index=a.obs.index, dtype=object, name="cell_label")
    # Non-malignant, non-myeloid, non-T_NK → keep by coarse celltype
    coarse = a.obs["celltype_coarse"].astype(str)
    keep_coarse = {"Endothelial", "Fibroblast", "B", "Plasma", "Mast"}
    coarse_mask = coarse.isin(keep_coarse)
    labels.loc[coarse_mask] = coarse[coarse_mask].values

    # Overwrite with malignant labels
    bc_mp = [bc for bc in mp_labels_dict if bc in labels.index]
    labels.loc[bc_mp] = pd.Series(mp_labels_dict)[bc_mp].values
    # Overwrite with myeloid
    bc_my = [bc for bc in keep_my if bc in labels.index]
    labels.loc[bc_my] = pd.Series(keep_my)[bc_my].values
    # Overwrite with T_NK
    bc_tnk = [bc for bc in keep_tnk if bc in labels.index]
    labels.loc[bc_tnk] = pd.Series(keep_tnk)[bc_tnk].values

    a.obs["cell_label"] = labels.values
    log("label counts (pre-downsample):")
    vc = a.obs["cell_label"].value_counts()
    log(vc.to_string())

    # --- 3. Downsample per label ---
    log(f"downsampling to max {CAP_PER_LABEL} per label")
    rng = np.random.default_rng(0)
    keep_idx = []
    for lbl, sub in a.obs.groupby("cell_label", observed=True):
        if lbl == "_drop_":
            continue
        if len(sub) <= CAP_PER_LABEL:
            keep_idx.extend(sub.index.tolist())
        else:
            keep_idx.extend(rng.choice(sub.index.values, size=CAP_PER_LABEL,
                                         replace=False).tolist())
    log(f"  retained cells: {len(keep_idx)}")

    ad_sub = a[keep_idx].copy()
    del a; gc.collect()
    log(f"  subset shape: {ad_sub.shape}")
    log("post-downsample label counts:")
    log(ad_sub.obs["cell_label"].value_counts().to_string())

    # --- 4. Use lognorm layer as X (LIANA wants log1p-normalized) ---
    log("set X = layers['lognorm']")
    ad_sub.X = ad_sub.layers["lognorm"].astype(np.float32).copy()

    # --- 5. Inject PDCD1 from marker_counts into var ---
    log("injecting PDCD1 (and CD8A/CD8B/FOXP3) from marker_counts into X")
    mc = ad_sub.obsm["marker_counts"]
    rescue_genes = [g for g in ["PDCD1", "CD8A", "CD8B", "FOXP3"]
                    if g in mc.columns and g not in ad_sub.var_names]
    log(f"  rescuing: {rescue_genes}")
    if rescue_genes:
        # log1p the raw marker counts for compatibility with lognorm X
        extra = np.log1p(mc[rescue_genes].values.astype(np.float32))
        from scipy.sparse import hstack, csr_matrix, issparse
        X_old = ad_sub.X
        X_extra = csr_matrix(extra)
        if issparse(X_old):
            X_new = hstack([X_old, X_extra]).tocsr()
        else:
            X_new = np.hstack([X_old, extra]).astype(np.float32)
        new_var = pd.concat([
            ad_sub.var,
            pd.DataFrame(index=rescue_genes,
                          data={c: pd.NA for c in ad_sub.var.columns}),
        ])
        # Need a fresh AnnData to rebuild (can't reassign X with different n_vars)
        ad_sub = ad.AnnData(
            X=X_new,
            obs=ad_sub.obs,
            var=new_var,
            obsm={k: v for k, v in ad_sub.obsm.items()
                  if k in ("X_umap", "X_pca")},  # keep minimal
            uns={},
        )
        log(f"  new shape: {ad_sub.shape}")

    # --- 5b. Sanitize X: replace NaN/Inf with 0 ---
    log("sanitizing X for non-finite values")
    from scipy.sparse import issparse
    if issparse(ad_sub.X):
        data = ad_sub.X.data
        bad = ~np.isfinite(data)
        if bad.any():
            log(f"  NaN/Inf entries: {int(bad.sum())} / {data.size}")
            ad_sub.X.data = np.where(bad, 0.0, data).astype(np.float32)
            ad_sub.X.eliminate_zeros()
    else:
        bad = ~np.isfinite(ad_sub.X)
        if bad.any():
            log(f"  NaN/Inf entries: {int(bad.sum())} / {ad_sub.X.size}")
            ad_sub.X = np.where(bad, 0.0, ad_sub.X).astype(np.float32)
    log("  sanitize done")

    # --- 6. Run LIANA ---
    log("running LIANA rank_aggregate (5 methods)")
    import liana
    import liana.mt as lmt
    log(f"  liana version: {liana.__version__}")
    lmt.rank_aggregate(ad_sub,
                       groupby="cell_label",
                       resource_name="consensus",
                       use_raw=False,
                       expr_prop=0.05,       # require ≥5% cells expressing each
                       verbose=True,
                       seed=0)

    res = ad_sub.uns["liana_res"]
    log(f"  LIANA result: {res.shape}")
    res.to_csv(RES / "step12_liana_all.csv", index=False)

    # --- 7. Extract MP3 focused subsets ---
    log("extracting MP3-centric subsets")
    src_col = "source"
    tgt_col = "target"
    # LIANA column naming may vary; check
    log(f"  columns: {list(res.columns)[:15]}")

    # MP3 → Neutrophil
    mp3_neut = res[(res[src_col] == "Malignant_MP3") &
                    (res[tgt_col] == "Neutrophil")].copy()
    mp3_neut.sort_values("magnitude_rank", inplace=True) \
        if "magnitude_rank" in mp3_neut.columns else None
    mp3_neut.to_csv(RES / "step12_liana_mp3_neutrophil.csv", index=False)
    log(f"  MP3 → Neutrophil pairs: {len(mp3_neut)}")

    # 4 MPs → Neutrophil comparison
    mp_neut = res[(res[src_col].str.startswith("Malignant_MP")) &
                   (res[tgt_col] == "Neutrophil")].copy()
    mp_neut.to_csv(RES / "step12_liana_mp_comparison.csv", index=False)

    # MP3 → Macro_SPP1
    mp3_spp1 = res[(res[src_col] == "Malignant_MP3") &
                    (res[tgt_col] == "Macro_SPP1")].copy() \
               if "Macro_SPP1" in res[tgt_col].unique() else pd.DataFrame()
    if len(mp3_spp1):
        mp3_spp1.to_csv(RES / "step12_liana_mp3_macroSPP1.csv", index=False)

    # Reverse: Macro_SPP1 → Malignant_MP3
    spp1_mp3 = res[(res[src_col] == "Macro_SPP1") &
                    (res[tgt_col] == "Malignant_MP3")].copy() \
               if "Macro_SPP1" in res[src_col].unique() else pd.DataFrame()
    if len(spp1_mp3):
        spp1_mp3.to_csv(RES / "step12_liana_macroSPP1_mp3.csv", index=False)

    # --- 8. Source-side chemokine ligand comparison (independent of LR pairing) ---
    log("source-side chemokine ligand expression per MP label")
    # LUAD-like malignant labels only
    malignant_labels = [l for l in ad_sub.obs["cell_label"].unique()
                        if str(l).startswith("Malignant_")]
    lig_present = [g for g in CHEMOKINE_LIGANDS if g in ad_sub.var_names]
    log(f"  ligands present: {len(lig_present)}/{len(CHEMOKINE_LIGANDS)} -> {lig_present}")
    rows = []
    for lbl in malignant_labels:
        mask = (ad_sub.obs["cell_label"] == lbl).values
        if not mask.any(): continue
        X = ad_sub[mask, lig_present].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        for i, g in enumerate(lig_present):
            rows.append({"label": lbl, "ligand": g,
                         "mean_log1p": float(X[:, i].mean()),
                         "pct_expressing": float((X[:, i] > 0).mean()),
                         "n_cells": int(mask.sum())})
    lig_df = pd.DataFrame(rows)
    lig_pivot = lig_df.pivot(index="ligand", columns="label", values="mean_log1p")
    lig_df.to_csv(FIG / "malignant_chemokine_ligand_expression.csv", index=False)
    lig_pivot.to_csv(RES / "step12_mp_chemokine_ligand_mean.csv")
    log("  mean log1p ligand expression per MP:")
    log(lig_pivot.round(2).to_string())

    # --- 9. LIANA figure-ready tables ---
    log("figure-ready exports")
    # dotplot: top N LR pairs per source (restricted to malignant→immune)
    dot = res.copy()
    if "magnitude_rank" in dot.columns:
        dot["rank_key"] = dot["magnitude_rank"]
    elif "specificity_rank" in dot.columns:
        dot["rank_key"] = dot["specificity_rank"]
    else:
        dot["rank_key"] = dot.index

    immune_targets = {"Neutrophil", "Macro_SPP1", "Macro_FOLR2", "Macro_C1QC",
                       "Macro_FCN1", "CD8_Effector", "CD8_Exhausted", "Treg",
                       "NK", "cDC_LAMP3", "cDC2", "cDC1"}
    malignant_sources = {f"Malignant_MP{i}" for i in range(1, 5)}
    sel = dot[dot[src_col].isin(malignant_sources) & dot[tgt_col].isin(immune_targets)]
    sel.sort_values("rank_key").to_csv(FIG / "liana_dotplot_data.csv", index=False)

    # MP3-specific (rank lower in MP3 than in other MPs for same ligand-target)
    # Wide format for comparison
    pivot_key = ["ligand_complex", "receptor_complex", "target"] \
                if "ligand_complex" in dot.columns else ["ligand", "receptor", "target"]
    pivot_rank = dot[dot[src_col].isin(malignant_sources)].pivot_table(
        index=pivot_key, columns=src_col, values="rank_key", aggfunc="min")
    pivot_rank = pivot_rank.rename_axis(columns=None).reset_index()
    pivot_rank.to_csv(FIG / "liana_mp3_specific.csv", index=False)

    # Circle: aggregate strength per source→target = count of LR pairs with
    # rank_key <= top 10% as "strong"
    thresh = dot["rank_key"].quantile(0.10)
    strong = dot[dot["rank_key"] <= thresh]
    circle = strong.groupby([src_col, tgt_col]).size().reset_index(name="n_strong_pairs")
    circle.to_csv(FIG / "liana_circle_data.csv", index=False)

    # --- 10. Summary md ---
    log("writing summary")
    with open(RES / "step12_summary.md", "w", encoding="utf-8") as f:
        f.write("# Step 12 — LIANA cell-cell communication\n\n")
        f.write(f"- Retained cells: {ad_sub.n_obs}\n")
        f.write(f"- Labels: {ad_sub.obs['cell_label'].nunique()}\n")
        f.write(f"- LIANA consensus LR result: {len(res)} rows\n")
        f.write(f"- MP3 → Neutrophil LR pairs passing filters: {len(mp3_neut)}\n\n")

        f.write("## Key gene-coverage caveats\n\n")
        f.write("Missing receptors (prevent LR inference): CXCR1, CXCR2, CCR2, "
                "CCR5, CSF2RA, CSF3R, PDCD1LG2. Rescued from marker_counts: "
                f"{rescue_genes}. Working axes: SPP1→CD44/ITGAV/ITGB1, "
                "TGFB1→TGFBR1/TGFBR2, CSF1→CSF1R, CD274(PD-L1)→PDCD1.\n\n")

        f.write("## Malignant-by-MP chemokine ligand source expression "
                "(independent of receptor matching)\n\n")
        f.write(lig_pivot.round(3).to_markdown() + "\n\n")

        if len(mp3_neut):
            topN = (mp3_neut.sort_values("rank_key").head(15)
                    if "rank_key" in mp3_neut.columns else mp3_neut.head(15))
            f.write("## Top LIANA pairs Malignant_MP3 → Neutrophil\n\n")
            cols_show = [c for c in ["ligand", "receptor", "ligand_complex",
                                      "receptor_complex", "magnitude_rank",
                                      "specificity_rank", "lrscore", "expr_prod"]
                         if c in topN.columns]
            f.write(topN[cols_show].round(4).to_markdown(index=False) + "\n\n")

        f.write("## Notes\n")
        f.write("- CXCL1/2/8 ligand expression from malignant cells compared "
                "across MP1-4 gives the cleanest bulk-consistent signal (see "
                "step12_mp_chemokine_ligand_mean.csv); receptor side on "
                "neutrophils cannot be cross-checked because CXCR1/CXCR2 "
                "were dropped in the 9881-HVG filter.\n")

    log(f"DONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
