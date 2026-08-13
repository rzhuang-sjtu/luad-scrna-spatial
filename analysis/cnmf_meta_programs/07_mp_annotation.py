"""
Step 7: MP functional annotation + cell scoring.

Inputs:
    results/step6_gep_mp_assignment.csv
    results/step6_mp_signatures_top100.csv
    data/cnmf_output/gep_pool_zscore.csv
    data/processed/luad_copykat.h5ad

Outputs:
    results/step7_mp_hallmark_gsea.csv     (GSEA prerank, Hallmark 2020)
    results/step7_mp_enrichr_hallmark.csv  (Enrichr overrep on top-100)
    results/step7_mp_enrichr_gobp.csv      (Enrichr overrep GO BP 2023)
    results/step7_mp_cell_scores.csv       (per-cell MP scores + dominant)
    results/step7_mp_umap.pdf              (malignant-only UMAP panels)
    results/step7_mp_dataset_composition.csv
    results/step7_mp_summary.md
    data/processed/luad_malignant_scored.h5ad   (subset + scores + UMAP)
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
os.environ["MPLBACKEND"] = "Agg"

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("${PROJECT_ROOT}")
RESULTS = ROOT / "results"
PROCESSED = ROOT / "data" / "processed"
CNMF_OUTPUT = ROOT / "data" / "cnmf_output"

ASSIGN_CSV = RESULTS / "step6_gep_mp_assignment.csv"
SIGS_CSV = RESULTS / "step6_mp_signatures_top100.csv"
POOL_CSV = CNMF_OUTPUT / "gep_pool_zscore.csv"
H5_IN = PROCESSED / "luad_copykat.h5ad"

OUT_GSEA = RESULTS / "step7_mp_hallmark_gsea.csv"
OUT_ENR_H = RESULTS / "step7_mp_enrichr_hallmark.csv"
OUT_ENR_GO = RESULTS / "step7_mp_enrichr_gobp.csv"
OUT_SCORES = RESULTS / "step7_mp_cell_scores.csv"
OUT_UMAP = RESULTS / "step7_mp_umap.pdf"
OUT_COMP = RESULTS / "step7_mp_dataset_composition.csv"
OUT_SUMMARY = RESULTS / "step7_mp_summary.md"
OUT_H5 = PROCESSED / "luad_malignant_scored.h5ad"

IMMUNE_MARKER_LABELS = {"T_NK", "Myeloid", "Mast", "Plasma"}
DOUBLET_SCORE_THRESHOLD = 0.25

N_TOP_FOR_SCORING = 50
N_TOP_FOR_ENRICHR = 100

HALLMARK_GMT = ROOT / "data" / "gmt" / "MSigDB_Hallmark_2020.gmt"
GOBP_GMT = ROOT / "data" / "gmt" / "GO_Biological_Process_2023.gmt"
HALLMARK_ENRICHR = "MSigDB_Hallmark_2020"
GOBP_ENRICHR = "GO_Biological_Process_2023"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mp_full_ranking(pool: pd.DataFrame, assign: pd.DataFrame) -> dict:
    """Per MP: mean z-score across member GEPs for every gene -> ranked Series."""
    ranks = {}
    for mp, sub in assign.groupby("MP"):
        member = sub["gep_id"].tolist()
        mean_scores = pool[member].mean(axis=1)
        mean_scores = mean_scores.sort_values(ascending=False)
        ranks[mp] = mean_scores
    return ranks


def run_prerank(mp_ranks: dict) -> pd.DataFrame:
    import gseapy as gp
    rows = []
    for mp, series in mp_ranks.items():
        log(f"  GSEA prerank: {mp} ({len(series)} genes)")
        df = pd.DataFrame({"gene": series.index, "rank": series.values})
        try:
            pre = gp.prerank(
                rnk=df,
                gene_sets=str(HALLMARK_GMT),
                outdir=None,
                min_size=5,
                max_size=1000,
                permutation_num=1000,
                seed=42,
                no_plot=True,
                verbose=False,
            )
            r = pre.res2d.copy()
            r["MP"] = mp
            rows.append(r)
        except Exception as e:
            log(f"    WARN {mp}: {e}")
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    keep = ["MP", "Term", "ES", "NES", "NOM p-val", "FDR q-val",
            "Gene %", "Lead_genes"]
    keep = [c for c in keep if c in out.columns]
    return out[keep]


def run_enrichr(sigs: pd.DataFrame, gene_set: str, label: str) -> pd.DataFrame:
    import gseapy as gp
    rows = []
    for mp, sub in sigs.groupby("MP"):
        genes = sub.sort_values("rank").head(N_TOP_FOR_ENRICHR)["gene"].tolist()
        log(f"  Enrichr {label}: {mp} ({len(genes)} genes)")
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_set,
                organism="human",
                outdir=None,
                no_plot=True,
                verbose=False,
            )
            r = enr.results.copy()
            r["MP"] = mp
            rows.append(r)
        except Exception as e:
            log(f"    WARN {mp}: {e}")
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    keep = ["MP", "Term", "Overlap", "P-value", "Adjusted P-value",
            "Odds Ratio", "Combined Score", "Genes"]
    keep = [c for c in keep if c in out.columns]
    return out[keep]


def subset_malignant(ad):
    mal = ad.obs["malignant"].astype(str) == "Malignant"
    immune_hit = ad.obs["celltype_marker"].astype(str).isin(IMMUNE_MARKER_LABELS)
    doublet_hit = ad.obs["doublet_score"].astype(float) > DOUBLET_SCORE_THRESHOLD
    keep = mal & ~(immune_hit | doublet_hit)
    log(f"Malignant subset: {int(keep.sum()):,} / {ad.shape[0]:,} cells")
    return ad[keep].copy()


def score_cells(ad_mal, sigs: pd.DataFrame):
    """sc.tl.score_genes per MP using top-N genes from consensus signatures."""
    for mp, sub in sigs.groupby("MP"):
        genes = sub.sort_values("rank").head(N_TOP_FOR_SCORING)["gene"].tolist()
        genes = [g for g in genes if g in ad_mal.var_names]
        log(f"  score {mp}: {len(genes)} genes")
        sc.tl.score_genes(
            ad_mal, gene_list=genes, score_name=f"{mp}_score",
            random_state=42, n_bins=25, ctrl_size=min(50, len(genes)),
        )
    import re
    mp_cols = sorted(
        [c for c in ad_mal.obs.columns if re.fullmatch(r"MP\d+_score", c)],
        key=lambda c: int(c.replace("MP", "").replace("_score", "")),
    )
    score_mat = ad_mal.obs[mp_cols].values
    dominant_idx = np.argmax(score_mat, axis=1)
    ad_mal.obs["dominant_MP"] = [mp_cols[i].replace("_score", "")
                                  for i in dominant_idx]
    ad_mal.obs["dominant_MP_score"] = score_mat[np.arange(len(score_mat)),
                                                dominant_idx]


def build_malignant_umap(ad_mal):
    """PCA -> Harmony on dataset -> UMAP, in-place on ad_mal."""
    log("Normalising (log1p) from layers['counts'] ...")
    import anndata as ad_mod
    ad2 = ad_mod.AnnData(
        X=ad_mal.layers["counts"].copy(),
        obs=ad_mal.obs.copy(),
        var=ad_mal.var.copy(),
    )
    sc.pp.normalize_total(ad2, target_sum=1e4)
    sc.pp.log1p(ad2)
    log("HVG selection (seurat_v3 on counts) ...")
    sc.pp.highly_variable_genes(ad2, n_top_genes=3000, flavor="seurat_v3",
                                 layer=None, subset=False,
                                 batch_key="dataset")
    ad2_hvg = ad2[:, ad2.var["highly_variable"]].copy()
    log(f"HVG subset: {ad2_hvg.shape}")
    sc.pp.scale(ad2_hvg, max_value=10)
    log("PCA ...")
    sc.tl.pca(ad2_hvg, n_comps=50, random_state=42)
    log("Harmony on dataset (CPU) ...")
    import harmonypy as hm
    ho = hm.run_harmony(
        ad2_hvg.obsm["X_pca"],
        ad2_hvg.obs,
        vars_use=["dataset"],
        max_iter_harmony=20,
        device="cpu",
    )
    # harmonypy 0.2.0 returns Z_corr as (n_cells, n_pcs) already — do NOT transpose
    Zc = np.asarray(ho.Z_corr)
    if Zc.shape[0] != ad2_hvg.shape[0]:
        Zc = Zc.T
    ad2_hvg.obsm["X_pca_harmony"] = Zc.astype(np.float32)
    log("Neighbours + UMAP ...")
    sc.pp.neighbors(ad2_hvg, use_rep="X_pca_harmony", n_neighbors=30,
                    random_state=42)
    sc.tl.umap(ad2_hvg, random_state=42, min_dist=0.3)
    ad_mal.obsm["X_umap_mal"] = ad2_hvg.obsm["X_umap"]
    ad_mal.obsm["X_pca_mal_harmony"] = ad2_hvg.obsm["X_pca_harmony"]


def seurat_v3_hvg_may_fail(ad2):
    try:
        sc.pp.highly_variable_genes(ad2, n_top_genes=3000, flavor="seurat_v3",
                                     batch_key="dataset")
        return True
    except Exception:
        return False


def plot_umap(ad_mal, mp_cols, out_pdf):
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf) as pdf:
        # Page 1: UMAP coloured by dominant MP
        fig, ax = plt.subplots(figsize=(7, 6))
        sc.pl.embedding(ad_mal, basis="X_umap_mal", color="dominant_MP",
                        ax=ax, show=False, frameon=False, legend_loc="right margin")
        ax.set_title("Malignant UMAP — dominant MP")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2..: one per MP score
        ncols = 3
        nrows = int(np.ceil(len(mp_cols) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
        axes = np.array(axes).ravel()
        for i, col in enumerate(mp_cols):
            sc.pl.embedding(ad_mal, basis="X_umap_mal", color=col,
                            ax=axes[i], show=False, frameon=False,
                            cmap="viridis", colorbar_loc=None)
            axes[i].set_title(col)
        for ax in axes[len(mp_cols):]:
            ax.set_visible(False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3: dominant-MP composition by dataset
        comp = (
            ad_mal.obs.groupby(["dataset", "dominant_MP"], observed=True)
            .size()
            .unstack(fill_value=0)
        )
        comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100
        fig, ax = plt.subplots(figsize=(9, 5))
        comp_pct.plot.bar(stacked=True, ax=ax, width=0.85, colormap="tab10")
        ax.set_ylabel("% Malignant cells")
        ax.set_title("Dominant-MP composition per dataset")
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def main():
    log(f"Reading {ASSIGN_CSV}")
    assign = pd.read_csv(ASSIGN_CSV)
    log(f"Reading {SIGS_CSV}")
    sigs = pd.read_csv(SIGS_CSV)
    log(f"Reading {POOL_CSV}")
    pool = pd.read_csv(POOL_CSV, index_col=0)

    mp_list = sorted(assign["MP"].unique(), key=lambda m: int(m.replace("MP", "")))
    log(f"MPs: {mp_list}")

    # 1. GSEA prerank (Hallmark)
    log("=" * 70)
    log("GSEA prerank vs Hallmark 2020")
    ranks = mp_full_ranking(pool, assign)
    gsea = run_prerank(ranks)
    if len(gsea):
        gsea.to_csv(OUT_GSEA, index=False)
        log(f"Wrote {OUT_GSEA} ({len(gsea)} rows)")
    else:
        log("GSEA prerank produced no results")

    # 2. Enrichr overrepresentation
    log("=" * 70)
    if OUT_ENR_H.exists():
        enr_h = pd.read_csv(OUT_ENR_H)
        log(f"Reused {OUT_ENR_H} ({len(enr_h)} rows)")
    else:
        log("Enrichr overrep on top-100 genes (Hallmark 2020)")
        enr_h = run_enrichr(sigs, HALLMARK_ENRICHR, "Hallmark")
        if len(enr_h):
            enr_h.to_csv(OUT_ENR_H, index=False)
            log(f"Wrote {OUT_ENR_H} ({len(enr_h)} rows)")

    if OUT_ENR_GO.exists():
        enr_go = pd.read_csv(OUT_ENR_GO)
        log(f"Reused {OUT_ENR_GO} ({len(enr_go)} rows)")
    else:
        log("Enrichr overrep on top-100 genes (GO BP 2023)")
        enr_go = run_enrichr(sigs, GOBP_ENRICHR, "GO_BP")
        if len(enr_go):
            enr_go.to_csv(OUT_ENR_GO, index=False)
            log(f"Wrote {OUT_ENR_GO} ({len(enr_go)} rows)")

    # 3. Cell scoring on Malignant subset
    log("=" * 70)
    log(f"Reading {H5_IN}")
    ad = sc.read_h5ad(H5_IN)
    ad_mal = subset_malignant(ad)
    del ad
    if "counts" not in ad_mal.layers:
        raise RuntimeError("layers['counts'] missing on malignant subset")

    log("Computing per-cell MP scores ...")
    # Normalise X for scoring
    log("Normalise X (log1p total-count) for score_genes ...")
    ad_mal.layers["lognorm"] = ad_mal.layers["counts"].copy()
    sc.pp.normalize_total(ad_mal, layer="lognorm", target_sum=1e4)
    sc.pp.log1p(ad_mal, layer="lognorm")
    ad_mal.X = ad_mal.layers["lognorm"]
    score_cells(ad_mal, sigs)

    # 4. Malignant-only UMAP
    log("=" * 70)
    log("Building malignant-only UMAP (PCA + harmony + UMAP) ...")
    build_malignant_umap(ad_mal)

    # 5. Export cell-level scores table (compact)
    log("Writing cell-score table ...")
    score_cols = [f"{mp}_score" for mp in mp_list]
    out_tab = ad_mal.obs[
        ["dataset", "patient_id", "sample_id", "tissue_type",
         "celltype_marker", "doublet_score",
         "dominant_MP", "dominant_MP_score"] + score_cols
    ].copy()
    out_tab.index.name = "cell_id"
    out_tab.to_csv(OUT_SCORES)
    log(f"Wrote {OUT_SCORES} ({len(out_tab)} rows)")

    comp = (
        ad_mal.obs.groupby(["dataset", "dominant_MP"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    comp.to_csv(OUT_COMP)
    log(f"Wrote {OUT_COMP}")

    # 6. Plots
    log("Plotting UMAP ...")
    plot_umap(ad_mal, score_cols, OUT_UMAP)
    log(f"Wrote {OUT_UMAP}")

    # 7. Save annotated h5ad (lighter: drop lognorm layer)
    log("Writing scored h5ad ...")
    if "lognorm" in ad_mal.layers:
        del ad_mal.layers["lognorm"]
    ad_mal.write_h5ad(OUT_H5, compression="gzip")
    log(f"Wrote {OUT_H5}")

    # 8. Markdown summary
    log("Writing summary ...")
    with open(OUT_SUMMARY, "w") as f:
        f.write("# Step 7 — MP functional annotation\n\n")
        f.write(f"- Malignant cells: {ad_mal.shape[0]:,}\n")
        f.write(f"- MPs: {', '.join(mp_list)}\n\n")
        f.write("## Dominant-MP distribution\n\n```\n")
        f.write(ad_mal.obs["dominant_MP"].value_counts().to_string())
        f.write("\n```\n\n")
        f.write("## Dataset × dominant-MP counts\n\n```\n")
        f.write(comp.to_string())
        f.write("\n```\n\n")
        if len(gsea):
            f.write("## Top Hallmark GSEA hits per MP (FDR < 0.25)\n\n")
            for mp in mp_list:
                sub = gsea[(gsea["MP"] == mp) &
                           (gsea["FDR q-val"].astype(float) < 0.25)]
                sub = sub.sort_values("NES", ascending=False).head(10)
                f.write(f"### {mp}\n\n```\n")
                if len(sub):
                    f.write(sub[["Term", "NES", "FDR q-val"]].to_string(index=False))
                else:
                    f.write("(no Hallmark term with FDR<0.25)")
                f.write("\n```\n\n")
        if len(enr_h):
            f.write("## Top Enrichr Hallmark per MP (by Combined Score)\n\n")
            for mp in mp_list:
                sub = enr_h[enr_h["MP"] == mp].sort_values(
                    "Combined Score", ascending=False).head(10)
                f.write(f"### {mp}\n\n```\n")
                f.write(sub[["Term", "Adjusted P-value", "Combined Score"]]
                        .to_string(index=False))
                f.write("\n```\n\n")
    log(f"Wrote {OUT_SUMMARY}")

    # 9. Print terse summary
    log("=" * 70)
    log("Dominant-MP distribution:")
    log(ad_mal.obs["dominant_MP"].value_counts().to_string())
    log("")
    if len(gsea):
        log("Top Hallmark per MP (by NES, FDR<0.25):")
        for mp in mp_list:
            sub = gsea[(gsea["MP"] == mp) &
                       (gsea["FDR q-val"].astype(float) < 0.25)]
            sub = sub.sort_values("NES", ascending=False).head(5)
            if len(sub):
                log(f"  {mp}: " + " | ".join(
                    f"{t} (NES={n:.2f}, q={q:.3g})"
                    for t, n, q in zip(sub["Term"], sub["NES"],
                                        sub["FDR q-val"].astype(float))))
            else:
                log(f"  {mp}: (no FDR<0.25 hits)")


if __name__ == "__main__":
    main()
