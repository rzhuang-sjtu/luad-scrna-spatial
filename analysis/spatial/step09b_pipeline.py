"""
Step 9b: full Fig 7 pipeline on Okamura LUAD No.1-5 cohort.

Reuses step02 inf_aver.csv (LUAD-derived 28 cell-type signatures).
Stages in one process:
  1. cell2location joint deconvolution (full-batch GPU, 30k epochs)
  2. MP1-5 spot scoring (signature score)
  3. PROGENy 14 pathways
  4. COMMOT OSM + IL1
  5. ROI definition (NFkB-high & Neutrophil-high) + ROI vs non-ROI stats
  6. Cell-type and pathway spatial plots per section

MISTy (R) runs in a separate script (step09c_misty.R).

Outputs:
  ${DATA_ROOT}/ST/results/step09_okamura_validation/
    cohort.h5ad                       (already saved by step09a)
    all_sections_c2l.h5ad
    cohort_with_mp.h5ad
    cohort_with_progeny.h5ad
    cohort_with_roi.h5ad
    section_h5ad/<sample>.h5ad        (already saved by step09a)
    commot/<sample>.h5ad
    spatial_plots/<sample>_*.png
    *.csv tables
"""
from __future__ import annotations
import os, time, gc, traceback, json
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("${DATA_ROOT}/ST/results/step09_okamura_validation")
SEC_DIR = ROOT / "section_h5ad"
PLOTS = ROOT / "spatial_plots"
COMMOT_DIR = ROOT / "commot"
for d in (PLOTS, COMMOT_DIR):
    d.mkdir(parents=True, exist_ok=True)
LOG = ROOT / "pipeline.log"
def log(m):
    s=f"[{time.strftime('%H:%M:%S')}] {m}"; print(s,flush=True)
    open(LOG,"a").write(s+"\n")

REF_CSV = Path("${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv")
MP_SIGS = Path(os.path.expanduser("~/luad/results/step6_mp_signatures_top100.csv"))

PATHWAYS_COMMOT = ["OSM", "IL1"]
PATHWAYS_PROGENY_HIGHLIGHT = ["NFkB","JAK-STAT","TNFa","TGFb","Hypoxia","MAPK"]
NEU_COLS = ["Neu_Inflammatory","Neu_OSM_priming","Neu_OSM_low","Neu_IFN_response",
            "Neu_Angiogenic","Neu_Metastatic","Neu_ECM_remodeling"]
TF_GENES = ["ATF3","FOSB","JUN","JUNB","NFKBIA","FOS"]
LIGAND_GENES = ["IL1B","OSM","IL1A","TNF","TGFB1"]


def stage_c2l(cohort_in: ad.AnnData, out_h5: Path) -> ad.AnnData:
    if out_h5.exists():
        log(f"[1/5] c2l output exists: {out_h5}, loading")
        return sc.read_h5ad(str(out_h5))
    log(f"[1/5] cell2location joint deconvolution on {cohort_in.shape}")
    inf_aver = pd.read_csv(REF_CSV, index_col=0)
    common = sorted(set(cohort_in.var_names) & set(inf_aver.index))
    cohort = cohort_in[:, common].copy()
    sub_inf = inf_aver.loc[common].copy()
    log(f"   intersect: {len(common)} genes")

    if not sp.issparse(cohort.X):
        cohort.X = sp.csr_matrix(cohort.X)
    cohort.X = cohort.X.astype("float32")

    import torch
    from cell2location.models import Cell2location
    use_gpu = torch.cuda.is_available()
    log(f"   GPU: {use_gpu}")

    Cell2location.setup_anndata(adata=cohort, batch_key="sample")
    mod = Cell2location(cohort, cell_state_df=sub_inf,
                        N_cells_per_location=30, detection_alpha=20)

    accel = "gpu" if use_gpu else "cpu"
    # 16-sample cohort (53k spots) — full 30k epochs ETA ~48h on RTX 3080.
    # 10k epochs is plenty for ELBO convergence (c2l docs: plateau by 10–15k).
    # tf32 matmul precision adds ~30% throughput.
    if use_gpu:
        torch.set_float32_matmul_precision("high")
    train_kwargs = dict(max_epochs=10000, batch_size=None, train_size=1.0, lr=0.002)
    t0 = time.time()
    try:
        mod.train(accelerator=accel, **train_kwargs)
    except torch.cuda.OutOfMemoryError as e:
        log(f"   [oom] full-batch failed; retrying batch_size=2500")
        torch.cuda.empty_cache(); gc.collect()
        Cell2location.setup_anndata(adata=cohort, batch_key="sample")
        mod = Cell2location(cohort, cell_state_df=sub_inf, N_cells_per_location=30, detection_alpha=20)
        train_kwargs["batch_size"] = 2500
        mod.train(accelerator=accel, **train_kwargs)
    log(f"   training: {(time.time()-t0)/60:.1f} min")

    try:
        fig, ax = plt.subplots(figsize=(5, 3))
        mod.plot_history(1000); fig.tight_layout()
        fig.savefig(ROOT / "elbo.png", dpi=150); plt.close(fig)
    except Exception:
        pass

    cohort = mod.export_posterior(cohort, sample_kwargs={"num_samples": 1000, "batch_size": cohort.n_obs})
    mod.save(str(ROOT / "c2l_model"), overwrite=True)

    # Strip prefix from obsm DataFrames
    import re
    if "q05_cell_abundance_w_sf" in cohort.obsm:
        df = cohort.obsm["q05_cell_abundance_w_sf"].copy()
        df.columns = [re.sub(r"^q05cell_abundance_w_sf_", "", c) for c in df.columns]
        cohort.obsm["q05_cell_abundance"] = df
    if "means_cell_abundance_w_sf" in cohort.obsm:
        df = cohort.obsm["means_cell_abundance_w_sf"].copy()
        df.columns = [re.sub(r"^meanscell_abundance_w_sf_", "", c) for c in df.columns]
        cohort.obsm["mean_cell_abundance"] = df

    cohort.write_h5ad(str(out_h5), compression="gzip")
    log(f"   saved {out_h5}")
    return cohort


def stage_mp(cohort: ad.AnnData, out_h5: Path) -> ad.AnnData:
    if out_h5.exists():
        log(f"[2/5] MP output exists: {out_h5}, loading")
        return sc.read_h5ad(str(out_h5))
    log(f"[2/5] MP1-5 spot scoring")
    cohort.layers["counts"] = cohort.X.copy()
    sc.pp.normalize_total(cohort, target_sum=1e4)
    sc.pp.log1p(cohort)
    sig_df = pd.read_csv(MP_SIGS)
    score_cols = []
    for mp in sorted(sig_df["MP"].unique()):
        genes = sig_df.loc[sig_df["MP"] == mp, "gene"].tolist()
        present = [g for g in genes if g in cohort.var_names]
        col = f"{mp}_score"
        sc.tl.score_genes(cohort, gene_list=present, score_name=col, ctrl_size=200, random_state=0)
        score_cols.append(col)
    main_cols = [f"MP{i}_score" for i in (1,2,3,4)]
    cohort.obs["dominant_MP_4"] = cohort.obs[main_cols].idxmax(axis=1).str.replace("_score","")
    cohort.obs["dominant_MP_4_score"] = cohort.obs[main_cols].max(axis=1)
    cohort.write_h5ad(str(out_h5), compression="gzip")
    means = cohort.obs.groupby("sample", observed=False)[score_cols].mean().round(3)
    means.to_csv(ROOT / "per_sample_mean_mp.csv")
    log(f"   per-sample mean MP:\n{means.to_string()}")
    return cohort


def stage_progeny(cohort: ad.AnnData, out_h5: Path) -> ad.AnnData:
    if out_h5.exists():
        log(f"[3/5] PROGENy output exists: {out_h5}, loading")
        return sc.read_h5ad(str(out_h5))
    log(f"[3/5] PROGENy 14 pathways")
    import decoupler as dc
    progeny = dc.op.progeny(organism="human", top=500)
    res = dc.mt.mlm(data=cohort, net=progeny, verbose=False)
    if isinstance(res, tuple) and len(res) == 2:
        estimate, pvals = res
    else:
        estimate = cohort.obsm.get("score_mlm")
        pvals    = cohort.obsm.get("padj_mlm")
        if estimate is None:
            raise RuntimeError(f"unexpected dc.mt.mlm return: {type(res)} and no obsm['score_mlm']")
    if not isinstance(estimate, pd.DataFrame):
        estimate = pd.DataFrame(estimate, index=cohort.obs_names)
    cohort.obsm["progeny_mlm"] = estimate
    if pvals is not None:
        cohort.obsm["progeny_mlm_pvals"] = pvals
    for col in estimate.columns:
        cohort.obs[f"progeny_{col}"] = estimate[col].values
    estimate.to_csv(ROOT / "spot_progeny_scores.csv")
    means = estimate.copy(); means["sample"] = cohort.obs["sample"].values
    means.groupby("sample").mean().round(3).to_csv(ROOT / "per_sample_mean_progeny.csv")
    cohort.write_h5ad(str(out_h5), compression="gzip")
    log(f"   {estimate.shape[1]} pathways scored")
    return cohort


def stage_commot(samples_to_run):
    """Run COMMOT per section since it's spatial-graph-local."""
    log(f"[4/5] COMMOT OSM + IL1 per section")
    import commot as ct
    df_lr_full = ct.pp.ligand_receptor_database(species="human", signaling_type=None, database="CellChat")
    df_lr = df_lr_full[df_lr_full.iloc[:, 2].isin(PATHWAYS_COMMOT)].copy()
    log(f"   filtered LR pairs: {df_lr.shape}")
    summary = []
    for s in samples_to_run:
        out = COMMOT_DIR / f"{s}.h5ad"
        if out.exists():
            log(f"   {s}: cached")
        else:
            t0 = time.time()
            sec = sc.read_h5ad(str(SEC_DIR / f"{s}.h5ad"))
            sec.var_names_make_unique()
            sec.layers["counts"] = sec.X.copy()
            sc.pp.normalize_total(sec, target_sum=1e4)
            sc.pp.log1p(sec)
            ct.tl.spatial_communication(sec, database_name="cellchat", df_ligrec=df_lr,
                                         dis_thr=500.0, heteromeric=True, pathway_sum=True)
            for p in PATHWAYS_COMMOT:
                try: ct.tl.communication_direction(sec, database_name="cellchat", pathway_name=p, k=5)
                except Exception: pass
            # extract pathway-level sender/receiver
            sender_df = sec.obsm["commot-cellchat-sum-sender"]
            recv_df   = sec.obsm["commot-cellchat-sum-receiver"]
            for p in PATHWAYS_COMMOT:
                s_cols = [c for c in sender_df.columns if c == f"s-{p}" or c.startswith(f"s-{p}-") or
                          (c.startswith("s-") and c.split("-",2)[1] == p)]
                r_cols = [c for c in recv_df.columns   if c == f"r-{p}" or c.startswith(f"r-{p}-") or
                          (c.startswith("r-") and c.split("-",2)[1] == p)]
                sec.obs[f"s_{p}"] = sender_df[s_cols].sum(axis=1).values if s_cols else 0.0
                sec.obs[f"r_{p}"] = recv_df[r_cols].sum(axis=1).values    if r_cols else 0.0
                sec.obs[f"total_{p}"] = sec.obs[f"s_{p}"].fillna(0) + sec.obs[f"r_{p}"].fillna(0)
            sec.write_h5ad(str(out), compression="gzip")
            log(f"   {s}: COMMOT done in {(time.time()-t0)/60:.1f} min")
            del sec; gc.collect()
        # summary
        a = sc.read_h5ad(str(out), backed="r")
        row = {"sample": s, "n_spots": a.n_obs}
        for p in PATHWAYS_COMMOT:
            for k in (f"s_{p}", f"r_{p}"):
                if k in a.obs.columns:
                    row[f"{k}_mean"] = float(a.obs[k].mean())
        summary.append(row)
        a.file.close()
    pd.DataFrame(summary).to_csv(ROOT / "per_sample_pathway_summary.csv", index=False)
    log(f"   summary written")


def stage_roi(cohort: ad.AnnData, out_h5: Path):
    log(f"[5/5] ROI analysis (NFkB-high & Neutrophil-high)")
    abund = cohort.obsm["q05_cell_abundance"].copy()
    neu = [c for c in NEU_COLS if c in abund.columns]
    cohort.obs["neu_total"] = abund[neu].sum(axis=1).values
    cohort.obs["roi"] = False
    cohort.obs["z_nfkb"] = 0.0
    cohort.obs["z_neu"] = 0.0
    samples = sorted(cohort.obs["sample"].unique().tolist())
    for s in samples:
        m = cohort.obs["sample"].values == s
        nf = cohort.obs.loc[m, "progeny_NFkB"].values
        nu = cohort.obs.loc[m, "neu_total"].values
        zn = (nf - nf.mean()) / (nf.std() + 1e-9)
        zu = (nu - nu.mean()) / (nu.std() + 1e-9)
        cohort.obs.loc[m, "z_nfkb"] = zn
        cohort.obs.loc[m, "z_neu"]  = zu
        cohort.obs.loc[m, "roi"]    = (zn > 0.5) & (zu > 0.5)

    counts = cohort.obs.groupby("sample")["roi"].agg(["sum","mean","count"]).rename(
        columns={"sum":"n_roi_spots","mean":"frac_roi","count":"n_total_spots"})
    counts.to_csv(ROOT / "roi_summary.csv")
    log(f"   ROI counts:\n{counts.to_string()}")

    # ROI vs non-ROI stats per section
    rows = []
    gex_genes = [g for g in TF_GENES + LIGAND_GENES if g in cohort.var_names]
    X = cohort[:, gex_genes].X
    if sp.issparse(X): X = X.toarray()
    gex = pd.DataFrame(X, index=cohort.obs_names, columns=gex_genes)

    for s in samples:
        m = cohort.obs["sample"].values == s
        sub_idx = cohort.obs_names[m]
        in_roi = cohort.obs.loc[sub_idx, "roi"].values
        if in_roi.sum() < 5:
            continue
        for col in ["MP1_score","MP2_score","MP3_score","MP4_score"] + \
                   [f"progeny_{p}" for p in PATHWAYS_PROGENY_HIGHLIGHT]:
            if col not in cohort.obs.columns: continue
            v = cohort.obs.loc[sub_idx, col].values.astype(float)
            rows.append({"sample": s, "metric": col, "type": "obs",
                         "mean_roi": float(np.mean(v[in_roi])),
                         "mean_nonroi": float(np.mean(v[~in_roi])),
                         "delta": float(np.mean(v[in_roi]) - np.mean(v[~in_roi]))})
        for ct in NEU_COLS + ["Macro_SPP1","Macro_C1QC","Fibroblast","Endothelial","Malignant","T_NK","B"]:
            if ct in abund.columns:
                v = abund.loc[sub_idx, ct].values
                rows.append({"sample": s, "metric": ct, "type": "celltype",
                             "mean_roi": float(np.mean(v[in_roi])),
                             "mean_nonroi": float(np.mean(v[~in_roi])),
                             "delta": float(np.mean(v[in_roi]) - np.mean(v[~in_roi]))})
        for g in gex_genes:
            v = gex.loc[sub_idx, g].values
            rows.append({"sample": s, "metric": g, "type": "gene",
                         "mean_roi": float(np.mean(v[in_roi])),
                         "mean_nonroi": float(np.mean(v[~in_roi])),
                         "delta": float(np.mean(v[in_roi]) - np.mean(v[~in_roi]))})
    cmp_df = pd.DataFrame(rows)
    cmp_df.to_csv(ROOT / "roi_vs_nonroi_stats.csv", index=False)
    if len(cmp_df):
        agg = cmp_df.groupby(["metric","type"])[["mean_roi","mean_nonroi","delta"]].mean().sort_values("delta", ascending=False)
        agg.to_csv(ROOT / "roi_vs_nonroi_aggregate.csv")
        log(f"\nTop +delta:\n{agg.head(20).to_string()}")
        log(f"\nBottom delta:\n{agg.tail(8).to_string()}")
    cohort.write_h5ad(str(out_h5), compression="gzip")


def stage_plots(cohort: ad.AnnData):
    log(f"[plots] per-section spatial overlays")
    samples = sorted(cohort.obs["sample"].unique().tolist())
    abund = cohort.obsm["q05_cell_abundance"].copy()
    abund.index = cohort.obs_names
    panels = ["Malignant","Macro_SPP1","Neu_Inflammatory","Neu_OSM_priming",
              "Neu_OSM_low","Fibroblast","Endothelial"]
    for s in samples:
        try:
            sec = sc.read_h5ad(str(SEC_DIR / f"{s}.h5ad"))
            sec.var_names_make_unique()
            mask = cohort.obs["sample"].values == s
            sub_idx = cohort.obs_names[mask]
            section_bc = pd.Index([n[:-len("-"+s)] if n.endswith("-"+s) else n for n in sub_idx])
            for col in panels:
                if col in abund.columns:
                    sec.obs[col] = pd.Series(abund.loc[sub_idx, col].values, index=section_bc).reindex(sec.obs_names).values
            for col in ["MP1_score","MP2_score","MP3_score","MP4_score","progeny_NFkB","progeny_TGFb","roi","neu_total"]:
                if col in cohort.obs.columns:
                    sec.obs[col] = pd.Series(cohort.obs.loc[sub_idx, col].values, index=section_bc).reindex(sec.obs_names).values

            # Page 1: cell types
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))
            axes = axes.flatten()
            for i, ct in enumerate(panels):
                if ct in sec.obs.columns:
                    sc.pl.spatial(sec, color=ct, library_id=s, ax=axes[i], show=False,
                                  cmap="magma", size=1.4, frameon=False, title=ct, colorbar_loc=None)
                else:
                    axes[i].axis("off")
            sc.pl.spatial(sec, color="total_counts", library_id=s, ax=axes[-1], show=False,
                          cmap="viridis", size=1.4, frameon=False,
                          title=f"total_counts (n={sec.n_obs})", colorbar_loc=None)
            fig.suptitle(f"{s}  cell2location q05 abundance", fontsize=13)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}_celltypes.png", dpi=130, bbox_inches="tight")
            plt.close(fig)

            # Page 2: MP + ROI + NFkB + TGFb
            fig, axes = plt.subplots(2, 3, figsize=(13, 9))
            for i, col in enumerate(["MP1_score","MP2_score","MP3_score","MP4_score","progeny_NFkB","progeny_TGFb"]):
                ax = axes[i//3, i%3]
                if col in sec.obs.columns:
                    sc.pl.spatial(sec, color=col, library_id=s, ax=ax, show=False,
                                  cmap="RdBu_r", size=1.4, frameon=False, title=col, colorbar_loc="right")
                else:
                    ax.axis("off")
            fig.suptitle(f"{s}  MP1-4 + NFkB + TGFb", fontsize=13)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}_pathways.png", dpi=130, bbox_inches="tight")
            plt.close(fig)

            # Page 3: ROI definition (NFkB / Neu / ROI)
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            sc.pl.spatial(sec, color="progeny_NFkB", library_id=s, ax=axes[0], show=False,
                          cmap="RdBu_r", size=1.4, frameon=False, title="NFkB activity", colorbar_loc="right")
            sc.pl.spatial(sec, color="neu_total", library_id=s, ax=axes[1], show=False,
                          cmap="magma", size=1.4, frameon=False, title="Neutrophil total q05", colorbar_loc="right")
            sec.obs["roi_int"] = sec.obs["roi"].astype(int).fillna(0)
            sc.pl.spatial(sec, color="roi_int", library_id=s, ax=axes[2], show=False,
                          cmap="Reds", size=1.4, frameon=False, title="ROI", colorbar_loc="right")
            fig.suptitle(f"{s}  ROI definition", fontsize=13)
            fig.tight_layout()
            fig.savefig(PLOTS / f"{s}_roi.png", dpi=130, bbox_inches="tight")
            plt.close(fig)

            del sec; gc.collect()
        except Exception as e:
            log(f"   [plot fail] {s}: {type(e).__name__}: {e}")


def main():
    samples = sorted([p.stem for p in SEC_DIR.glob("*.h5ad")])
    log(f"samples: {samples}")
    cohort_in = sc.read_h5ad(str(ROOT / "cohort.h5ad"))
    log(f"cohort_in: {cohort_in.shape}")

    cohort = stage_c2l(cohort_in, ROOT / "all_sections_c2l.h5ad")
    cohort = stage_mp(cohort, ROOT / "cohort_with_mp.h5ad")
    cohort = stage_progeny(cohort, ROOT / "cohort_with_progeny.h5ad")
    stage_commot(samples)
    stage_roi(cohort, ROOT / "cohort_with_roi.h5ad")
    stage_plots(cohort)
    log("[done] step09 pipeline")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        raise
