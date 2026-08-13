"""
Step 3: cell2location spot deconvolution per LUAD tumor section.

Inputs:
  ${DATA_ROOT}/ST/results/step01_qc/section_h5ad/<sample>.h5ad   per-section QC'ed Visium AnnData
  ${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv          gene × cell-type signatures (9078 × 28)
  ${DATA_ROOT}/ST/results/step01_qc/qc_summary.csv               for is_luad_tumor flag

Outputs:
  ${DATA_ROOT}/ST/results/step03_deconvolution/section_h5ad/<sample>.h5ad   per-section result
  ${DATA_ROOT}/ST/results/step03_deconvolution/c2l_models/<sample>/         per-section trained model
  ${DATA_ROOT}/ST/results/step03_deconvolution/elbo/<sample>.png            per-section ELBO curve
  ${DATA_ROOT}/ST/results/step03_deconvolution/run.log                      progress log

Params (matching paper):
  N_cells_per_location=30, detection_alpha=20, max_epochs=30000,
  batch_size=None (full batch), lr=0.002, GPU.

Robust to per-section failures: try/except per section, log error, continue.
"""
from __future__ import annotations
import os, sys, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QC_DIR   = Path("${DATA_ROOT}/ST/results/step01_qc")
REF_CSV  = Path("${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv")
OUT      = Path("${DATA_ROOT}/ST/results/step03_deconvolution")
H5AD_OUT = OUT / "section_h5ad"
MODEL_OUT = OUT / "c2l_models"
ELBO_OUT = OUT / "elbo"
for d in (OUT, H5AD_OUT, MODEL_OUT, ELBO_OUT):
    d.mkdir(parents=True, exist_ok=True)

LOG = OUT / "run.log"

def log(msg: str):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def main():
    # 1) Resolve LUAD tumor sections from qc_summary
    qc = pd.read_csv(QC_DIR / "qc_summary.csv")
    luad = qc[qc["is_luad_tumor"] == True].copy().reset_index(drop=True)
    sections = sorted(luad["sample"].tolist())
    log(f"LUAD tumor sections to process: {len(sections)} -> {sections}")

    # 2) Load reference signatures
    inf_aver = pd.read_csv(REF_CSV, index_col=0)
    log(f"reference signatures: {inf_aver.shape[0]} genes × {inf_aver.shape[1]} cell types")

    # 3) Imports inside main so import errors are logged
    import torch
    from cell2location.models import Cell2location
    use_gpu = torch.cuda.is_available()
    log(f"CUDA available: {use_gpu}; device={torch.cuda.get_device_name(0) if use_gpu else 'cpu'}")

    overall_t0 = time.time()
    successes, failures = [], []

    for i, sample in enumerate(sections, 1):
        out_h5 = H5AD_OUT / f"{sample}.h5ad"
        if out_h5.exists():
            log(f"[{i}/{len(sections)}] {sample}: already done, skipping")
            successes.append(sample)
            continue

        t_start = time.time()
        try:
            log(f"[{i}/{len(sections)}] {sample}: loading section ...")
            adata_vis = sc.read_h5ad(QC_DIR / "section_h5ad" / f"{sample}.h5ad")
            adata_vis.var_names_make_unique()
            log(f"    loaded {adata_vis.shape}; X={type(adata_vis.X).__name__}")

            # Ensure sparse + integer-valued raw counts
            if not sp.issparse(adata_vis.X):
                adata_vis.X = sp.csr_matrix(adata_vis.X)
            adata_vis.X = adata_vis.X.astype("float32")

            # Subset to intersection genes with reference
            common = sorted(set(adata_vis.var_names) & set(inf_aver.index))
            adata_vis = adata_vis[:, common].copy()
            sub_inf = inf_aver.loc[common].copy()
            log(f"    intersect with reference: {len(common)} genes")
            log(f"    spots: {adata_vis.n_obs}")

            # Setup AnnData (single section -> no batch_key)
            Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

            mod = Cell2location(
                adata_vis,
                cell_state_df=sub_inf,
                N_cells_per_location=30,
                detection_alpha=20,
            )

            log(f"    train: max_epochs=30000, batch_size=None, lr=0.002, GPU")
            train_kwargs = dict(max_epochs=30000, batch_size=None, train_size=1.0, lr=0.002)
            try:
                mod.train(accelerator="gpu" if use_gpu else "cpu", **train_kwargs)
            except TypeError:
                mod.train(**train_kwargs)

            # ELBO plot
            try:
                fig, ax = plt.subplots(figsize=(5, 3))
                mod.plot_history(1000)
                fig.tight_layout()
                fig.savefig(ELBO_OUT / f"{sample}.png", dpi=150)
                plt.close(fig)
            except Exception as e:
                log(f"    [warn] elbo plot failed: {type(e).__name__}: {e}")

            # Posterior
            log(f"    export_posterior ...")
            adata_vis = mod.export_posterior(
                adata_vis,
                sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs},
            )

            # Promote q05_cell_abundance_w_sf into obs columns for easy plotting
            key_q05 = "q05_cell_abundance_w_sf"
            key_means = "means_cell_abundance_w_sf"
            if key_q05 in adata_vis.obsm:
                df_q05 = adata_vis.obsm[key_q05].copy()
                # column names like "q05_cell_abundance_w_sf_<celltype>"
                df_q05.columns = [c.replace(key_q05 + "_", "") for c in df_q05.columns]
                adata_vis.obsm["q05_cell_abundance"] = df_q05
            if key_means in adata_vis.obsm:
                df_m = adata_vis.obsm[key_means].copy()
                df_m.columns = [c.replace(key_means + "_", "") for c in df_m.columns]
                adata_vis.obsm["mean_cell_abundance"] = df_m

            # Save model
            mod.save(str(MODEL_OUT / sample), overwrite=True)
            # Save AnnData
            adata_vis.write_h5ad(str(out_h5), compression="gzip")

            elapsed = (time.time() - t_start) / 60
            log(f"    {sample} done in {elapsed:.1f} min  ({i}/{len(sections)})")
            successes.append(sample)

            del mod, adata_vis
            gc.collect()
            if use_gpu:
                torch.cuda.empty_cache()
        except Exception as e:
            tb = traceback.format_exc()
            log(f"[ERROR] {sample}: {type(e).__name__}: {e}\n{tb}")
            failures.append((sample, f"{type(e).__name__}: {e}"))
            try:
                if use_gpu:
                    torch.cuda.empty_cache()
            except Exception:
                pass

    total_min = (time.time() - overall_t0) / 60
    log(f"\nALL DONE in {total_min:.1f} min  | success: {len(successes)}/{len(sections)}  fail: {len(failures)}")
    if failures:
        for s, msg in failures:
            log(f"  [fail] {s}: {msg}")

    # 3.2 Concatenate successes into a single cohort AnnData
    if successes:
        log("\n[merge] concatenating per-section results ...")
        adatas = {}
        for s in successes:
            a = sc.read_h5ad(H5AD_OUT / f"{s}.h5ad")
            # drop heavy uns to keep merge slim, but keep scalefactors for plotting later
            if "spatial" in a.uns and isinstance(a.uns["spatial"], dict):
                for k, v in list(a.uns["spatial"].items()):
                    a.uns["spatial"][k] = {
                        kk: v[kk] for kk in ("scalefactors",) if kk in v
                    }
            adatas[s] = a
        cohort = ad.concat(adatas, label="sample", join="outer", index_unique="-",
                           merge="unique", uns_merge="unique")
        cohort.obs_names_make_unique()
        out_cohort = OUT / "all_sections_c2l.h5ad"
        cohort.write_h5ad(str(out_cohort), compression="gzip")
        log(f"[merge] cohort: {cohort.n_obs} spots × {cohort.n_vars} genes -> {out_cohort}")
        log(f"[merge] file size: {out_cohort.stat().st_size/1e9:.2f} GB")


if __name__ == "__main__":
    main()
