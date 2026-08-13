"""
Step 3 (joint version): cell2location spot deconvolution on the merged 12-section
LUAD tumor cohort in a single training run, with batch_key='sample'.

This is more efficient than per-section runs because:
  - the model amortizes inference across spots,
  - 28k spots × 9078 genes still fits comfortably on a 20 GB GPU.

Inputs:
  ${DATA_ROOT}/ST/results/step01_qc/section_h5ad/<sample>.h5ad   per-section QC'ed h5ad
  ${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv          gene × cell-type signatures (9078 × 28)

Outputs:
  ${DATA_ROOT}/ST/results/step03_deconvolution/all_sections_c2l.h5ad   joint result AnnData
  ${DATA_ROOT}/ST/results/step03_deconvolution/c2l_model/              trained model dir
  ${DATA_ROOT}/ST/results/step03_deconvolution/elbo.png                training curve
  ${DATA_ROOT}/ST/results/step03_deconvolution/run.log                 progress log

Params (matching paper):
  N_cells_per_location=30, detection_alpha=20, max_epochs=30000,
  batch_size=None (full batch by default), lr=0.002, GPU.

Fallback: if OOM during train or export_posterior, retry with batch_size=2500.
"""
from __future__ import annotations
import os, sys, time, gc, traceback
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QC_DIR  = Path("${DATA_ROOT}/ST/results/step01_qc")
REF_CSV = Path("${DATA_ROOT}/ST/results/step02_reference/inf_aver.csv")
OUT     = Path("${DATA_ROOT}/ST/results/step03_deconvolution")
MODEL_OUT = OUT / "c2l_model"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_OUT.mkdir(parents=True, exist_ok=True)

LOG = OUT / "run.log"
def log(msg: str):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def main():
    qc = pd.read_csv(QC_DIR / "qc_summary.csv")
    luad = qc[qc["is_luad_tumor"] == True].copy()
    sections = sorted(luad["sample"].tolist())
    log(f"LUAD tumor sections: {len(sections)} -> {sections}")

    log("loading per-section h5ads ...")
    adatas = {}
    for s in sections:
        a = sc.read_h5ad(QC_DIR / "section_h5ad" / f"{s}.h5ad")
        a.var_names_make_unique()
        # keep only scalefactors in uns to avoid dragging huge images into the joint h5ad
        if "spatial" in a.uns and isinstance(a.uns["spatial"], dict):
            for k, v in list(a.uns["spatial"].items()):
                a.uns["spatial"][k] = {kk: v[kk] for kk in ("scalefactors",) if kk in v}
        adatas[s] = a

    cohort = ad.concat(adatas, label="sample", join="outer", index_unique="-",
                       merge="unique", uns_merge="unique")
    cohort.obs_names_make_unique()
    log(f"cohort merged: {cohort.shape}; samples: {cohort.obs['sample'].nunique()}")

    # Reference
    inf_aver = pd.read_csv(REF_CSV, index_col=0)
    log(f"reference signatures: {inf_aver.shape[0]} genes × {inf_aver.shape[1]} cell types")

    common = sorted(set(cohort.var_names) & set(inf_aver.index))
    cohort = cohort[:, common].copy()
    sub_inf = inf_aver.loc[common].copy()
    log(f"intersect with reference: {len(common)} genes; cohort -> {cohort.shape}")

    # Ensure sparse + integer
    if not sp.issparse(cohort.X):
        cohort.X = sp.csr_matrix(cohort.X)
    cohort.X = cohort.X.astype("float32")

    # Sanity: integer values
    sample_data = cohort.X[: min(2000, cohort.n_obs)].data if sp.issparse(cohort.X) else cohort.X[: min(2000, cohort.n_obs)]
    if hasattr(sample_data, "data"): sample_data = sample_data.data
    is_int = bool(np.all(np.equal(np.mod(sample_data[sample_data!=0], 1), 0)))
    log(f"X integer? {is_int}")

    import torch
    from cell2location.models import Cell2location
    use_gpu = torch.cuda.is_available()
    log(f"CUDA available: {use_gpu}; device={torch.cuda.get_device_name(0) if use_gpu else 'cpu'}")

    Cell2location.setup_anndata(adata=cohort, batch_key="sample")
    mod = Cell2location(
        cohort,
        cell_state_df=sub_inf,
        N_cells_per_location=30,
        detection_alpha=20,
    )
    log("Cell2location model built; cell_state_df rows={}, cols={}".format(*sub_inf.shape))

    # Train (full batch first; fallback to 2500 on OOM)
    train_kwargs_full = dict(max_epochs=30000, batch_size=None, train_size=1.0, lr=0.002)
    train_kwargs_mini = dict(max_epochs=30000, batch_size=2500, train_size=1.0, lr=0.002)
    accel = "gpu" if use_gpu else "cpu"

    t_start = time.time()
    try:
        log(f"train: max_epochs=30000, batch_size=None (full batch), lr=0.002, accel={accel}")
        mod.train(accelerator=accel, **train_kwargs_full)
    except torch.cuda.OutOfMemoryError as e:
        log(f"[oom] full-batch failed: {e}. retrying with batch_size=2500 ...")
        torch.cuda.empty_cache()
        gc.collect()
        # rebuild model since training state may be partial
        Cell2location.setup_anndata(adata=cohort, batch_key="sample")
        mod = Cell2location(cohort, cell_state_df=sub_inf, N_cells_per_location=30, detection_alpha=20)
        mod.train(accelerator=accel, **train_kwargs_mini)
    except TypeError:
        # very old API
        log("[fallback] train() did not accept accelerator; retrying without it")
        mod.train(**train_kwargs_full)

    log(f"training elapsed: {(time.time()-t_start)/60:.1f} min")

    # ELBO
    try:
        fig, ax = plt.subplots(figsize=(5, 3))
        mod.plot_history(1000)
        fig.tight_layout()
        fig.savefig(OUT / "elbo.png", dpi=150)
        plt.close(fig)
        log("elbo plot saved")
    except Exception as e:
        log(f"[warn] elbo plot failed: {type(e).__name__}: {e}")

    # Posterior (with OOM fallback on batch_size)
    log("export_posterior: num_samples=1000")
    t_p = time.time()
    try:
        cohort = mod.export_posterior(
            cohort,
            sample_kwargs={"num_samples": 1000, "batch_size": cohort.n_obs},
        )
    except torch.cuda.OutOfMemoryError as e:
        log(f"[oom] export_posterior full failed: {e}. retrying batch_size=2500 ...")
        torch.cuda.empty_cache(); gc.collect()
        cohort = mod.export_posterior(
            cohort,
            sample_kwargs={"num_samples": 1000, "batch_size": 2500},
        )
    log(f"posterior elapsed: {(time.time()-t_p)/60:.1f} min")

    # Promote q05/means cell abundance into easy-to-access obsm matrices
    for src_key, dst_key in [("q05_cell_abundance_w_sf", "q05_cell_abundance"),
                             ("means_cell_abundance_w_sf", "mean_cell_abundance")]:
        if src_key in cohort.obsm:
            df = cohort.obsm[src_key].copy()
            df.columns = [c.replace(src_key + "_", "") for c in df.columns]
            cohort.obsm[dst_key] = df

    # Save model + h5ad
    mod.save(str(MODEL_OUT), overwrite=True)
    out_h5 = OUT / "all_sections_c2l.h5ad"
    cohort.write_h5ad(str(out_h5), compression="gzip")
    log(f"saved: {out_h5}  ({out_h5.stat().st_size/1e9:.2f} GB)")
    log(f"cohort obs columns: {list(cohort.obs.columns)[-10:]}")
    log(f"cohort obsm keys: {list(cohort.obsm.keys())}")

    # Quick per-sample summary of mean abundances for sanity
    if "q05_cell_abundance" in cohort.obsm:
        df = cohort.obsm["q05_cell_abundance"].copy()
        df["sample"] = cohort.obs["sample"].values
        means = df.groupby("sample").mean()
        means.to_csv(OUT / "per_sample_mean_q05_abundance.csv")
        log(f"per-sample mean q05 abundance saved -> per_sample_mean_q05_abundance.csv  shape={means.shape}")
    log("[done]")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[FATAL] {type(e).__name__}: {e}\n{traceback.format_exc()}")
        sys.exit(1)
