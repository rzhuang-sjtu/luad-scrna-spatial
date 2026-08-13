"""
Export Figure S2 (cNMF QC) data to ${WORK_ROOT}/luad_figures/fig_s2/.

Files produced:
    s2_k_selection_summary.csv      — per-patient K / stability / error
    s2_k_error_curves.csv           — long: patient_key x K x (err, stab)
    s2_k_histogram.csv              — K x count
    s2_stability_vs_ncells.csv      — per-patient scatter data
    s2_gep_pool_info.csv            — summary of the gep pool (counts)
    s2_gep_top30_genes.csv          — top-30 gene list per GEP (long)
    s2_consensus_example_<key>.csv  — pairwise Pearson correlation of
                                       all per-iter spectra, subsampled
                                       to SUBSAMPLE_ITERS iterations
                                       (for 3 representative patients)
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("${PROJECT_ROOT}")
CNMF_INPUT = ROOT / "data" / "cnmf_input"
CNMF_OUTPUT = ROOT / "data" / "cnmf_output"
RESULTS = ROOT / "results"

OUT_DIR = Path("${WORK_ROOT}/luad_figures/fig_s2")
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_SUMMARY_CSV = RESULTS / "step5d_k_selection_summary.csv"
GEP_META_CSV = CNMF_OUTPUT / "gep_metadata.csv"
GEP_TOP100_CSV = CNMF_OUTPUT / "gep_top100_genes.csv"

SUBSAMPLE_ITERS = 50  # per-iter correlation matrix size = K * SUBSAMPLE_ITERS


def log(msg):
    print(msg, flush=True)


def derive_dataset(key: str) -> str:
    return key.split("__")[0]


def build_summary():
    df = pd.read_csv(K_SUMMARY_CSV)
    df["dataset"] = df["patient_key"].apply(derive_dataset)
    cols = [
        "patient_key", "dataset", "n_cells", "selected_k",
        "stability_at_k", "error_at_k", "elbow_distance",
        "k_max_constraint", "reason",
    ]
    out = df[cols].rename(columns={
        "stability_at_k": "stability",
        "error_at_k": "prediction_error",
    })
    p = OUT_DIR / "s2_k_selection_summary.csv"
    out.to_csv(p, index=False)
    log(f"  {p.name}  ({len(out)} rows)")
    return df


def build_k_error_curves():
    """Long format: (patient_key, K, prediction_error, silhouette)."""
    rows = []
    for pdir in sorted(CNMF_INPUT.iterdir()):
        if not pdir.is_dir():
            continue
        key = pdir.name
        if not (pdir / "DONE").exists():
            continue
        stats_path = pdir / "cnmf_out" / key / f"{key}.k_selection_stats.df.npz"
        if not stats_path.exists():
            continue
        with np.load(stats_path, allow_pickle=True) as z:
            df = pd.DataFrame(**z)
        df["k"] = df["k"].astype(int)
        df = df.rename(columns={
            "silhouette": "stability",
            "prediction_error": "prediction_error",
        })
        df["patient_key"] = key
        df["dataset"] = derive_dataset(key)
        rows.append(df[["patient_key", "dataset", "k", "stability",
                        "prediction_error", "local_density_threshold"]])
    long = pd.concat(rows, ignore_index=True)
    long = long.rename(columns={"k": "K"})
    p = OUT_DIR / "s2_k_error_curves.csv"
    long.to_csv(p, index=False)
    log(f"  {p.name}  ({len(long)} rows, {long['patient_key'].nunique()} patients)")


def build_histogram(summary_df: pd.DataFrame):
    hist = (
        summary_df["selected_k"]
        .value_counts()
        .rename_axis("K")
        .reset_index(name="count")
        .sort_values("K")
        .reset_index(drop=True)
    )
    p = OUT_DIR / "s2_k_histogram.csv"
    hist.to_csv(p, index=False)
    log(f"  {p.name}")


def build_stability_scatter(summary_df: pd.DataFrame):
    out = summary_df[
        ["patient_key", "dataset", "n_cells", "selected_k", "stability_at_k"]
    ].rename(columns={"stability_at_k": "stability"})
    p = OUT_DIR / "s2_stability_vs_ncells.csv"
    out.to_csv(p, index=False)
    log(f"  {p.name}")


def build_gep_pool_info():
    meta = pd.read_csv(GEP_META_CSV)
    meta["dataset"] = meta["patient_key"].apply(derive_dataset)

    summary_rows = [
        {"metric": "n_GEPs_total", "value": len(meta)},
        {"metric": "n_patients", "value": meta["patient_key"].nunique()},
        {"metric": "n_datasets", "value": meta["dataset"].nunique()},
        {"metric": "mean_K_per_patient", "value": meta["k_total"].mean()},
        {"metric": "median_K_per_patient", "value": meta["k_total"].median()},
        {"metric": "min_K_per_patient", "value": meta["k_total"].min()},
        {"metric": "max_K_per_patient", "value": meta["k_total"].max()},
        {"metric": "mean_stability", "value": meta["stability"].mean()},
    ]
    for ds, sub in meta.groupby("dataset"):
        summary_rows.append({
            "metric": f"n_GEPs_{ds}", "value": int(len(sub)),
        })
    info = pd.DataFrame(summary_rows)
    p = OUT_DIR / "s2_gep_pool_info.csv"
    info.to_csv(p, index=False)
    log(f"  {p.name}")


def build_gep_top30():
    top100 = pd.read_csv(GEP_TOP100_CSV)
    top30 = top100[top100["rank"] <= 30].copy()
    top30["dataset"] = top30["gep_id"].apply(
        lambda g: g.split("__")[0] if "__" in g else ""
    )
    p = OUT_DIR / "s2_gep_top30_genes.csv"
    top30.to_csv(p, index=False)
    log(f"  {p.name}  ({len(top30)} rows, "
        f"{top30['gep_id'].nunique()} GEPs)")


def build_consensus_example(key: str, rng: np.random.Generator):
    """Pairwise Pearson correlation of per-iter spectra for one patient."""
    pdir = CNMF_INPUT / key
    if not (pdir / "DONE").exists():
        log(f"  SKIP consensus example {key}: no DONE flag")
        return
    summary = pd.read_csv(K_SUMMARY_CSV)
    k_row = summary[summary["patient_key"] == key]
    if len(k_row) == 0:
        log(f"  SKIP {key}: not in summary")
        return
    k = int(k_row.iloc[0]["selected_k"])

    tmp = pdir / "cnmf_out" / key / "cnmf_tmp"
    pat = re.compile(rf"{re.escape(key)}\.spectra\.k_{k}\.iter_(\d+)\.df\.npz$")
    iter_files = {}
    for f in tmp.iterdir():
        m = pat.match(f.name)
        if m:
            iter_files[int(m.group(1))] = f
    if not iter_files:
        log(f"  WARN {key}: no iter files found for K={k}")
        return

    iter_ids = sorted(iter_files)
    if len(iter_ids) > SUBSAMPLE_ITERS:
        chosen = sorted(rng.choice(iter_ids, SUBSAMPLE_ITERS, replace=False))
    else:
        chosen = iter_ids

    spectra_list = []
    labels = []
    for it in chosen:
        with np.load(iter_files[it], allow_pickle=True) as z:
            df = pd.DataFrame(**z)  # K rows x genes cols
        spectra_list.append(df.values)
        for p in df.index:
            labels.append(f"iter{it}_P{int(p)}")
    spectra = np.vstack(spectra_list)  # (K*len(chosen)) x genes

    # Pearson correlation between spectra rows
    corr = np.corrcoef(spectra)
    corr_df = pd.DataFrame(corr, index=labels, columns=labels)
    p = OUT_DIR / f"s2_consensus_example_{key}.csv"
    corr_df.round(4).to_csv(p)
    log(f"  {p.name}  (K={k}, iters={len(chosen)}, "
        f"{corr_df.shape[0]}x{corr_df.shape[1]}, "
        f"{p.stat().st_size/1024/1024:.1f} MB)")


def pick_representative_patients(summary_df: pd.DataFrame):
    def closest(target_n):
        diff = (summary_df["n_cells"] - target_n).abs()
        return summary_df.iloc[int(diff.idxmin())]["patient_key"]
    small = summary_df[summary_df["n_cells"] < 200].sort_values(
        "n_cells", ascending=False).iloc[0]["patient_key"]
    mid = closest(400)
    large = summary_df[summary_df["n_cells"] >= 1000].sort_values(
        "n_cells").iloc[0]["patient_key"]
    return small, mid, large


def main():
    log(f"Output dir: {OUT_DIR}")
    log("Exporting Fig S2 data ...")

    summary_df = build_summary()
    build_k_error_curves()
    build_histogram(summary_df)
    build_stability_scatter(summary_df)
    build_gep_pool_info()
    build_gep_top30()

    rng = np.random.default_rng(42)
    small, mid, large = pick_representative_patients(summary_df)
    log(f"Representative patients for consensus example:")
    log(f"  small: {small}")
    log(f"  mid  : {mid}")
    log(f"  large: {large}")
    for key in (small, mid, large):
        build_consensus_example(key, rng)

    log("")
    log("File list:")
    for f in sorted(OUT_DIR.iterdir()):
        size_mb = f.stat().st_size / 1024 / 1024
        log(f"  {f.name:55s}  {size_mb:6.2f} MB")


if __name__ == "__main__":
    main()
