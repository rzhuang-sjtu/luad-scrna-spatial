"""
Step 5d (part 1): K selection via error-curve elbow (Kotliar 2019 style).

For each patient:
    - Read cnmf_out/<key>/<key>.k_selection_stats.df.npz
      (already contains one (silhouette, prediction_error) row per K)
    - Pick elbow K by max distance from the chord connecting
      (K_min, err_min-normalised) to (K_max, err_max-normalised).
    - Constraints:
        * K in [K_MIN_GLOBAL, K_MAX_GLOBAL] (from the cNMF run)
        * K <= n_cells // MIN_CELLS_PER_PROGRAM
          (i.e. each program should represent on average >=10 cells;
          the spec said ">=" but that is mathematically impossible for
          our K=8-20 window on mid-size patients, so we use the sane
          interpretation that matches the "10 cells per program" note.)
        * If silhouette at elbow K < STAB_MIN, drop K by one step (but
          stay within global bounds).

Outputs (no consensus yet):
    results/step5d_k_selection_summary.csv
    results/step5d_k_error_curves.pdf   # all error curves + red dots
    results/step5d_k_histogram.pdf      # selected-K distribution
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path("${PROJECT_ROOT}")
CNMF_INPUT = ROOT / "data" / "cnmf_input"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

OUT_CSV = RESULTS / "step5d_k_selection_summary.csv"
OUT_CURVES_PDF = RESULTS / "step5d_k_error_curves.pdf"
OUT_HIST_PDF = RESULTS / "step5d_k_histogram.pdf"

K_MIN_GLOBAL = 8
K_MAX_GLOBAL = 20
MIN_CELLS_PER_PROGRAM = 10
STAB_MIN = 0.3


def load_patient_stats(patient_dir: Path):
    if not (patient_dir / "DONE").exists():
        return None
    key = patient_dir.name
    stats_path = patient_dir / "cnmf_out" / key / f"{key}.k_selection_stats.df.npz"
    if not stats_path.exists():
        return None
    with np.load(stats_path, allow_pickle=True) as f:
        df = pd.DataFrame(**f)
    df["k"] = df["k"].astype(int)
    df = df.sort_values("k").reset_index(drop=True)
    return key, df


def load_n_cells(patient_dir: Path) -> int:
    """Cheap: use h5py to get obs length from counts.h5ad without loading X."""
    import h5py
    with h5py.File(patient_dir / "counts.h5ad", "r") as h:
        # anndata stores obs as a group; len = any column length
        obs = h["obs"]
        if "_index" in obs.attrs:
            idx_key = obs.attrs["_index"]
            return int(obs[idx_key].shape[0])
        for k in obs.keys():
            ds = obs[k]
            if hasattr(ds, "shape") and len(ds.shape) == 1:
                return int(ds.shape[0])
    raise RuntimeError("Couldn't infer n_cells")


def find_elbow_k(k_values: np.ndarray, errors: np.ndarray):
    """Max perpendicular distance from the chord p1->p2 (Kneedle-style)."""
    k_values = np.asarray(k_values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    if len(k_values) < 2:
        return int(k_values[0]), 0.0
    k_rng = k_values.max() - k_values.min()
    e_rng = errors.max() - errors.min()
    if k_rng == 0 or e_rng == 0:
        return int(k_values[int(np.argmin(errors))]), 0.0
    k_norm = (k_values - k_values.min()) / k_rng
    e_norm = (errors - errors.min()) / e_rng
    p1x, p1y = k_norm[0], e_norm[0]
    p2x, p2y = k_norm[-1], e_norm[-1]
    dx, dy = p2x - p1x, p2y - p1y
    line_len = np.hypot(dx, dy)
    # 2D cross product magnitude: |dx*(p1y-py) - dy*(p1x-px)| / ||p2-p1||
    dists = np.abs(dx * (p1y - e_norm) - dy * (p1x - k_norm)) / line_len
    elbow_idx = int(np.argmax(dists))
    return int(k_values[elbow_idx]), float(dists[elbow_idx])


def pick_k(df: pd.DataFrame, n_cells: int):
    k_all = df["k"].values.astype(int)
    err = df["prediction_error"].values.astype(float)
    stab = df["silhouette"].values.astype(float)

    k_max_cells = max(K_MIN_GLOBAL, n_cells // MIN_CELLS_PER_PROGRAM)
    k_max = min(K_MAX_GLOBAL, k_max_cells)
    mask = (k_all >= K_MIN_GLOBAL) & (k_all <= k_max)
    if mask.sum() < 2:
        mask = np.ones_like(k_all, dtype=bool)
        k_max_note = "FALLBACK_full_range"
    else:
        k_max_note = f"k_max={k_max}"

    kv = k_all[mask]
    ev = err[mask]
    sv = stab[mask]

    elbow_k, dist = find_elbow_k(kv, ev)
    idx = int(np.where(kv == elbow_k)[0][0])
    stab_at = float(sv[idx])
    err_at = float(ev[idx])

    adj_reason = "elbow"
    if stab_at < STAB_MIN:
        # Drop one step but keep within global bound
        new_k = max(elbow_k - 1, int(kv.min()))
        if new_k in kv and new_k != elbow_k:
            new_idx = int(np.where(kv == new_k)[0][0])
            elbow_k = new_k
            stab_at = float(sv[new_idx])
            err_at = float(ev[new_idx])
            adj_reason = "elbow_minus1_low_stab"

    return {
        "selected_k": int(elbow_k),
        "stability_at_k": stab_at,
        "error_at_k": err_at,
        "elbow_distance": float(dist),
        "k_max_constraint": k_max_note,
        "reason": adj_reason,
    }


def main():
    patient_dirs = sorted(p for p in CNMF_INPUT.iterdir() if p.is_dir())
    print(f"Scanning {len(patient_dirs)} patient dirs")

    records = []
    curves = {}  # key -> df
    selected = {}  # key -> elbow_k
    skipped = 0

    for pdir in patient_dirs:
        res = load_patient_stats(pdir)
        if res is None:
            skipped += 1
            continue
        key, df = res
        n_cells = load_n_cells(pdir)
        pick = pick_k(df, n_cells)

        curves[key] = df
        selected[key] = pick["selected_k"]
        records.append({
            "patient_key": key,
            "n_cells": n_cells,
            **pick,
        })

    print(f"Skipped {skipped}; loaded {len(records)}")

    summary = (
        pd.DataFrame(records)
        .sort_values("patient_key")
        .reset_index(drop=True)
    )
    summary.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")

    print()
    print("=" * 95)
    print("Per-patient elbow K:")
    print("=" * 95)
    cols = ["patient_key", "n_cells", "selected_k",
            "stability_at_k", "error_at_k",
            "k_max_constraint", "reason"]
    print(summary[cols].to_string(index=False))
    print()
    print("Selected-K distribution:")
    print(summary["selected_k"].value_counts().sort_index())
    print()
    print(f"Median selected K: {summary['selected_k'].median()}")
    print(f"Mean selected K: {summary['selected_k'].mean():.2f}")
    print(f"Mean stability at selected K: {summary['stability_at_k'].mean():.3f}")
    low_stab = summary[summary["stability_at_k"] < STAB_MIN]
    print(f"Patients with stability < {STAB_MIN}: {len(low_stab)}")
    if len(low_stab):
        print(low_stab[["patient_key", "selected_k", "stability_at_k"]].to_string(index=False))

    # Figure 1: all error curves with elbow points
    fig, ax = plt.subplots(figsize=(10, 6))
    for key, df in curves.items():
        e = df["prediction_error"].values
        e_norm = (e - e.min()) / (e.max() - e.min()) if e.max() > e.min() else e * 0
        ax.plot(df["k"], e_norm, alpha=0.25, linewidth=0.8, color="grey")
        sk = selected[key]
        row = df[df["k"] == sk].iloc[0]
        e_val = row["prediction_error"]
        e_val_norm = (e_val - e.min()) / (e.max() - e.min()) if e.max() > e.min() else 0.0
        ax.scatter([sk], [e_val_norm], color="red", s=15, zorder=5)
    ax.set_xlabel("K")
    ax.set_ylabel("Normalised prediction error (per-patient min-max)")
    ax.set_title(
        f"cNMF prediction-error curves ({len(curves)} patients) "
        f"— red dots = elbow"
    )
    plt.tight_layout()
    plt.savefig(OUT_CURVES_PDF, dpi=150)
    plt.close()
    print(f"Wrote {OUT_CURVES_PDF}")

    # Figure 2: histogram
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.arange(K_MIN_GLOBAL, K_MAX_GLOBAL + 2) - 0.5
    ax.hist(summary["selected_k"], bins=bins, edgecolor="black")
    ax.axvline(summary["selected_k"].median(), color="red", linestyle="--",
               label=f"median K={summary['selected_k'].median():.0f}")
    ax.axvline(15, color="grey", linestyle=":", label="HCC paper K=15")
    ax.set_xticks(range(K_MIN_GLOBAL, K_MAX_GLOBAL + 1))
    ax.set_xlabel("Elbow K")
    ax.set_ylabel("# patients")
    ax.set_title("Distribution of elbow K per patient")
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUT_HIST_PDF, dpi=150)
    plt.close()
    print(f"Wrote {OUT_HIST_PDF}")


if __name__ == "__main__":
    main()
