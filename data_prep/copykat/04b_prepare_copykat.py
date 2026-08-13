#!/usr/bin/env python3
"""
Step 4b: Prepare per-patient CopyKAT inputs.

- Per (dataset, patient_id) directory: copykat_input/{dataset}__{patient_id}/
- Outputs: counts.mtx (genes × cells, int), barcodes.tsv, genes.tsv, ref_cells.txt, metadata.json
- Reference strategy:
    paired        : patient has Normal_Lung / Adjacent_Normal → use those
    dataset_pool  : no normal in patient, but other patients in the same dataset have → sample 300 from pool
    ref_free      : no normal in the whole dataset
  (Normal_LN is not used as ref)
- Downsample: if a patient has > 15000 cells, randomly downsample to 15000 (keep own normal), random_state=42
- Write run_order.txt (ascending by cell count)
"""
import gc
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite

IN_H5AD = Path("${PROJECT_ROOT}/data/processed/luad_integrated.h5ad")
OUT_DIR = Path("${PROJECT_ROOT}/data/copykat_input")
RESULTS_DIR = Path("${PROJECT_ROOT}/results")
STATS_CSV = RESULTS_DIR / "step4_patient_stats.csv"
RUN_ORDER = OUT_DIR / "run_order.txt"

MAX_CELLS_PER_PATIENT = 15000
POOL_REF_N = 300          # number of normal cells sampled from other patients in dataset_pool mode
NORMAL_CAP_PAIRED = 5000  # normal cap in paired mode
NORMAL_FLOOR = 200        # keep all normals if below this count (avoid too few refs)
MIN_CELLS = 50            # minimum cells per patient (CopyKAT requirement)
RAND_STATE = 42

NORMAL_TT = {"Normal_Lung", "Adjacent_Normal"}   # exclude Normal_LN

OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def safe_name(s):
    """Sanitize patient_id for use as a directory name"""
    return str(s).replace("/", "_").replace(" ", "_")


def main():
    t0 = datetime.now()
    log("=" * 60)
    log(f"Step 4b Prepare started: {t0}")
    log("=" * 60)

    log(f"Reading {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    log(f"  shape: {adata.shape}")

    # Ensure X is raw counts (int)
    X = adata.X if not sparse.issparse(adata.X) else adata.X.tocsr()
    assert X.dtype in (np.float32, np.float64, np.int32, np.int64), f"Unexpected X dtype: {X.dtype}"

    obs = adata.obs.copy()
    obs["cell_idx"] = np.arange(adata.n_obs)
    # Normalize fields
    obs["dataset"] = obs["dataset"].astype(str)
    obs["patient_id"] = obs["patient_id"].astype(str)
    obs["tissue_type"] = obs["tissue_type"].astype(str)

    # Drop empty patient_id
    valid_mask = (obs["patient_id"] != "") & (obs["patient_id"].str.lower() != "nan")
    n_drop = (~valid_mask).sum()
    if n_drop > 0:
        log(f"Dropped {n_drop:,} cells without patient_id")
    obs = obs[valid_mask].copy()
    obs["barcode"] = obs.index.astype(str)
    obs = obs.reset_index(drop=True)

    # Build (dataset, patient) index
    obs["dp_key"] = obs["dataset"] + "__" + obs["patient_id"].map(safe_name)
    patients = sorted(obs["dp_key"].unique())
    log(f"Unique (dataset, patient) pairs: {len(patients)}")

    # --- Normal pool stats by dataset ---
    log("\nNormal cell pool by dataset (Normal_Lung + Adjacent_Normal):")
    dataset_normal_pool = {}     # dataset -> obs rows (DataFrame)
    for ds in obs["dataset"].unique():
        sub = obs[(obs["dataset"] == ds) & (obs["tissue_type"].isin(NORMAL_TT))]
        dataset_normal_pool[ds] = sub
        log(f"{ds}: normal cells {len(sub):,}")

    # --- Per-patient stats and decisions ---
    stats_rows = []
    written = 0
    skipped = []

    rng = np.random.default_rng(RAND_STATE)

    for dp in patients:
        sub = obs[obs["dp_key"] == dp]
        ds = sub["dataset"].iloc[0]
        pid = sub["patient_id"].iloc[0]
        n_total = len(sub)

        tt_counts = sub["tissue_type"].value_counts().to_dict()
        n_own_normal = sum(v for k, v in tt_counts.items() if k in NORMAL_TT)

        # Choose reference mode
        ds_pool = dataset_normal_pool[ds]
        ds_pool_other = ds_pool[ds_pool["dp_key"] != dp]   # normals from other patients in the same dataset

        if n_own_normal >= 10:
            ref_mode = "paired"
        elif len(ds_pool_other) >= 50:
            ref_mode = "dataset_pool"
        else:
            ref_mode = "ref_free"

        # ---- Assemble input cell set for this patient ----
        if n_total < MIN_CELLS and ref_mode == "ref_free":
            # Too few cells; skip
            skipped.append({"dp_key": dp, "reason": f"n_cells={n_total} < {MIN_CELLS} (ref_free)"})
            stats_rows.append({
                "dp_key": dp, "dataset": ds, "patient_id": pid,
                "n_cells_total": n_total, "n_own_normal": n_own_normal,
                "ref_mode": "SKIPPED", "n_ref_cells": 0,
                "n_cells_output": 0, "tissue_types": str(tt_counts),
            })
            continue

        own = sub
        if ref_mode == "paired":
            own_normal = own[own["tissue_type"].isin(NORMAL_TT)]
            own_non_normal = own[~own["tissue_type"].isin(NORMAL_TT)]

            # Normal downsample: keep all if < NORMAL_FLOOR; else cap at NORMAL_CAP_PAIRED
            n_norm_avail = len(own_normal)
            if n_norm_avail < NORMAL_FLOOR:
                n_norm_keep = n_norm_avail
            else:
                n_norm_keep = min(n_norm_avail, NORMAL_CAP_PAIRED)
            if n_norm_keep < n_norm_avail:
                norm_idx = rng.choice(own_normal["cell_idx"].values,
                                      size=n_norm_keep, replace=False)
                own_normal = own_normal[own_normal["cell_idx"].isin(norm_idx)]

            # Fill remaining budget with tumour cells
            tumor_budget = MAX_CELLS_PER_PATIENT - n_norm_keep
            if len(own_non_normal) > tumor_budget and tumor_budget > 0:
                idx = rng.choice(own_non_normal["cell_idx"].values,
                                 size=tumor_budget, replace=False)
                own_non_normal = own_non_normal[own_non_normal["cell_idx"].isin(idx)]
            out_cells = pd.concat([own_normal, own_non_normal])
            ref_ids = own_normal["barcode"].tolist()

        elif ref_mode == "dataset_pool":
            # Sample normals from other patients in the same dataset
            n_pool = min(POOL_REF_N, len(ds_pool_other))
            pool_idx = rng.choice(ds_pool_other.index.values, size=n_pool, replace=False)
            pool = ds_pool_other.loc[pool_idx]
            budget = MAX_CELLS_PER_PATIENT - n_pool
            if len(own) > budget:
                idx = rng.choice(own["cell_idx"].values, size=budget, replace=False)
                own = own[own["cell_idx"].isin(idx)]
            out_cells = pd.concat([own, pool])
            ref_ids = pool["barcode"].tolist()

        else:  # ref_free
            if len(own) > MAX_CELLS_PER_PATIENT:
                idx = rng.choice(own["cell_idx"].values, size=MAX_CELLS_PER_PATIENT, replace=False)
                own = own[own["cell_idx"].isin(idx)]
            out_cells = own
            ref_ids = []

        # ---- Export ----
        pdir = OUT_DIR / dp
        pdir.mkdir(parents=True, exist_ok=True)

        cell_idx = out_cells["cell_idx"].values
        barcodes = out_cells["barcode"].tolist()
        # Subset X by cell_idx → genes × cells
        sub_X = X[cell_idx, :].T.tocoo()  # (genes, cells) COO
        # Cast to int (CopyKAT requires integer counts)
        sub_X = sparse.coo_matrix(
            (sub_X.data.astype(np.int32), (sub_X.row, sub_X.col)),
            shape=sub_X.shape,
        )

        mmwrite(str(pdir / "counts.mtx"), sub_X, field="integer")
        with open(pdir / "barcodes.tsv", "w") as f:
            f.write("\n".join(barcodes) + "\n")
        with open(pdir / "genes.tsv", "w") as f:
            f.write("\n".join(adata.var_names) + "\n")
        with open(pdir / "ref_cells.txt", "w") as f:
            if ref_ids:
                f.write("\n".join(ref_ids) + "\n")
            # Empty file = ref_free

        meta = {
            "dp_key": dp,
            "dataset": ds,
            "patient_id": pid,
            "n_cells_input": int(len(out_cells)),
            "n_cells_own_total": int(n_total),
            "n_cells_own_normal": int(n_own_normal),
            "n_ref_cells": int(len(ref_ids)),
            "reference_mode": ref_mode,
            "tissue_types": tt_counts,
        }
        (pdir / "metadata.json").write_text(json.dumps(meta, indent=2))

        stats_rows.append({
            "dp_key": dp, "dataset": ds, "patient_id": pid,
            "n_cells_total": n_total, "n_own_normal": n_own_normal,
            "ref_mode": ref_mode, "n_ref_cells": len(ref_ids),
            "n_cells_output": len(out_cells),
            "tissue_types": str(tt_counts),
        })
        written += 1
        if written % 10 == 0 or written == len(patients):
            log(f"Export progress: {written}/{len(patients)}")

    # ---- Write stats and run_order ----
    stats_df = pd.DataFrame(stats_rows).sort_values("n_cells_output")
    stats_df.to_csv(STATS_CSV, index=False)
    log(f"\n Stats table: {STATS_CSV}")

    runnable = stats_df[stats_df["ref_mode"] != "SKIPPED"].sort_values("n_cells_output")
    with open(RUN_ORDER, "w") as f:
        f.write("# dp_key\tn_cells\tref_mode\n")
        for _, r in runnable.iterrows():
            f.write(f"{r['dp_key']}\t{r['n_cells_output']}\t{r['ref_mode']}\n")
    log(f"Run order: {RUN_ORDER} ({len(runnable)} patients)")

    # ---- Summary ----
    log("\n" + "=" * 60)
    log("Summary")
    log("=" * 60)
    log(f"Total patients: {len(patients)}")
    log(f"Exported: {written}")
    log(f"Skipped: {len(skipped)}")
    for s in skipped:
        log(f"  - {s['dp_key']}: {s['reason']}")
    log("\nreference_mode distribution:")
    log(stats_df["ref_mode"].value_counts().to_string())
    log("\nn_cells_output quantiles:")
    log(stats_df["n_cells_output"].describe().to_string())

    elapsed = (datetime.now() - t0).total_seconds() / 60
    log(f"\nElapsed: {elapsed:.1f} min")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n Exception: {type(e).__name__}: {e}", file=sys.stderr)
        raise
