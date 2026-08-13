"""Step 17: export CopyKAT CNA matrices for Fig S2-CNA panels.

Sources:
  ~/luad/data/copykat_input/<dataset>__<patient>/<...>_copykat_CNA_results.txt
  ~/luad/data/processed/luad_copykat.h5ad         (malignant flag + dominant_MP)
  ~/luad/results/step7_mp_cell_scores.csv         (per-cell dominant_MP)

Outputs → ${WORK_ROOT}/luad_figures/fig_s2_cna/:
  cna_by_sample.csv.gz       bins × patients (mean CNA over malignant cells)
  cna_by_mp.csv.gz           bins × {MP1,MP2,MP3,MP4} (mean CNA over malignant cells)
  cna_sample_metadata.csv    per-patient stats (n_cells, n_malignant, majority MP)
  cna_patient_dist.csv       optional: pairwise Euclidean distance between patient CNA profiles
"""
from __future__ import annotations
import os, glob, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

CNA_ROOT = Path.home()/"luad/data/copykat_input"
H5AD = Path.home()/"luad/data/processed/luad_copykat.h5ad"
MP_CSV = Path.home()/"luad/results/step7_mp_cell_scores.csv"
OUT = Path("${WORK_ROOT}/luad_figures/fig_s2_cna")
OUT.mkdir(parents=True, exist_ok=True)


def log(msg): print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    t0 = time.time()

    log("loading h5ad obs (malignant flag + dominant MP)")
    a = sc.read_h5ad(H5AD, backed="r")
    obs = a.obs[["dataset", "patient_id", "malignant"]].copy()
    obs["malignant_bool"] = obs["malignant"].astype(str) == "Malignant"
    log(f"  h5ad cells: {len(obs)}; malignant: {int(obs['malignant_bool'].sum())}")

    # Per-cell dominant_MP from step7 (only for malignant cells)
    mp = pd.read_csv(MP_CSV, index_col=0)["dominant_MP"].astype(str).to_dict()
    log(f"  cells with dominant_MP: {len(mp)}")

    cna_files = sorted(glob.glob(str(CNA_ROOT/"*"/"*_copykat_CNA_results.txt")))
    log(f"  CNA result files: {len(cna_files)}")

    first = cna_files[0]
    bin_meta = pd.read_csv(first, sep="\t", usecols=[0, 1, 2])
    bin_meta.columns = ["chr", "chrompos", "abspos"]
    log(f"  shared bin grid: {len(bin_meta)} bins; first 3:\n{bin_meta.head(3).to_string()}")

    log("processing each patient's CNA matrix")
    sample_means = {}        # patient_key → np.array(n_bins,)
    sample_meta_rows = []
    mp_sums = {f"MP{i}": np.zeros(len(bin_meta), dtype=np.float64) for i in range(1, 5)}
    mp_counts = {f"MP{i}": 0 for i in range(1, 5)}

    for i, f in enumerate(cna_files, 1):
        patient_key = Path(f).parent.name  # e.g. "GSE164789__LHL"
        try:
            dataset_id, patient_id = patient_key.split("__", 1)
        except ValueError:
            log(f"  [{i}/{len(cna_files)}] skip {patient_key} (bad name)")
            continue

        # Read header to learn columns; columns 4+ are cells (R-formatted with dots)
        try:
            cna = pd.read_csv(f, sep="\t")
        except Exception as e:
            log(f"  [{i}/{len(cna_files)}] read fail {patient_key}: {e}")
            continue

        # Verify bin alignment
        if len(cna) != len(bin_meta):
            log(f"  WARN {patient_key}: nbins={len(cna)} != {len(bin_meta)}; skipping")
            continue

        cell_cols = list(cna.columns[3:])  # skip chrom/chrompos/abspos
        # Normalize: R dots → original dashes (CopyKAT col name has '.' but our barcode has '-')
        normalized = [c.replace(".", "-") for c in cell_cols]
        col_to_norm = dict(zip(cell_cols, normalized))

        # Find malignant cells in obs that match this patient
        sub_obs = obs[(obs["dataset"] == dataset_id)
                       & (obs["patient_id"] == patient_id)
                       & (obs["malignant_bool"])]
        if len(sub_obs) == 0:
            log(f"  [{i}/{len(cna_files)}] {patient_key}: no malignant cells in h5ad — skip")
            continue

        malignant_barcodes = set(sub_obs.index)
        # Match: select cna columns whose normalized name is in malignant_barcodes
        matched_cols = [c for c in cell_cols if col_to_norm[c] in malignant_barcodes]

        if len(matched_cols) < 5:
            log(f"  [{i}/{len(cna_files)}] {patient_key}: only {len(matched_cols)} matched malignant cells — skip")
            continue

        # Mean CNA over malignant cells
        mean_cna = cna[matched_cols].astype(np.float32).mean(axis=1).values
        sample_means[patient_key] = mean_cna

        # Determine majority MP for this patient
        mp_for_patient = [mp.get(col_to_norm[c]) for c in matched_cols]
        mp_for_patient = [m for m in mp_for_patient if m is not None]
        if mp_for_patient:
            from collections import Counter
            ctr = Counter(mp_for_patient)
            top_mp, top_n = ctr.most_common(1)[0]
            mp_breakdown = ";".join(f"{k}:{v}" for k, v in ctr.most_common())
        else:
            top_mp, top_n, mp_breakdown = "Unknown", 0, ""

        # Update MP-level mean accumulator (cell-weighted within MP)
        # For per-MP mean across all patients, we need per-cell CNA × MP labels
        # Cheaper: per-patient compute MP-specific mean → weight by cell count
        for mp_lab in ("MP1","MP2","MP3","MP4"):
            mp_cols = [c for c in matched_cols if mp.get(col_to_norm[c]) == mp_lab]
            if len(mp_cols) >= 5:
                # Sum (not mean) of cells in this MP; divide later
                mp_sums[mp_lab] += cna[mp_cols].astype(np.float32).sum(axis=1).values
                mp_counts[mp_lab] += len(mp_cols)

        sample_meta_rows.append({
            "sample_id": patient_key,
            "dataset": dataset_id,
            "patient_id": patient_id,
            "n_cells_total": len(cell_cols),
            "n_malignant_in_h5ad": len(sub_obs),
            "n_malignant_in_CNA": len(matched_cols),
            "dominant_MP_majority": top_mp,
            "MP_breakdown": mp_breakdown,
        })
        if i % 10 == 0:
            log(f"  [{i}/{len(cna_files)}] processed {patient_key} ({len(matched_cols)} malig)")

    log(f"\nfinished {len(sample_means)} patients with malignant CNA")

    log("building bin × patient matrix")
    sample_df = pd.DataFrame(sample_means, index=bin_meta.index)
    sample_df = pd.concat([bin_meta.rename(columns={"chrompos":"start","abspos":"abspos"})
                           .assign(end=bin_meta["chrompos"]),
                           sample_df], axis=1)
    # Re-order columns: chr, start, end, abspos, then patient cols
    cols_first = ["chr", "start", "end", "abspos"]
    cols_pat = [c for c in sample_df.columns if c not in cols_first]
    sample_df = sample_df[cols_first + cols_pat]
    sample_df.to_csv(OUT/"cna_by_sample.csv.gz", index=False, compression="gzip")
    log(f"  cna_by_sample.csv.gz {sample_df.shape}")

    log("building bin × MP matrix (cell-weighted across patients)")
    mp_means = {}
    for mp_lab in ("MP1","MP2","MP3","MP4"):
        if mp_counts[mp_lab] > 0:
            mp_means[mp_lab] = mp_sums[mp_lab] / mp_counts[mp_lab]
        else:
            mp_means[mp_lab] = np.full(len(bin_meta), np.nan)
    log(f"  cells per MP used: {mp_counts}")
    mp_df = pd.DataFrame(mp_means, index=bin_meta.index)
    mp_df = pd.concat([bin_meta.rename(columns={"chrompos":"start","abspos":"abspos"})
                        .assign(end=bin_meta["chrompos"]),
                        mp_df], axis=1)
    mp_df = mp_df[cols_first + list(mp_means.keys())]
    mp_df.to_csv(OUT/"cna_by_mp.csv.gz", index=False, compression="gzip")
    log(f"  cna_by_mp.csv.gz {mp_df.shape}")

    meta_df = pd.DataFrame(sample_meta_rows)
    meta_df.to_csv(OUT/"cna_sample_metadata.csv", index=False)
    log(f"  cna_sample_metadata.csv {meta_df.shape}")

    log("computing pairwise patient CNA distance (euclidean)")
    pat_matrix = np.array(list(sample_means.values()))  # (n_patients, n_bins)
    pat_keys = list(sample_means.keys())
    from scipy.spatial.distance import pdist, squareform
    d = squareform(pdist(pat_matrix, metric="euclidean"))
    dist_df = pd.DataFrame(d, index=pat_keys, columns=pat_keys)
    dist_df.to_csv(OUT/"cna_patient_dist.csv")
    log(f"  cna_patient_dist.csv {dist_df.shape}")

    log(f"\nDONE in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
