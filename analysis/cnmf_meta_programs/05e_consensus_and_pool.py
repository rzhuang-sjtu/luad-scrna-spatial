"""
Step 5e: per-patient cNMF consensus at patient-specific K, then pool GEPs.

- Reads K per patient from results/step5d_k_selection_summary.csv
- For each patient:
    cnmf_obj.consensus(k=selected_k, density_threshold=DENSITY_THRESHOLD,
                       show_clustering=False, close_clustergram_fig=True)
    (u, spectra_scores, spectra_tpm, top_genes) = cnmf_obj.load_results(...)
    Saves per-patient:
        cnmf_out/<key>/<key>_gep_spectra_scores_k{K}.csv
        cnmf_out/<key>/<key>_gep_top100_k{K}.csv
- Aggregates all GEP spectra into:
    data/cnmf_output/gep_pool.csv            # rows = gene, cols = GEP_id
    data/cnmf_output/gep_pool_zscore.csv     # same layout, spectra_scores (z-scored)
    data/cnmf_output/gep_metadata.csv        # one row per GEP
    data/cnmf_output/gep_top100_genes.csv    # long table of top-100 genes per GEP
- Writes a CONSENSUS_DONE flag per patient for resume.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MPLBACKEND"] = "Agg"

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

ROOT = Path("${PROJECT_ROOT}")
CNMF_INPUT = ROOT / "data" / "cnmf_input"
CNMF_OUTPUT = ROOT / "data" / "cnmf_output"
CNMF_OUTPUT.mkdir(parents=True, exist_ok=True)
RESULTS = ROOT / "results"

K_SUMMARY = RESULTS / "step5d_k_selection_summary.csv"
DENSITY_THRESHOLD = 0.1
N_TOP_GENES = 100


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_consensus_for_patient(pdir: Path, k: int):
    """Run cNMF consensus + load results. Returns (spectra_scores, spectra_tpm, usage, top_genes)."""
    from cnmf import cNMF
    key = pdir.name
    output_dir = pdir / "cnmf_out"
    obj = cNMF(output_dir=str(output_dir), name=key)
    obj.consensus(
        k=k,
        density_threshold=DENSITY_THRESHOLD,
        show_clustering=False,
        close_clustergram_fig=True,
        refit_usage=True,
    )
    usage, spectra_scores, spectra_tpm, top_genes = obj.load_results(
        K=k, density_threshold=DENSITY_THRESHOLD, n_top_genes=N_TOP_GENES
    )
    return usage, spectra_scores, spectra_tpm, top_genes


def save_per_patient(pdir: Path, k: int, spectra_scores, spectra_tpm, top_genes):
    key = pdir.name
    out = pdir / "cnmf_out" / key
    # spectra_scores / spectra_tpm: gene rows x program cols
    spectra_scores.to_csv(out / f"{key}_spectra_scores_k{k}.csv")
    spectra_tpm.to_csv(out / f"{key}_spectra_tpm_k{k}.csv")
    top_genes.to_csv(out / f"{key}_top{N_TOP_GENES}_k{k}.csv")


def main():
    if not K_SUMMARY.exists():
        raise FileNotFoundError(K_SUMMARY)
    summary = pd.read_csv(K_SUMMARY)
    log(f"Loaded {len(summary)} patients from {K_SUMMARY.name}")

    gep_score_frames = []   # list of (K x gene) DataFrames, one per patient
    gep_tpm_frames = []
    gep_meta_rows = []
    top_gene_rows = []

    n_ok = 0
    n_fail = 0
    for i, row in summary.iterrows():
        key = row["patient_key"]
        k = int(row["selected_k"])
        pdir = CNMF_INPUT / key
        if not (pdir / "DONE").exists():
            log(f"[{i+1}/{len(summary)}] SKIP {key}: no DONE flag")
            continue

        flag = pdir / f"CONSENSUS_DONE_k{k}"
        try:
            if flag.exists():
                # Reuse saved results
                from cnmf import cNMF
                obj = cNMF(output_dir=str(pdir / "cnmf_out"), name=key)
                usage, spectra_scores, spectra_tpm, top_genes = obj.load_results(
                    K=k, density_threshold=DENSITY_THRESHOLD, n_top_genes=N_TOP_GENES
                )
                log(f"[{i+1}/{len(summary)}] {key} K={k} (reused)")
            else:
                t0 = time.time()
                usage, spectra_scores, spectra_tpm, top_genes = run_consensus_for_patient(pdir, k)
                save_per_patient(pdir, k, spectra_scores, spectra_tpm, top_genes)
                flag.touch()
                log(f"[{i+1}/{len(summary)}] {key} K={k} OK in {time.time()-t0:.1f}s")
            n_ok += 1
        except Exception as e:
            log(f"[{i+1}/{len(summary)}] FAIL {key}: {e}")
            n_fail += 1
            continue

        # spectra_scores, spectra_tpm: rows = gene, cols = program idx (1..K)
        # Transpose to (programs x genes) and label with global GEP ids.
        prog_ids = [f"{key}__K{k}__P{int(p)}" for p in spectra_scores.columns]
        ss = spectra_scores.T.copy()
        ss.index = prog_ids
        ss.index.name = "gep_id"
        gep_score_frames.append(ss)

        st = spectra_tpm.T.copy()
        st.index = prog_ids
        st.index.name = "gep_id"
        gep_tpm_frames.append(st)

        # top_genes: rows = rank 0..N-1, cols = program idx
        for rank, row_genes in top_genes.iterrows():
            for col, gene in row_genes.items():
                top_gene_rows.append({
                    "gep_id": f"{key}__K{k}__P{int(col)}",
                    "rank": int(rank) + 1,
                    "gene": gene,
                })

        for p in spectra_scores.columns:
            gep_meta_rows.append({
                "gep_id": f"{key}__K{k}__P{int(p)}",
                "patient_key": key,
                "k_total": k,
                "program_idx": int(p),
                "n_cells": int(row["n_cells"]),
                "stability": float(row["stability_at_k"]),
            })

    log(f"consensus finished: {n_ok} OK, {n_fail} FAIL")

    # Pool: union of genes across patients, wide matrix (gene rows x GEP cols).
    log("Pooling GEPs into gene x GEP matrices ...")
    all_genes = sorted(set().union(*[df.columns for df in gep_score_frames]))
    log(f"Union gene count: {len(all_genes)}")

    def pool(frames):
        parts = []
        for df in frames:
            reindexed = df.reindex(columns=all_genes)
            parts.append(reindexed)
        pool_wide = pd.concat(parts, axis=0)  # GEPs x genes
        return pool_wide.T  # genes x GEPs

    pool_scores = pool(gep_score_frames)
    pool_tpm = pool(gep_tpm_frames)

    log(f"GEP pool (spectra_scores): {pool_scores.shape[0]} genes x {pool_scores.shape[1]} GEPs")

    pool_scores.to_csv(CNMF_OUTPUT / "gep_pool_zscore.csv")
    pool_tpm.to_csv(CNMF_OUTPUT / "gep_pool_tpm.csv")
    # Convenience: default "gep_pool.csv" == z-scores (used for Spearman clustering)
    pool_scores.to_csv(CNMF_OUTPUT / "gep_pool.csv")

    meta = pd.DataFrame(gep_meta_rows)
    # Derive dataset from patient_key robustly
    meta["dataset"] = meta["patient_key"].str.split("__").str[0]
    meta = meta[["gep_id", "patient_key", "dataset", "k_total", "program_idx",
                 "n_cells", "stability"]]
    meta.to_csv(CNMF_OUTPUT / "gep_metadata.csv", index=False)

    top_df = pd.DataFrame(top_gene_rows)
    top_df.to_csv(CNMF_OUTPUT / "gep_top100_genes.csv", index=False)

    log("Wrote:")
    for f in ["gep_pool.csv", "gep_pool_zscore.csv", "gep_pool_tpm.csv",
              "gep_metadata.csv", "gep_top100_genes.csv"]:
        p = CNMF_OUTPUT / f
        log(f"  {p}  ({p.stat().st_size/1024/1024:.1f} MB)")

    log(f"Total GEPs pooled: {len(meta)}  from {meta['patient_key'].nunique()} patients")
    log(f"Dataset distribution:")
    log(meta.groupby("dataset")["gep_id"].count().to_string())


if __name__ == "__main__":
    main()
