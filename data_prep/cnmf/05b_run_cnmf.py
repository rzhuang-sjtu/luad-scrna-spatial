"""
Step 5b: Run cNMF full pipeline for a single patient.

Usage:
    python 05b_run_cnmf.py <patient_dir> <n_workers>

cNMF 1.7.x parallelism note:
    cNMF.factorize() does NOT auto-parallelize. Jobs are split via
    (worker_i, total_workers). To use all cores we spawn `n_workers`
    subprocesses, each running factorize() with its own worker_i, and
    each constrained to 1 BLAS thread via env vars below.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import time
import subprocess
from pathlib import Path

import numpy as np

K_RANGE = np.arange(8, 21)
N_ITER = 200
SEED = 42
NUM_HIGHVAR_GENES = 3000
BETA_LOSS = "frobenius"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_prepare(output_dir, name, counts_path):
    from cnmf import cNMF
    obj = cNMF(output_dir=output_dir, name=name)
    obj.prepare(
        counts_fn=counts_path,
        components=K_RANGE,
        n_iter=N_ITER,
        seed=SEED,
        num_highvar_genes=NUM_HIGHVAR_GENES,
        beta_loss=BETA_LOSS,
    )
    return obj


def run_factorize_worker(output_dir, name, worker_i, total_workers):
    """Entrypoint for a single factorize worker (invoked as subprocess)."""
    from cnmf import cNMF
    obj = cNMF(output_dir=output_dir, name=name)
    obj.factorize(worker_i=worker_i, total_workers=total_workers,
                  skip_completed_runs=True)


def spawn_factorize_workers(output_dir, name, n_workers):
    procs = []
    for wi in range(n_workers):
        cmd = [
            sys.executable, __file__, "__worker__",
            output_dir, name, str(wi), str(n_workers),
        ]
        env = os.environ.copy()
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OMP_NUM_THREADS"] = "1"
        env["NUMEXPR_NUM_THREADS"] = "1"
        p = subprocess.Popen(cmd, env=env)
        procs.append(p)
    failed = []
    for i, p in enumerate(procs):
        ret = p.wait()
        if ret != 0:
            failed.append((i, ret))
    if failed:
        raise RuntimeError(f"Factorize workers failed: {failed}")


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "__worker__":
        # Subprocess entry point for a factorize worker.
        _, _, output_dir, name, worker_i, total_workers = sys.argv
        run_factorize_worker(output_dir, name, int(worker_i), int(total_workers))
        return

    if len(sys.argv) != 3:
        print("Usage: python 05b_run_cnmf.py <patient_dir> <n_workers>")
        sys.exit(1)

    patient_dir = Path(sys.argv[1]).resolve()
    n_workers = int(sys.argv[2])

    patient_key = patient_dir.name
    counts_path = patient_dir / "counts.h5ad"
    output_dir = patient_dir / "cnmf_out"

    if not counts_path.exists():
        raise FileNotFoundError(f"{counts_path} not found")

    log(f"=== {patient_key} | workers={n_workers} ===")
    t_start = time.time()

    done_flag = patient_dir / "DONE"
    if done_flag.exists():
        log(f"Already DONE, skipping: {patient_key}")
        return

    log("Step 1/4: prepare")
    run_prepare(str(output_dir), patient_key, str(counts_path))

    log(f"Step 2/4: factorize (spawn {n_workers} workers)")
    spawn_factorize_workers(str(output_dir), patient_key, n_workers)

    log("Step 3/4: combine")
    from cnmf import cNMF
    obj = cNMF(output_dir=str(output_dir), name=patient_key)
    obj.combine()

    log("Step 4/4: k_selection_plot")
    try:
        obj.k_selection_plot(close_fig=True)
    except TypeError:
        obj.k_selection_plot()

    done_flag.touch()
    dt = time.time() - t_start
    log(f"Done: {patient_key} in {dt/60:.1f} min")


if __name__ == "__main__":
    main()
