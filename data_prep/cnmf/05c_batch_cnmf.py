"""
Step 5c: Batch-schedule per-patient cNMF on 9950X (16C32T).

Strategy:
    Serial per-patient. Each patient uses 8 factorize workers. No CCD
    splitting, no small/large distinction. 05b internally sets 1 BLAS
    thread per worker so 8 workers ~= 8 cores.

Resume:
    Skip any patient whose cnmf_input/<key>/DONE flag exists.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("${PROJECT_ROOT}")
CNMF_INPUT = ROOT / "data" / "cnmf_input"
RUN_ORDER = CNMF_INPUT / "run_order.txt"
SCRIPT_05B = Path(__file__).resolve().parent / "05b_run_cnmf.py"

N_WORKERS = 8


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def read_run_order():
    if not RUN_ORDER.exists():
        raise FileNotFoundError(f"{RUN_ORDER} missing; run 05a first.")
    patients = []
    with open(RUN_ORDER) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            key, n_cells = parts[0], int(parts[1])
            patients.append((key, n_cells))
    return patients


def already_done(patient_dir: Path) -> bool:
    return (patient_dir / "DONE").exists()


def spawn(patient_dir: Path, n_workers: int):
    cmd = [
        sys.executable, str(SCRIPT_05B),
        str(patient_dir), str(n_workers),
    ]
    log(f"  spawn: {SCRIPT_05B.name} {patient_dir.name} {n_workers}")
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    log_path = patient_dir / "run.log"
    lf = open(log_path, "a")
    p = subprocess.Popen(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)
    p._log_fh = lf
    return p


def wait_one(p):
    ret = p.wait()
    try:
        p._log_fh.close()
    except Exception:
        pass
    return ret


def main():
    all_patients = read_run_order()
    log(f"Loaded {len(all_patients)} patients from run_order.txt")

    todo = []
    skipped = 0
    for key, n in all_patients:
        pdir = CNMF_INPUT / key
        if already_done(pdir):
            skipped += 1
            continue
        todo.append((key, n, pdir))
    log(f"Skipping {skipped} already-DONE patients; {len(todo)} remaining")

    total = len(todo)
    t0 = time.time()
    completed = 0
    for key, n, pdir in todo:
        log(f"[{completed+1}/{total}] {key} ({n} cells)")
        p = spawn(pdir, N_WORKERS)
        r = wait_one(p)
        if r != 0:
            log(f"  FAIL {key}: exit {r}")
        completed += 1
        elapsed = time.time() - t0
        eta = elapsed / completed * (total - completed) if completed else 0
        log(
            f"  progress {completed}/{total}  elapsed {elapsed/60:.1f}m  "
            f"ETA {eta/60:.1f}m"
        )

    log(f"Batch finished. Total wall: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
