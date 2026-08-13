#!/usr/bin/env python3
"""
Step 4d: Batch-schedule per-patient CopyKAT; serial runs with resume support.

- Order from run_order.txt (smaller patients first)
- Skip completed patients (copykat_prediction.csv present)
- Log failures to failed_patients.txt
"""
import json
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

INPUT_DIR = Path("${PROJECT_ROOT}/data/copykat_input")
RUN_ORDER = INPUT_DIR / "run_order.txt"
FAILED_LOG = INPUT_DIR / "failed_patients.txt"
R_SCRIPT = Path("${PROJECT_ROOT}/data_prep/copykat/04c_run_copykat.R")

N_CORES = 8


def log(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def read_run_order():
    items = []
    for line in RUN_ORDER.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        items.append({"dp_key": parts[0], "n_cells": int(parts[1]), "ref_mode": parts[2]})
    return items


def already_done(pdir: Path) -> bool:
    pred = pdir / "copykat_prediction.csv"
    return pred.exists() and pred.stat().st_size > 0


def fmt_eta(seconds: float) -> str:
    if seconds < 0 or seconds != seconds:  # NaN
        return "?"
    return str(timedelta(seconds=int(seconds)))


def main():
    if not RUN_ORDER.exists():
        log(f"{RUN_ORDER} not found; run 04b_prepare_copykat.py first")
        sys.exit(1)

    items = read_run_order()
    total = len(items)
    log(f"Total patients: {total}")

    # Pre-scan done / todo
    done, todo = [], []
    for it in items:
        pdir = INPUT_DIR / it["dp_key"]
        (done if already_done(pdir) else todo).append(it)
    log(f"Done: {len(done)}  Pending: {len(todo)}")

    if not todo:
        log("All patients already completed; exit.")
        return

    # Batch run
    t_start = time.time()
    n_ok = 0
    n_fail = 0
    cum_cells = 0
    cum_cells_total = sum(it["n_cells"] for it in todo)

    for i, it in enumerate(todo, 1):
        dp = it["dp_key"]
        pdir = INPUT_DIR / dp
        log(f"--- [{i}/{len(todo)}] {dp}  n_cells={it['n_cells']}  mode={it['ref_mode']} ---")

        t_patient = time.time()
        try:
            proc = subprocess.run(
                ["Rscript", str(R_SCRIPT), str(pdir), str(N_CORES)],
                check=False,
                capture_output=True,
                text=True,
            )
        except Exception as e:
            log(f"Failed to start Rscript: {e}")
            n_fail += 1
            with FAILED_LOG.open("a") as f:
                f.write(f"{dp}\tstartup_error\t{e}\n")
            continue

        dt = time.time() - t_patient

        # Print key R stdout lines
        if proc.stdout:
            for line in proc.stdout.splitlines():
                if any(tok in line for tok in ("START", "DONE", "FAILED", "mode:", "matrix:",
                                                "prediction rows", "label distribution",
                                                "", "", "aneuploid", "diploid", "not.defined")):
                    print(f"    {line}", flush=True)

        if proc.returncode == 0 and already_done(pdir):
            n_ok += 1
            cum_cells += it["n_cells"]
            elapsed = time.time() - t_start
            rate = cum_cells / max(elapsed, 1e-6)   # cells / sec
            remain_cells = sum(x["n_cells"] for x in todo[i:])
            eta = remain_cells / max(rate, 1e-6)
            log(f"Done ({dt/60:.1f} min). Cumulative ok={n_ok} fail={n_fail}.  ETA={fmt_eta(eta)}")
        else:
            n_fail += 1
            stderr_tail = "\n".join(proc.stderr.splitlines()[-20:]) if proc.stderr else ""
            log(f"Failed (rc={proc.returncode}, {dt/60:.1f} min)")
            if stderr_tail:
                log("stderr tail:")
                for line in stderr_tail.splitlines():
                    print(f"    {line}", flush=True)
            with FAILED_LOG.open("a") as f:
                f.write(f"{dp}\trc={proc.returncode}\t{stderr_tail[:300]}\n")

    total_elapsed = time.time() - t_start
    log("=" * 60)
    log(f"Batch finished. ok={n_ok}  fail={n_fail}  skip={len(done)}")
    log(f"Total elapsed: {fmt_eta(total_elapsed)}")
    if n_fail > 0:
        log(f"Failure details: {FAILED_LOG}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted by user; re-run this script to resume")
        sys.exit(130)
