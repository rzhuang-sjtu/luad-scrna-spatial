#!/usr/bin/env bash
# G3 — launch parallel inferCNV with three-layer memory guardrails.
#
# Why guardrails: largest patient has up to 9,185 epithelial cells; inferCNV HMM peaks at 30–40 GB.
# Concurrent large patients can trip the system OOM killer, which kills processes at random,
#
# Three-layer guardrails:
#   1) Tier by cell count: small patients in parallel; large patients serial and last to avoid bunching large jobs
#   2) Each R process self-caps virtual memory with ulimit -v and exits on overrun,
#      rather than letting the system OOM killer kill other jobs
#   3) Watchdog every 60 s: if free memory is below threshold, pause dispatching new patients (leave running ones alone);
#      resume automatically when free memory recovers
#
# Resume: skip patients that already have summary/<patient>_cells.csv; re-run this script after interrupt.
#
# Usage: bash G3_launch_infercnv.sh
# Overridable: NPROC parallel slots / MEM_PER_JOB_GB per-process cap / MIN_FREE_GB pause threshold /
#         BIG_CELLS large-patient cutoff
set -uo pipefail
# inferCNV was run on a rented GPU node; ROOT there was a scratch path.
ROOT=${INFERCNV_ROOT:-${WORK_ROOT:-.}/infercnv}
IN=$ROOT/input; SUM=$ROOT/summary; LOG=$ROOT/logs
mkdir -p "$SUM" "$LOG"

# System R is 4.1.2 / Bioconductor 3.14; full dependency stack would need building from source.
# Use conda binaries instead: infercnv 1.26.0 + R 4.5.3 with rjags.
# Point RSCRIPT at the Rscript of an environment holding infercnv 1.26.0 and rjags;
# the distribution R on that node was too old for the dependency stack.
RSCRIPT=${RSCRIPT:-Rscript}

TOTAL_GB=$(free -g | awk '/^Mem:/{print $2}')
# Measured (32 cores / 503G machine, 13 concurrent jobs):
#   each job averages only 1.2 cores — num_threads=2 applies only in a few parallel sections; the main path is single-threaded.
#   Memory far below estimate: largest patient (9,685 cells) RSS only 3.5G; 13 jobs total 21.6G.
# Therefore set concurrency to cores/1.2; memory caps need only a few-fold headroom.
NPROC=${NPROC:-20}            # concurrent slots for small patients
NPROC_BIG=${NPROC_BIG:-8}     # concurrent slots for large patients (runs alongside the small pool; 8 covers all large ones)
MEM_PER_JOB_GB=${MEM_PER_JOB_GB:-20}
MEM_BIG_GB=${MEM_BIG_GB:-60}
MIN_FREE_GB=${MIN_FREE_GB:-40}
BIG_CELLS=${BIG_CELLS:-4000}

# Scheduling: large patients are long jobs and must start first (longest first); otherwise they serialise at the end and dominate wall time.
# Two pools run together: large NPROC_BIG + small NPROC, each counted at num_threads=2.
echo "Total memory ${TOTAL_GB}G"
echo "Large patients (>${BIG_CELLS} cells) concurrency ${NPROC_BIG}, per-process cap ${MEM_BIG_GB}G"
echo "Small patients concurrency ${NPROC}, per-process cap ${MEM_PER_JOB_GB}G"
echo "Pause new dispatches when free memory < ${MIN_FREE_GB}G (leave running jobs alone)"

LIST=$LOG/patient_order.txt
: > "$LIST"
for d in "$IN"/*/; do
  p=$(basename "$d")
  [ -f "$SUM/${p}_cells.csv" ] && continue
  n=$(wc -l < "$d/cells.tsv" 2>/dev/null || echo 0)
  printf '%s\t%s\n' "$n" "$p" >> "$LIST"
done
sort -n -o "$LIST" "$LIST"
NTODO=$(wc -l < "$LIST"); NBIG=$(awk -v b="$BIG_CELLS" '$1>b' "$LIST" | wc -l)
echo "${NTODO} patients pending, of which ${NBIG} large"
[ "$NTODO" -eq 0 ] && { echo "All done"; exit 0; }

PAUSE=$LOG/PAUSE; rm -f "$PAUSE"
(
  while :; do
    f=$(free -g | awk '/^Mem:/{print $7}')
    if [ "${f:-999}" -lt "$MIN_FREE_GB" ]; then
      [ -f "$PAUSE" ] || { touch "$PAUSE"; echo "$(date +%H:%M:%S) free ${f}G pause dispatch" >> "$LOG/watchdog.log"; }
    else
      [ -f "$PAUSE" ] && { rm -f "$PAUSE"; echo "$(date +%H:%M:%S) free ${f}G resume dispatch" >> "$LOG/watchdog.log"; }
    fi
    sleep 60
  done
) & WD=$!
trap 'kill $WD 2>/dev/null' EXIT

run_one () {
  local p=$1 lim=$2
  ( ulimit -v $(( lim * 1024 * 1024 ))
    nice -n 5 "$RSCRIPT" "$ROOT/G2_run_infercnv.R" --single "$p" ) > "$LOG/${p}.log" 2>&1
  [ -f "$SUM/${p}_cells.csv" ] || echo "$(date +%H:%M:%S) $p failed" >> "$LOG/failed.log"
}

# Large-patient pool: start from the largest
(
  awk -v b="$BIG_CELLS" '$1>b{print $2}' "$LIST" | tac | while read -r p; do
    while [ -f "$PAUSE" ] || [ "$(jobs -rp | wc -l)" -ge "$NPROC_BIG" ]; do sleep 15; done
    echo "$(date +%H:%M:%S) [large] dispatch $p"
    run_one "$p" "$MEM_BIG_GB" &
  done
  wait
) & BIGPOOL=$!

# Small-patient pool: start from the smallest to raise throughput first
(
  awk -v b="$BIG_CELLS" '$1<=b{print $2}' "$LIST" | while read -r p; do
    while [ -f "$PAUSE" ] || [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 10; done
    run_one "$p" "$MEM_PER_JOB_GB" &
  done
  wait
) & SMALLPOOL=$!

wait $BIGPOOL $SMALLPOOL

echo "All finished $(date +%H:%M:%S)"
echo "Done $(find "$SUM" -maxdepth 1 -type f -name "*_cells.csv" | wc -l) / 89 个"
[ -f "$LOG/failed.log" ] && { echo "Failure log:"; cat "$LOG/failed.log"; }
