#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="${ROOT:-$BENCH_ROOT}"
DATA_ROOT="${DATA_ROOT:-$ROOT/data}"
JOBS="${JOBS:-64}"
MAX_DATASETS="${MAX_DATASETS:-12}"
REPEATS="${REPEATS:-3}"
SCALES="${SCALES:-1000,5000,10000,25000,50000,all}"
ENGINES="${ENGINES:-both,accelerated}"
ABLATIONS="${ABLATIONS:-no_accel_kernel,no_accel_pathway,no_accel_aggregate}"
ACCEL_ALGORITHMS="${ACCEL_ALGORITHMS:-dense}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export BLIS_NUM_THREADS="${BLIS_NUM_THREADS:-1}"

cd "$ROOT"

python3 code/01_build_manifest.py \
  --data-root "$DATA_ROOT" \
  --out results/data_manifest.tsv \
  --candidates results/dataset_candidates.tsv

python3 code/02_make_experiment_grid.py \
  --manifest results/data_manifest.tsv \
  --out results/experiment_grid.csv \
  --max-datasets "$MAX_DATASETS" \
  --repeats "$REPEATS" \
  --scales "$SCALES" \
  --engines "$ENGINES" \
  --ablations "$ABLATIONS" \
  --accel-algorithms "$ACCEL_ALGORITHMS"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is not available. Install or activate R, then rerun this launcher." >&2
  exit 2
fi

N_TASKS="$(python3 - <<'PY'
import csv
with open("results/experiment_grid.csv", newline="") as fh:
    print(sum(1 for _ in csv.DictReader(fh)))
PY
)"

echo "Launching $N_TASKS tasks with JOBS=$JOBS"
seq 1 "$N_TASKS" | xargs -P "$JOBS" -I{} \
  Rscript code/run_grid.R \
    --root "$ROOT" \
    --grid "$ROOT/results/experiment_grid.csv" \
    --task-id {}

Rscript code/summarize_results.R

