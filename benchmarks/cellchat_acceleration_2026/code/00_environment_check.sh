#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="${ROOT:-$BENCH_ROOT}"
DATA_ROOT="${DATA_ROOT:-$ROOT/data}"
OUT="$ROOT/results/environment"
mkdir -p "$OUT"

{
  echo "timestamp: $(date -Is)"
  echo "host: $(hostname)"
  echo "kernel: $(uname -a)"
  echo "cpu_threads: $(nproc 2>/dev/null || echo NA)"
  echo
  echo "memory:"
  free -h || true
  echo
  echo "disk:"
  df -h "$ROOT" "$DATA_ROOT" || true
  echo
  echo "executables:"
  command -v R || true
  command -v Rscript || true
  command -v python3 || true
  command -v git || true
} | tee "$OUT/system.txt"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available in PATH." | tee "$OUT/R_missing.txt"
  exit 2
fi

Rscript - <<'RSCRIPT' | tee "$OUT/R_packages.txt"
cat("R.version.string:", R.version.string, "\n")
pkgs <- c(
  "Seurat", "CellChat", "CellChatAccelRcpp", "Matrix", "dplyr",
  "data.table", "ggplot2", "peakRAM"
)
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(p)) else "MISSING"
  cat(sprintf("%-22s %s\n", p, ver))
}
cat("\nSession info:\n")
print(sessionInfo())
RSCRIPT

echo "Environment check complete: $OUT"

