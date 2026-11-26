#!/usr/bin/env bash
# Non-spinal exclusion at rho = 1e-2 (gamma_ns = 0.99 S*(N)) for one N in 5..16
# (main certificate table, right part).  Per-branch time limit 7200 s.
#   scripts/nonspinal.sh N [WORKERS]
# At N = 16 three branches reach the time limit; they are closed by scripts/nonspinal_n16_split.sh.
set -euo pipefail
source "$(dirname "$0")/common.sh"
N="$1"; W="${2:-7}"
"$PY" "$SRC/ns_check.py" "$N" --rel 1e-2 --workers "$W" --outdir "$OUT/nonspinal"
CSV="$OUT/nonspinal/ns_N${N}_rel0.01.csv"
[[ "$N" == 16 ]] || "$PY" "$SRC/ns_audit.py" "$N" "$CSV"
