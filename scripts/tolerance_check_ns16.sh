#!/usr/bin/env bash
# Stability of the N = 16 non-spinal exclusion under tighter numerical settings: re-solve the 14 t0 sub-branches
# of Appendix B (sub-branch table), which close the three main branches left open at 7,200 s, with
# FeasibilityTol = OptimalityTol = IntFeasTol = 1e-9 and NumericFocus = 3 (same model, same threshold).
# Every sub-branch must close again (CUTOFF, INFEASIBLE or TIE).  Reference at default settings: 4,779 s wall,
# 16,244 thread-s.
#   scripts/tolerance_check_ns16.sh [WORKERS]
set -euo pipefail
source "$(dirname "$0")/common.sh"
W="${1:-11}"
R="$REPO/results"
"$PY" "$SRC/ns_check.py" 16 --rel 1e-2 \
    --only "1-2-2-2-2-2-3-2:7,1-2-2-2-2-3-3-1:7,1-2-2-2-3-3-2-1:6" \
    --split-t0 --workers "$W" --threads 1 --time-limit 172800 \
    --param FeasibilityTol=1e-9 --param OptimalityTol=1e-9 --param NumericFocus=3 \
    --outdir "$OUT/tolerance_ns16"
"$PY" "$SRC/ns_audit.py" 16 "$R/nonspinal/ns_N16_rel0.01.csv" \
    --split "$OUT/tolerance_ns16/ns_N16_rel0.01_only_t0.csv"
