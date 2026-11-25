#!/usr/bin/env bash
# N = 16: the three branches left open by scripts/nonspinal.sh 16, partitioned by the lowest impure layer t0
# (paper: t0-split proposition; Appendix B, N = 16 sub-branch table).  14 sub-branches: 3 TIE (closed without
# computation), 11 CUTOFF.  Reference: 4,779 s wall, 16,244 thread-s.
#   scripts/nonspinal_n16_split.sh [WORKERS]
set -euo pipefail
source "$(dirname "$0")/common.sh"
W="${1:-11}"
"$PY" "$SRC/ns_check.py" 16 --rel 1e-2 \
    --only "1-2-2-2-2-2-3-2:7,1-2-2-2-2-3-3-1:7,1-2-2-2-3-3-2-1:6" \
    --split-t0 --workers "$W" --threads 1 --time-limit 86400 --outdir "$OUT/nonspinal_n16_split"
