#!/usr/bin/env bash
# Positive reverse controls (Section 6.1): with the cutoff lowered, the same models must find feasible stacks.
#  (a) N = 13, profile (1,2,2,2,3,2,1) with family SP, cutoff 1.9844089 -> a spinal stack of overhang ~1.9845089
#  (b) N = 7, non-spinal branches with gamma = 1.38 -> not closed; a non-spinal stack of overhang ~1.3842371
set -euo pipefail
source "$(dirname "$0")/common.sh"
"$PY" "$SRC/profile_run.py" 13 "${CERT_OPTS[@]}" --cutoff 1.9844089 --only 1-2-2-2-3-2-1 \
    --time-limit 3600 --outdir "$OUT/reverse_n13"
"$PY" "$SRC/ns_check.py" 7 --gamma 1.38 --workers 4 --outdir "$OUT/reverse_n7" || true
