#!/usr/bin/env bash
# Stability of the certificates under tighter numerical settings: re-solve the three most expensive profiles
# of the tight certificates at N = 14 and N = 15 with FeasibilityTol = OptimalityTol = IntFeasTol = 1e-9 and
# NumericFocus = 3 (same model, same cutoff).  Every profile must close again (CUTOFF or INFEASIBLE).
#   scripts/tolerance_check.sh [WORKERS]   (N = 16 non-spinal sub-branches: scripts/tolerance_check_ns16.sh)
# Optional, long: N = 16, most expensive profile (12,475 s at default settings):
#   EXTRA16=1 scripts/tolerance_check.sh
set -euo pipefail
source "$(dirname "$0")/common.sh"
W="${1:-3}"
TOL=(--param FeasibilityTol=1e-9 --param OptimalityTol=1e-9 --param NumericFocus=3 --ptag TOL)
declare -A HARD=( [14]="1-2-2-2-3-3-1,1-2-2-2-2-2-2-1,1-2-2-3-3-2-1"
                  [15]="1-2-2-2-2-2-3-1,1-2-2-2-3-2-2-1,1-2-2-2-2-2-2-2" )
[[ -n "${EXTRA16:-}" ]] && HARD[16]="1-2-2-2-2-2-3-2"
for N in $(echo "${!HARD[@]}" | tr ' ' '\n' | sort -n); do
  "$PY" "$SRC/profile_run.py" "$N" "${CERT_OPTS[@]}" "${TOL[@]}" --cutoff "${CUT[$N]}" --only "${HARD[$N]}" \
      --workers "$W" --time-limit 172800 --outdir "$OUT/tolerance_N$N"
done
