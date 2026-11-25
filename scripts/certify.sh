#!/usr/bin/env bash
# Tight global-optimality certificate for one N in 5..16 (main certificate table, left part):
# spinal-cut layer-profile decomposition over P_2(N) at cutoff gamma_N, then certificate assembly/audit.
#   scripts/certify.sh N [WORKERS]          e.g.  scripts/certify.sh 13 7
# Reference effort (single-thread subproblems): N <= 12 seconds to minutes, N = 15 about 2.4 thread-hours,
# N = 16 about 20 thread-hours.  Re-running with the same arguments plus --resume continues a stopped run:
#   EXTRA="--resume" scripts/certify.sh 16 7
set -euo pipefail
source "$(dirname "$0")/common.sh"
N="$1"; W="${2:-4}"
[[ -n "${CUT[$N]:-}" ]] || { echo "N must be in 5..16"; exit 1; }
D="$OUT/cert_N$N"
"$PY" "$SRC/profile_run.py" "$N" "${CERT_OPTS[@]}" --cutoff "${CUT[$N]}" --workers "$W" \
    --time-limit 259200 --nodefile-start 1.5 --soft-mem-limit 6 --outdir "$D" ${EXTRA:-}
CSV=$(ls "$D"/*PDcert_N${N}_s42_profiles.csv)
"$PY" "$SRC/pd_certificate.py" "$N" "$CSV" --cutoff "${CUT[$N]}" --out "$D/certificate_N$N"
