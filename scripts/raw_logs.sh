#!/usr/bin/env bash
# Raw Gurobi logs and model files (MPS, gzip) for key instances, with the settings of the archived runs:
#   * the tight certificates and the non-spinal exclusions at N = 5 and N = 7 (every solved profile and branch),
#   * the 11 t0 sub-branches of the N = 16 non-spinal exclusion that need the solver (sub-branch table of Appendix B).
# Every record also stores the closure kind (analytic or solver), the raw Gurobi status code, ObjBound and the
# cutoff.  The audits at the end must pass again.  Reference: N = 5, 7 take seconds; N = 16 about 4,800 s wall
# with 11 workers (16,244 thread-s in the archived run).
#   scripts/raw_logs.sh [WORKERS]
set -euo pipefail
source "$(dirname "$0")/common.sh"
W="${1:-11}"
R="$REPO/results"
D="$OUT/raw_logs"
mkdir -p "$D"
for N in 5 7; do
    "$PY" "$SRC/profile_run.py" "$N" "${CERT_OPTS[@]}" --cutoff "${CUT[$N]}" --workers 4 \
        --time-limit 259200 --keep-logs --write-models --outdir "$D/cert_N$N"
    PREV="$( [[ $N == 5 ]] && ls "$R"/anchor_n4/*_result.json || echo "$R/certificates/cert_N$((N - 1)).json")"
    "$PY" "$SRC/pd_certificate.py" "$N" "$(ls "$D"/cert_N$N/*PDcert_N${N}_s42_profiles.csv)" --cutoff "${CUT[$N]}" \
        --prev-cert "$PREV" --out "$D/cert_N$N/certificate_N$N"
    "$PY" "$SRC/ns_check.py" "$N" --rel 1e-2 --workers 4 --logdir "$D/ns_N$N/logs" --outdir "$D/ns_N$N"
    "$PY" "$SRC/ns_audit.py" "$N" "$D/ns_N$N/ns_N${N}_rel0.01.csv" --prev-cert "$PREV" --out "$D/ns_N$N/ns_audit_N$N.json"
done
"$PY" "$SRC/ns_check.py" 16 --rel 1e-2 \
    --only "1-2-2-2-2-2-3-2:7,1-2-2-2-2-3-3-1:7,1-2-2-2-3-3-2-1:6" \
    --split-t0 --workers "$W" --threads 1 --time-limit 86400 \
    --logdir "$D/ns_N16_t0/logs" --outdir "$D/ns_N16_t0"
"$PY" "$SRC/ns_audit.py" 16 "$R/nonspinal/ns_N16_rel0.01.csv" \
    --split "$D/ns_N16_t0/ns_N16_rel0.01_only_t0.csv" --out "$D/ns_N16_t0/ns_audit_N16.json"
echo "RAW LOGS WRITTEN TO $D"
