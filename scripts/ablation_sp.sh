#!/usr/bin/env bash
# Matched comparison of the spinal cuts (Section 6.2): the certificate runs of N = 11, 12, 13 repeated with
# the same options, cutoff, energy filter, tolerances and one thread per subproblem, but without family SP.
# Time limit per profile: 3,600 s at N = 11, 12 and 14,400 s at N = 13 (the paper used 7 workers).
#   scripts/ablation_sp.sh [WORKERS]
set -euo pipefail
source "$(dirname "$0")/common.sh"
W="${1:-7}"
NOSP_OPTS=("${CERT_OPTS[@]}")
NOSP_OPTS=("${NOSP_OPTS[@]/--spine}")                  # drop family SP, keep everything else
NOSP_OPTS=("${NOSP_OPTS[@]/SP1/NOSPABL}")
declare -A TL=( [11]=3600 [12]=3600 [13]=14400 )
for N in 11 12 13; do
  "$PY" "$SRC/profile_run.py" "$N" ${NOSP_OPTS[@]} --cutoff "${CUT[$N]}" --workers "$W" \
      --time-limit "${TL[$N]}" --outdir "$OUT/ablation_sp/N$N"
  "$PY" "$SRC/pd_certificate.py" "$N" "$(ls "$OUT"/ablation_sp/N$N/*_profiles.csv)" --cutoff "${CUT[$N]}" \
      --out "$OUT/ablation_sp/audit_N$N" | sed -n 3,6p || true
done
