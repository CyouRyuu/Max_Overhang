#!/usr/bin/env bash
# Re-audit the archived records in results/ without any solver (Python and mpmath only):
# every tight certificate N = 5..16 and every non-spinal exclusion N = 5..16 (with the t0 split at N = 16).
set -euo pipefail
source "$(dirname "$0")/common.sh"
R="$REPO/results"; mkdir -p "$OUT/audit"
for N in $(seq 5 16); do
  CSV=$(ls "$R"/certificates/*PDcert_N${N}_s42_profiles.csv)
  "$PY" "$SRC/pd_certificate.py" "$N" "$CSV" --cutoff "${CUT[$N]}" --out "$OUT/audit/cert_N$N" | sed -n 3,6p
done
B="$R/baselines/nosp_N13"
"$PY" "$SRC/pd_certificate.py" 13 "$(ls "$B"/*_profiles.csv)" --race "$B/race" --split "$B/split_1-2-2-2-3-2-1" \
    --cutoff 1.9855055567 --allow-missing-run-record --out "$OUT/audit/nosp_cert_N13" | sed -n 3,6p
for N in $(seq 5 15); do "$PY" "$SRC/ns_audit.py" "$N" "$R/nonspinal/ns_N${N}_rel0.01.csv" | tail -1; done
"$PY" "$SRC/ns_audit.py" 16 "$R/nonspinal/ns_N16_rel0.01.csv" --split "$R/nonspinal/ns_N16_rel0.01_only_t0.csv" | tail -1
# Matched comparison without family SP (Section 6.2): N = 11 is a complete certificate; N = 12 and 13 leave
# 2 and 4 profiles open at the time limit, so their audit reports "certified: False" by design.
for N in 11 12 13; do
  "$PY" "$SRC/pd_certificate.py" "$N" "$(ls "$R"/ablation_sp/N$N/*_profiles.csv)" --cutoff "${CUT[$N]}" \
      --out "$OUT/audit/ablation_N$N" | sed -n 3,6p || true
done
# Reverse controls (Section 6.2): N = 8 non-spinal branches closed at threshold 1.60 (the n_2 = 1 step for
# N = 9) and N = 9 closed at rho = 3e-2; the N = 9 run at 1.60 is open by design.
"$PY" "$SRC/ns_audit.py" 8 "$R/reverse_controls/n8/ns_N8_g1.6000000.csv" --gamma 1.60 | tail -1
"$PY" "$SRC/ns_audit.py" 9 "$R/reverse_controls/n9/ns_N9_rel0.03.csv" --rel 3e-2 | tail -1
