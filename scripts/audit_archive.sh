#!/usr/bin/env bash
# Re-audit the archived records in results/ without any solver (Python and mpmath only).
#  * tight certificates N = 5..16, each with the base-singleton (S2+) condition checked against the audited
#    certificate of N - 1 (the M0 anchor for N = 4); analytic closures (energy test, no candidate layer) are
#    reported separately from solver closures;
#  * non-spinal exclusions N = 5..16 (t0 split at N = 16), each with the condition for stacks with n_2 = 1
#    checked against the audited certificate and exclusion of N - 1;
#  * the baselines, the matched comparison, the reverse controls and the exact N = 9 non-spinal stack;
#  * the negative tests of both auditors (tampered copies of archived records must be rejected);
#  * the exact re-check of every analytic screening decision (energy test, dropped candidate layers, bound conflicts
#    of non-spinal branches) in rational arithmetic (results/screening/);
#  * a summary table of the whole chain (results of this run in $OUT/audit/).
set -euo pipefail
source "$(dirname "$0")/common.sh"
R="$REPO/results"; A="$OUT/audit"; mkdir -p "$A"
ANCHOR="$R/anchor_n4/M0cert_N4_s42_result.json"
FAIL=0
for N in $(seq 5 16); do
  CSV=$(ls "$R"/certificates/*PDcert_N${N}_s42_profiles.csv)
  PREV=$([ "$N" = 5 ] && echo "$ANCHOR" || echo "$A/cert_N$((N - 1)).json")
  "$PY" "$SRC/pd_certificate.py" "$N" "$CSV" --cutoff "${CUT[$N]}" --prev-cert "$PREV" --out "$A/cert_N$N" \
      | sed -n 3,9p || FAIL=1
done
for N in $(seq 5 16); do
  PREVC=$([ "$N" = 5 ] && echo "$ANCHOR" || echo "$A/cert_N$((N - 1)).json")
  PREVN=(); [ "$N" -ge 6 ] && PREVN=(--prev-ns "$A/ns_N$((N - 1)).json")
  SPLIT=(); [ "$N" = 16 ] && SPLIT=(--split "$R/nonspinal/ns_N16_rel0.01_only_t0.csv")
  "$PY" "$SRC/ns_audit.py" "$N" "$R/nonspinal/ns_N${N}_rel0.01.csv" "${SPLIT[@]}" --prev-cert "$PREVC" \
      "${PREVN[@]}" --out "$A/ns_N$N.json" | tail -2 || FAIL=1
done
B="$R/baselines/nosp_N13"
"$PY" "$SRC/pd_certificate.py" 13 "$(ls "$B"/*_profiles.csv)" --race "$B/race" --split "$B/split_1-2-2-2-3-2-1" \
    --cutoff 1.9855055567 --allow-missing-run-record --prev-cert "$A/cert_N12.json" --out "$A/nosp_cert_N13" \
    | sed -n 3,9p || FAIL=1
# Matched comparison without family SP (Section 6.2): N = 11 is a complete certificate; N = 12 and 13 leave
# 2 and 4 profiles open at the time limit, so their audit reports "certified: False" by design.
for N in 11 12 13; do
  "$PY" "$SRC/pd_certificate.py" "$N" "$(ls "$R"/ablation_sp/N$N/*_profiles.csv)" --cutoff "${CUT[$N]}" \
      --prev-cert "$A/cert_N$((N - 1)).json" --out "$A/ablation_N$N" | sed -n 3,6p || true
done
# Reverse controls (Section 6.2): N = 8 non-spinal branches closed at threshold 1.60, which also excludes the N = 9
# stacks with n_2 = 1 above 1.60 + 1/18; N = 9 closed at rho = 3e-2; the N = 9 run at 1.60 is open by design and
# its stack, made exactly balanced in rational arithmetic, is non-spinal in both senses.
"$PY" "$SRC/ns_audit.py" 8 "$R/reverse_controls/n8/ns_N8_g1.6000000.csv" --gamma 1.60 \
    --prev-cert "$A/cert_N7.json" --out "$A/rc_ns_N8.json" | tail -2 || FAIL=1
"$PY" "$SRC/ns_audit.py" 9 "$R/reverse_controls/n9/ns_N9_rel0.03.csv" --rel 3e-2 \
    --prev-ns "$A/rc_ns_N8.json" --out "$A/rc_ns_N9.json" | tail -2 || FAIL=1
"$PY" "$SRC/exact_witness.py" --check "$R/reverse_controls/n9/ns_N9_exact_stack.json" || FAIL=1
"$PY" "$SRC/test_auditors.py" | tail -1 || FAIL=1
"$PY" "$SRC/check_screening.py" --out "$A/screening" | tail -1 || FAIL=1
for f in "$REPO"/results/screening/*; do cmp -s "$f" "$A/screening/$(basename "$f")" || { echo "screening output differs: $f"; FAIL=1; }; done
"$PY" "$SRC/chain_summary.py" "$A" || FAIL=1
exit $FAIL
