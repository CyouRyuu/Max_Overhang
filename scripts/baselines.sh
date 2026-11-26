#!/usr/bin/env bash
# Methodological-progression table (Section 6.2), N = 12 and N = 13, Gurobi 13.0.1.  These runs are LONG.
#   scripts/baselines.sh m0_12        monolithic M0, N = 12, 4 threads, 72 h
#   scripts/baselines.sh m0_13        monolithic M0, N = 13, 4 threads, 48 h (open at gap 2.9e-3)
#   scripts/baselines.sh pd_12        profile decomposition without SP, N = 12, 4 workers (about 11 min)
#   scripts/baselines.sh pd_13        profile decomposition without SP, N = 13, 4 workers; the profile
#                                     (1,2,2,2,3,2,1) dominates the effort
#   scripts/baselines.sh pd_13_race   the same profile on 3 threads in a separate process (about 40 min)
#   scripts/baselines.sh pd_13_split  independent cross-check: rightmost-layer split of (1,2,2,2,3,2,1)
#   scripts/baselines.sh pd_13_cert   assemble the N = 13 certificate without SP (after pd_13 + pd_13_race)
# The cutoffs are those of the original runs: LB (1 + 5e-4) with the feasible values available at the time.
set -euo pipefail
source "$(dirname "$0")/common.sh"
NOSP=(--model M1u --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --param Heuristics=0.0 --ptag NOSP)
CUT12=1.9260577476        # 1.9250952 * 1.0005
CUT13=1.9855055567        # 1.9845133 * 1.0005
case "${1:-}" in
  m0_12) "$PY" "$SRC/certify_m0.py" 12 --threads 4 --time-limit 259200 --outdir "$OUT/m0" ;;
  m0_13) "$PY" "$SRC/certify_m0.py" 13 --threads 4 --time-limit 172800 --outdir "$OUT/m0" ;;
  pd_12) "$PY" "$SRC/profile_run.py" 12 "${NOSP[@]}" --cutoff $CUT12 --workers 4 --time-limit 14400 \
             --outdir "$OUT/nosp_N12" ;;
  pd_13) "$PY" "$SRC/profile_run.py" 13 "${NOSP[@]}" --cutoff $CUT13 --workers 4 --time-limit 259200 \
             --keep-logs --nodefile-start 1.5 --soft-mem-limit 6 --outdir "$OUT/nosp_N13" ${EXTRA:-} ;;
  pd_13_race)
         CSV=$(ls "$OUT"/nosp_N13/*_profiles.csv)
         "$PY" "$SRC/pd_race.py" 13 "$CSV" --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer \
             --param Heuristics=0.0 --cutoff $CUT13 --threads 3 --seed 1 --outdir "$OUT/nosp_N13_race" ;;
  pd_13_split)
         "$PY" "$SRC/pd_split.py" 13 1-2-2-2-3-2-1 --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer \
             --param Heuristics=0.0 --cutoff $CUT13 --threads-budget 3 --outdir "$OUT/nosp_N13_split" ;;
  pd_13_cert)
         CSV=$(ls "$OUT"/nosp_N13/*_profiles.csv)
         # the main run is stopped by hand once the race has closed the hard profile, so it has no result record
         "$PY" "$SRC/pd_certificate.py" 13 "$CSV" --race "$OUT/nosp_N13_race" --cutoff $CUT13 \
             --split "$OUT/nosp_N13_split" --split-profile 1-2-2-2-3-2-1 --allow-missing-run-record \
             --out "$OUT/nosp_N13/certificate_N13" ;;
  *) sed -n '2,12p' "$0"; exit 1 ;;
esac
