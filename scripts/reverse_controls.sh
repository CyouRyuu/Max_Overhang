#!/usr/bin/env bash
# Positive reverse controls (Section 6): with the threshold lowered, the same models must find feasible stacks.
# Every feasible stack found is saved with positions, contact forces and contact points (a witness) and
# re-checked without a solver by src/check_witness.py.
#  (a) N = 13, profile (1,2,2,2,3,2,1) with family SP, cutoff 1.9844089 -> a spinal stack of overhang ~1.9845089.
#      Run with the default Gurobi heuristics, IntFeasTol = 1e-6 and relative gap 1e-4 (one thread): with the
#      certificate settings (Heuristics = 0, IntFeasTol = 1e-9) the solver finds no incumbent within an hour,
#      while its dual bound stalls at S^_5(13) = 1.9845089, which is also consistent but yields no witness.
#  (b) N = 7, non-spinal branches with gamma = 1.38 -> not closed; a stack of overhang ~1.3842371 that is
#      non-spinal under the closed resting relation (it carries force through a point contact)
#  (c) N = 9, non-spinal branches with gamma = 1.60 -> not closed; a stack of overhang ~1.63189 that is
#      non-spinal under the closed resting relation (with Gurobi 13.0.1 it also uses a point contact)
#  (d) N = 8, non-spinal branches with gamma = 1.60 -> all closed (1.60 > zbar(7) + 1/16 = 1.5925553); this
#      excludes N = 9 non-spinal stacks with n_2 = 1 above 1.60 + 1/18 < 0.97 S*(9) (Corollary, condition (56))
#  (e) N = 9, non-spinal branches at rho = 3e-2 (gamma = 0.97 S*(9) ~ 1.66381) -> all closed
set -euo pipefail
source "$(dirname "$0")/common.sh"
D="$OUT/reverse_controls"; mkdir -p "$D"
"$PY" "$SRC/profile_run.py" 13 --model M1u --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --energy --spine \
    --param IntFeasTol=1e-6 --gap 1e-4 --ptag NEG --cutoff 1.9844089 --only 1-2-2-2-3-2-1 \
    --time-limit 3600 --outdir "$D/n13_spinal"
"$PY" "$SRC/check_witness.py" "$(ls "$D"/n13_spinal/*_result.json | head -1)"
"$PY" "$SRC/ns_check.py" 7 --gamma 1.38 --workers 4 --outdir "$D/n7" || true
"$PY" "$SRC/check_witness.py" "$D/n7/ns_N7_g1.3800000_witnesses.json" 1
"$PY" "$SRC/ns_check.py" 9 --gamma 1.60 --workers 4 --outdir "$D/n9" || true
"$PY" "$SRC/check_witness.py" "$D/n9/ns_N9_g1.6000000_witnesses.json" 1
"$PY" "$SRC/ns_check.py" 8 --gamma 1.60 --workers 4 --outdir "$D/n8"
"$PY" "$SRC/ns_check.py" 9 --rel 3e-2 --workers 4 --outdir "$D/n9"
