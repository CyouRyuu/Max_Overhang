#!/usr/bin/env bash
# Exact spinal optima (Appendix C); no optimization solver needed.
#  1. standard-library interval verifier (src/exact_spines/): S(N) = S*(N) for N = 4, 6, 8..19, the stacks for
#     N = 5, 7, all margins and the exception faces; its CSV files must equal results/exact_spines/
#  2. optional SymPy cross-check of that verifier: balance identities in Q(x) and the root polynomials
#  3. second verifier: face bounds at N = 5, 7 (src/spine_faces.py) and the algebraic stacks over Q(x)
#     (src/check_algebraic_spines.py, N = 4..NMAX, default 19)
# Everything takes well under a minute.
# Needs Python with sympy and mpmath.
set -euo pipefail
source "$(dirname "$0")/common.sh"
E="$OUT/exact_spines"
mkdir -p "$E"
"$PY" "$SRC/exact_spines/overhang_spine_check.py" > "$E/original_stdout.txt"
cmp -s "$E/original_stdout.txt" "$REPO/results/exact_spines/original_stdout.txt" \
    && echo "overhang_spine_check.py: output identical to the archived one" \
    || { echo "WARNING: overhang_spine_check.py output differs"; exit 1; }
"$PY" "$SRC/exact_spines/audit_spines.py" --out "$E"
"$PY" "$SRC/exact_spines/symbolic_check.py" --out "$E"
for f in "$REPO"/results/exact_spines/*.csv; do
    cmp -s "$f" "$E/$(basename "$f")" || { echo "WARNING: $(basename "$f") differs from the archived file"; exit 1; }
done
echo "all CSV files identical to results/exact_spines/"
"$PY" "$SRC/spine_faces.py" --leaves "$E/face_leaves" | tee "$E/spine_faces.log"
for f in "$REPO"/results/exact_spines/face_leaves/*.csv; do
    cmp -s "$f" "$E/face_leaves/$(basename "$f")" || { echo "WARNING: $(basename "$f") differs from the archived file"; exit 1; }
done
"$PY" "$SRC/check_face_leaves.py" "$E/face_leaves"
"$PY" "$SRC/check_algebraic_spines.py" --nmax "${NMAX:-19}" | tee "$E/check_algebraic_spines.log"
