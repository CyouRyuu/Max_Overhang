#!/usr/bin/env bash
# Exact-arithmetic checks; no optimization solver needed (Python + mpmath only).
#  * builds the 13 optimal stacks N = 4..16 in exact rational arithmetic, verifies equilibrium with zero residual,
#    and writes their fractions to data/rational_stacks.json (Appendix A),
#  * re-checks them with an independent checker and compares with S*(N),
#  * checks the base-singleton chain (S2+ table) and the non-spinal chain margins (Section 6.1).
set -euo pipefail
source "$(dirname "$0")/common.sh"
"$PY" "$SRC/verify_rational_stacks.py" --dump-json "$OUT/rational_stacks.json"
cmp -s "$OUT/rational_stacks.json" "$REPO/data/rational_stacks.json" \
    && echo "rational_stacks.json identical to data/rational_stacks.json" \
    || { echo "WARNING: $OUT/rational_stacks.json differs from data/rational_stacks.json"; exit 1; }
"$PY" "$SRC/check_exact_stacks.py"
