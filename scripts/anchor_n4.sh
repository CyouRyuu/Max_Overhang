#!/usr/bin/env bash
# Base anchor of the S2+ chain: monolithic M0 at N = 4, relative gap 1e-9 (zbar(4) = 1.1678933).
set -euo pipefail
source "$(dirname "$0")/common.sh"
"$PY" "$SRC/certify_run.py" 4 --model M0 --gap 1e-9 --threads 1 --time-limit 600 --outdir "$OUT/anchor_n4"
