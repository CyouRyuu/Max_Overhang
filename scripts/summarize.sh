#!/usr/bin/env bash
# Print a summary of all certificate runs found under runs/ (or the directory given).
set -euo pipefail
source "$(dirname "$0")/common.sh"
for f in $(ls "${1:-$OUT}"/cert_N*/certificate_N*.md 2>/dev/null); do head -20 "$f"; echo; done
cat "${1:-$OUT}"/nonspinal*/ns_summary.txt 2>/dev/null || true
