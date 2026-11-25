# Shared settings for the reproduction scripts (sourced, not executed).
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO/src"
PY="${PYTHON:-python3}"
OUT="${OUT:-$REPO/runs}"            # all new outputs go here; the archived results are in results/
mkdir -p "$OUT"

# Tight cutoffs gamma_N = z_Q(N) + 1e-5, rounded upward at the 7th decimal (paper, Section 6.1).
declare -A CUT=( [5]=1.3045649 [6]=1.4367103 [7]=1.5300553 [8]=1.6315169 [9]=1.7152829 [10]=1.7871359
                 [11]=1.8587896 [12]=1.9251014 [13]=1.9845190 [14]=2.0382270 [15]=2.0929064 [16]=2.1438528 )

# Options of the spinal-cut layer-profile certificate (identical for every N).
CERT_OPTS=(--model M1u --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --energy --spine
           --param Heuristics=0.0 --param IntFeasTol=1e-9 --ptag SP1)
