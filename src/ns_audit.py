"""Audit of the non-spinal exclusion records written by ns_check.py.

Checks, for one N:
  * the expected branches (n, l*) are exactly those ns_check.py enumerates (profiles of P_2(N) that pass the
    energy filter, l* = 2..H): no branch missing, no unexpected branch, repeated records reported;
  * every branch is closed (INFEASIBLE, CUTOFF, EMPTY, or TIE), or, if it is open in the main run, it is
    covered by a t0 split (--split) in which every expected sub-branch t0 in J_{l*} is present and closed;
  * the threshold: gamma_ns(N) = ceil_7((1 - rel) S*(N)) is computed so that rounding errors cannot change the
    seventh decimal (the computed value must lie at least 1e-5 grid units from a grid point), and every closed
    record (CUTOFF, INFEASIBLE) carries the threshold its run used, which must not exceed gamma_ns(N); the solver
    threshold of ns_check.py (double precision) is checked in the same way.  Needs no solver.

  python3 src/ns_audit.py 16 results/nonspinal/ns_N16_rel0.01.csv --split results/nonspinal/ns_N16_rel0.01_only_t0.csv
"""
import argparse
import csv
import math
import sys
from collections import Counter

import profiles as pf

CLOSED = {"INFEASIBLE", "CUTOFF", "EMPTY", "TIE"}


def read(path, with_t0):
    recs, cnt = {}, Counter()
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = (r["profile"], int(r["lstar"])) + ((int(r["t0"]),) if with_t0 and r.get("t0") not in (None, "") else ())
            cnt[key] += 1
            recs[key] = r
    return recs, {k: c for k, c in cnt.items() if c > 1}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("N", type=int)
    ap.add_argument("csv", help="main ns_check.py record (ns_N<N>_rel0.01.csv)")
    ap.add_argument("--split", default=None, help="t0-split record of the branches left open (optional)")
    ap.add_argument("--rel", type=float, default=1e-2)
    ap.add_argument("--gamma", type=float, default=None,
                    help="absolute threshold of the run (ns_check.py --gamma); overrides --rel")
    a = ap.parse_args()

    expected = set()
    J = {}
    for p in pf.enumerate_profiles(a.N, 2):
        ok, _, _ = pf.energy_profile_check(a.N, p)
        if not ok:
            continue
        ptxt = "-".join(map(str, p))
        for l in range(2, len(p) + 1):
            expected.add((ptxt, l))
            J[(ptxt, l)] = [t for t in range(3, l + 1) if p[t - 2] >= 2]

    recs, dup = read(a.csv, with_t0=False)
    missing = sorted(expected - set(recs))
    unexpected = sorted(set(recs) - expected)
    open_main = sorted(k for k in expected & set(recs) if recs[k]["status"] not in CLOSED)

    split_missing, split_open, split_unexpected, split_dup = [], [], [], {}
    if a.split:
        srecs, split_dup = read(a.split, with_t0=True)
        want = {k + (t,) for k in open_main for t in J[k]} | {k for k in open_main if not J[k]}
        split_missing = sorted(want - set(srecs))
        split_unexpected = sorted(set(srecs) - want)
        split_open = sorted(k for k in want & set(srecs) if srecs[k]["status"] not in CLOSED)
        still_open = []
    else:
        want, still_open = set(), open_main

    if a.gamma is not None:                  # absolute threshold, used as given (reverse controls, Section 6.2)
        solver_gamma = gamma_ns = a.gamma
        robust = True
    else:
        solver_gamma = pf.spinal_opt(a.N) * (1.0 - a.rel)
        scaled = solver_gamma * 1e7
        frac = scaled - math.floor(scaled)
        robust = 1e-5 < frac < 1.0 - 1e-5    # float error (~1e-8 grid units) cannot move the value across a grid point
        gamma_ns = math.ceil(scaled) / 1e7
    # every CUTOFF record must carry the finite threshold of its run; INFEASIBLE records may carry it or -inf
    used, bad_bound = [], []
    for r in list(recs.values()) + (list(srecs.values()) if a.split else []):
        if r["status"] not in ("CUTOFF", "INFEASIBLE"):
            continue
        raw = (r.get("bound") or "").strip()
        if r["status"] == "INFEASIBLE" and raw in ("", "-inf"):
            continue
        try:
            b = float(raw)
        except ValueError:
            b = float("nan")
        if math.isfinite(b):
            used.append(b)
        else:
            bad_bound.append((r["profile"], r["lstar"], r.get("t0", ""), r["status"], raw))
    above = [b for b in used if b > gamma_ns]
    gamma_ok = robust and solver_gamma <= gamma_ns and not above and not bad_bound

    ok = not (missing or unexpected or still_open or split_missing or split_open or split_unexpected) and gamma_ok
    print(f"N={a.N}: expected branches {len(expected)}, records {len(recs)}, missing {len(missing)}, "
          f"unexpected {len(unexpected)}, repeated {len(dup)}, open in main run {len(open_main)}")
    if a.split:
        print(f"  t0 split: expected sub-branches {len(want)}, missing {len(split_missing)}, "
              f"unexpected {len(split_unexpected)}, repeated {len(split_dup)}, open {len(split_open)}; "
              f"status counts {dict(Counter(srecs[k]['status'] for k in want & set(srecs)))}")
    print(f"  solver threshold {solver_gamma:.15f} <= gamma_ns(N) = {gamma_ns:.7f} (rounding robust: {robust}); "
          f"closed records with a recorded threshold: {len(used)}, max {max(used) if used else float('nan'):.15f}, "
          f"above gamma_ns: {len(above)}")
    for name, v in [("missing", missing), ("unexpected", unexpected), ("repeated", dup),
                    ("open (no split given)", still_open), ("split missing", split_missing),
                    ("split unexpected", split_unexpected), ("split open", split_open),
                    ("closed records without a finite threshold", bad_bound)]:
        if v:
            print(f"  {name}: {v}")
    print(f"  NON-SPINAL EXCLUSION {'CERTIFIED' if ok else 'NOT CERTIFIED'} above gamma_ns({a.N}) = {gamma_ns:.7f}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
