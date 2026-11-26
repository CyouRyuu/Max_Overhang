"""Audit of the non-spinal exclusion records written by ns_check.py.

Checks, for one N:
  * the expected branches (n, l*) are exactly those ns_check.py enumerates (profiles of P_2(N) that pass the
    energy filter, l* = 2..H): no branch missing, no unexpected branch, repeated records reported;
  * every branch is closed (INFEASIBLE, CUTOFF, EMPTY, or TIE), or, if it is open in the main run, it is
    covered by a t0 split (--split) in which every expected sub-branch t0 in J_{l*} is present and closed;
  * the analytic statuses are checked mechanically: a main record is EMPTY exactly when J_{l*} is empty, a
    split record is TIE exactly when t0 = l*, and a recorded closure kind (newer runs) agrees with the status;
    repeated records with conflicting closed statuses fail the audit;
  * the threshold: gamma_ns(N) = ceil_7((1 - rel) S*(N)) is computed so that rounding errors cannot change the
    seventh decimal (the computed value must lie at least 1e-5 grid units from a grid point), and every closed
    record (CUTOFF, INFEASIBLE) carries the threshold its run used, which must not exceed gamma_ns(N); the solver
    threshold of ns_check.py (double precision) is checked in the same way;
  * with --prev-cert / --prev-ns, the condition of the corollary that excludes the stacks with n_2 = 1:
    zbar(N-1) + 1/(2N) < gamma_ns(N), or gamma_ns(N-1) + 1/(2N) < gamma_ns(N) with the exclusion at N - 1
    certified (its JSON from --out).  Without these options only the branch coverage is audited, and the
    report says so.  Needs no solver.

  python3 src/ns_audit.py 16 results/nonspinal/ns_N16_rel0.01.csv --split results/nonspinal/ns_N16_rel0.01_only_t0.csv
"""
import argparse
import csv
import json
import math
import sys
from collections import Counter

import profiles as pf

CLOSED = {"INFEASIBLE", "CUTOFF", "EMPTY", "TIE"}


KIND = {"EMPTY": "analytic: J empty", "TIE": "analytic: t0 = l*"}


def read(path, with_t0):
    recs, cnt, hist = {}, Counter(), {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = (r["profile"], int(r["lstar"])) + ((int(r["t0"]),) if with_t0 and r.get("t0") not in (None, "") else ())
            cnt[key] += 1
            hist.setdefault(key, []).append(r["status"])
            recs[key] = r
    conflicts = {k: h for k, h in hist.items() if len(h) > 1 and
                 (len({x for x in h if x in CLOSED}) > 1 or (h[-1] not in CLOSED and any(x in CLOSED for x in h)))}
    return recs, {k: c for k, c in cnt.items() if c > 1}, conflicts


def status_errors(recs, J, split):
    """Mechanical check of the analytic statuses (see the module docstring)."""
    bad = []
    for k, r in recs.items():
        st, kind = r["status"], (r.get("closure") or "").strip()
        jk = J.get(k[:2])
        if jk is None:
            continue                                   # unexpected record, reported elsewhere
        if split:
            t0 = k[2] if len(k) > 2 else None
            ok = (st == "TIE") == (t0 is not None and t0 == k[1]) and st != "EMPTY"
        else:
            ok = (st == "EMPTY") == (not jk) and st != "TIE"
        if kind:
            ok = ok and (KIND.get(st) == kind if st in KIND else not kind.startswith("analytic"))
        if not ok:
            bad.append(k + (st, kind))
    return sorted(bad)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("N", type=int)
    ap.add_argument("csv", help="main ns_check.py record (ns_N<N>_rel0.01.csv)")
    ap.add_argument("--split", default=None, help="t0-split record of the branches left open (optional)")
    ap.add_argument("--rel", type=float, default=1e-2)
    ap.add_argument("--gamma", type=float, default=None,
                    help="absolute threshold of the run (ns_check.py --gamma); overrides --rel")
    ap.add_argument("--prev-cert", default=None,
                    help="certificate JSON of N - 1 (pd_certificate.py output, or the M0 anchor result for N = 5)")
    ap.add_argument("--prev-ns", default=None, help="this script's JSON output (--out) for N - 1")
    ap.add_argument("--out", default=None, help="write a JSON summary")
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

    recs, dup, conflicts = read(a.csv, with_t0=False)
    status_bad = status_errors(recs, J, split=False)
    missing = sorted(expected - set(recs))
    unexpected = sorted(set(recs) - expected)
    open_main = sorted(k for k in expected & set(recs) if recs[k]["status"] not in CLOSED)

    split_missing, split_open, split_unexpected, split_dup = [], [], [], {}
    if a.split:
        srecs, split_dup, split_conflicts = read(a.split, with_t0=True)
        status_bad += status_errors(srecs, J, split=True)
        conflicts = {**conflicts, **split_conflicts}
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

    ok = not (missing or unexpected or still_open or split_missing or split_open or split_unexpected
              or status_bad or conflicts) and gamma_ok

    # Stacks with n_2 = 1 are not enumerated; the corollary excludes them above gamma_ns(N) if
    #   zbar(N-1) + 1/(2N) < gamma_ns(N)          (general condition, zbar from the certificate of N - 1), or
    #   gamma_ns(N-1) + 1/(2N) < gamma_ns(N)      (sharper condition, needs the exclusion at N - 1).
    chain = {"checked": False}
    if a.prev_cert or a.prev_ns:
        chain["checked"] = True
        if a.prev_cert:
            pc = json.load(open(a.prev_cert))
            if "certified" in pc:
                ub, cok = pc["cutoff_upper_bound"], bool(pc["certified"]) and pc.get("N") == a.N - 1
            else:
                ub, cok = math.ceil(pc["bound"] * 1e7 - 1e-6) / 1e7, pc.get("status") == "OPTIMAL"
            chain["general_margin"] = gamma_ns - (ub + 0.5 / a.N)
            chain["general_ok"] = cok and chain["general_margin"] > 0
        if a.prev_ns:
            pn = json.load(open(a.prev_ns))
            chain["sharper_margin"] = gamma_ns - (pn["gamma_ns"] + 0.5 / a.N)
            chain["sharper_ok"] = bool(pn["certified"]) and pn.get("N") == a.N - 1 and chain["sharper_margin"] > 0
        chain["ok"] = bool(chain.get("general_ok") or chain.get("sharper_ok"))
        ok = ok and chain["ok"]
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
                    ("closed records without a finite threshold", bad_bound),
                    ("analytic status not matching its condition", status_bad),
                    ("conflicting repeated records", conflicts)]:
        if v:
            print(f"  {name}: {v}")
    if chain["checked"]:
        parts = [f"{k.split('_')[0]} margin {chain[k]:+.6f} ({'holds' if chain[k.replace('margin', 'ok')] else 'fails'})"
                 for k in ("general_margin", "sharper_margin") if k in chain]
        print(f"  stacks with n_2 = 1 (condition of the corollary): {'; '.join(parts)}")
    else:
        print("  stacks with n_2 = 1: condition not checked here (pass --prev-cert and/or --prev-ns)")
    print(f"  analytic statuses checked: EMPTY iff J empty, TIE iff t0 = l*; violations {len(status_bad)}, "
          f"conflicting repeated records {len(conflicts)}")
    print(f"  NON-SPINAL EXCLUSION {'CERTIFIED' if ok else 'NOT CERTIFIED'} above gamma_ns({a.N}) = {gamma_ns:.7f}"
          + ("" if chain["checked"] else " (branch coverage only; n_2 = 1 condition not checked)"))
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump({"N": a.N, "certified": ok, "gamma_ns": gamma_ns, "solver_gamma": solver_gamma,
                       "branches": len(expected), "open_main": len(open_main), "split_subbranches": len(want),
                       "n2_equals_1_condition": chain, "chain_checked": chain["checked"],
                       "analytic_status_violations": [list(map(str, x)) for x in status_bad],
                       "conflicting_repeated_records": {str(k): v for k, v in conflicts.items()}}, f, indent=2)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
