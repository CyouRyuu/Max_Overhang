#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pd_certificate.py -- assemble and check a layer-profile-decomposition certificate.

Reads the profile CSV of profile_run.py plus the result JSON files of pd_race.py / pd_split.py runs, checks
that EVERY n_2 = 2 profile of N blocks is closed (CUTOFF / INFEASIBLE, or OPTIMAL with bound <= cutoff,
or all rightmost-layer branches of pd_split.py closed), and writes <out>.json and a readable <out>.md.

  python3 pd_certificate.py 13 results_pd13/<tag>_profiles.csv \
      --race results_pd13_race --split results_pd13_split/1-2-2-2-3-2-1 --out certificate_N13
(--split DIR may be repeated; DIR must contain the L<l>_result.json files of one profile, and the profile is
 taken from the directory name or from --split-profile.)

Checks: every expected profile has a record and no record lies outside P_2(N); repeated records are reported
(the last one is used); every closed record (CUTOFF) carries a bound <= the certificate cutoff; the
cutoff stored in the run's *_result.json (required unless --allow-missing-run-record) is <= the certificate
cutoff; a profile closed by --split has a closed L<l>_result.json for every candidate layer l (layers whose
rightmost block can exceed the cutoff by the base envelope and the energy layer bounds); with --prev-cert, the
S2+ condition zbar(N-1) + 1/(2N) < lower bound against the audited certificate of N - 1 (part of the pass
condition).  Profiles closed without a solver run are listed with their analytic reason (energy test, or no
candidate layer), and each such record (no model built: vars = 0) is checked against the analytic reason
recomputed here, in both directions; repeated records with conflicting closed statuses fail the check.  Needs no
solver.
"""
import os
import math
import sys
import csv
import json
import glob
import argparse
from datetime import datetime

import profiles as pf

CLOSED = ("CUTOFF", "INFEASIBLE")


def split_candidate_layers(N, profile, cutoff, eta=0.0, energy=True):
    """Layers l >= 2 whose rightmost block can exceed the cutoff: X_first(l) + W <= x1ub + W + 2w (l - 1) with the
    tightened base envelope x1ub = (N - 3) w / N - eps, and X_first(l) + W <= energy layer bound.  Every stack of
    the profile above the cutoff has its principal block in one of these layers."""
    W = 0.5; w = W * (1.0 - eta); eps = W * eta
    x1ub = (N - 3) * w / N - eps
    ok, _, ubs = pf.energy_profile_check(N, profile, eta)
    if not energy:
        ok, ubs = True, {}
    if not ok:
        return []
    return [l for l in range(2, len(profile) + 1)
            if x1ub + W + 2 * w * (l - 1) > cutoff and (l not in ubs or ubs[l] > cutoff)]


def closure_reason(N, profile, cutoff, spine, eta=0.0, energy=True):
    """Analytic reason why a profile needs no solver run, or None.  'energy': the second-moment energy test
    rejects the profile (Delta E < 0); 'no candidate layer': no layer can hold a principal block above the cutoff
    (base envelope, energy layer bounds and, with family SP, the pure-chain bound 2w S^_l(N) of the spinal cuts)."""
    ok, dE, _ = pf.energy_profile_check(N, profile, eta)
    if energy and not ok:
        return "energy", dE
    W = 0.5; w = W * (1.0 - eta)
    cands = split_candidate_layers(N, profile, cutoff, eta, energy)
    if spine:
        cands = [l for l in cands if [t for t in range(3, l + 1) if profile[t - 2] >= 2]
                 or 2.0 * w * pf.spinal_bound(l, N) + 1e-9 > cutoff]
    if not cands:
        return "no candidate layer", None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("csv")
    ap.add_argument("--race", action="append", default=[], help="directory with P<profile>_s<seed>_result.json")
    ap.add_argument("--split", action="append", default=[], help="pd_split.py output directory of one profile")
    ap.add_argument("--split-profile", action="append", default=[], help="profile for each --split (optional)")
    ap.add_argument("--gap", type=float, default=5e-4)
    ap.add_argument("--cutoff", type=float, default=None, help="absolute cutoff used by the runs (overrides gap)")
    ap.add_argument("--lb", type=float, default=None, help="exactly verified lower bound (default KNOWN_LB)")
    ap.add_argument("--allow-missing-run-record", action="store_true",
                    help="accept a profile CSV without its *_result.json (e.g. a run stopped by hand)")
    ap.add_argument("--prev-cert", default=None,
                    help="certificate JSON of N - 1 (this script's output) or the M0 anchor result JSON for N = 5; "
                         "its upper bound enters the S2+ condition zbar(N-1) + 1/(2N) < lower bound, which must hold")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    lb = pf.KNOWN_LB[a.N] if a.lb is None else a.lb
    cutoff = lb * (1 + a.gap) if a.cutoff is None else a.cutoff
    rows, seen, hist, csv_meta = {}, {}, {}, {}
    with open(a.csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            # a run continued with --resume may record a profile again; the last record is the one used
            seen[r["profile"]] = seen.get(r["profile"], 0) + 1
            hist.setdefault(r["profile"], []).append(r["status"])
            csv_meta[r["profile"]] = (r["status"], (r.get("vars") or "").strip(), (r.get("closure") or "").strip())
            rows[r["profile"]] = {"status": r["status"], "bound": float(r["bound"]), "nodes": int(float(r["nodes"])),
                                  "runtime_sec": float(r["runtime_sec"]), "source": os.path.basename(a.csv)}
    duplicates = {p: c for p, c in seen.items() if c > 1}
    # repeated records conflict if two of them are closed with different statuses, or if a closed record is
    # followed by an open one (a resumed run only re-solves profiles that were open)
    conflicts = {p: h for p, h in hist.items() if len(h) > 1 and
                 (len({x for x in h if x in CLOSED}) > 1 or (h[-1] not in CLOSED and any(x in CLOSED for x in h)))}
    # the cutoff actually used by the run is stored in its *_result.json next to the CSV
    run_json = a.csv.replace("_profiles.csv", "_result.json")
    run_rec = json.load(open(run_json)) if os.path.exists(run_json) else {}
    run_cutoff = run_rec.get("cutoff")
    spine = "'SP'" in str(run_rec.get("method", ""))
    # the energy test runs only with option EN (in the run tag, e.g. ...-RB-EN-...)
    energy = "-EN-" in os.path.basename(a.csv) or "'EN'" in str(run_rec.get("method", ""))
    for d in a.race:
        for p in glob.glob(os.path.join(d, "P*_result.json")):
            r = json.load(open(p))
            prof = r["profile"]
            if prof not in rows or rows[prof]["status"] not in CLOSED:
                rows[prof] = {"status": r["status"], "bound": r["bound"], "nodes": r["nodes"],
                              "runtime_sec": r["runtime_sec"], "threads": r.get("threads"),
                              "seed": r.get("seed"), "source": os.path.relpath(p)}
    split_issues = {}
    for k, d in enumerate(a.split):
        prof = a.split_profile[k] if k < len(a.split_profile) else os.path.basename(os.path.normpath(d))
        if prof.startswith("split_"):
            prof = prof[len("split_"):]
        br = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(d, "L*_result.json")))]
        need = split_candidate_layers(a.N, [int(x) for x in prof.split("-")], cutoff)
        have = {b["layer"]: b for b in br}
        bad = sorted(l for l in need if l not in have or have[l]["status"] not in CLOSED
                     or have[l]["bound"] > cutoff + 1e-12)
        if bad:
            split_issues[prof] = bad
            continue
        if br and (prof not in rows or rows[prof]["status"] not in CLOSED):
            rows[prof] = {"status": "CUTOFF", "bound": cutoff, "nodes": sum(b["nodes"] for b in br),
                          "runtime_sec": sum(b["runtime_sec"] for b in br),
                          "source": f"pd_split {os.path.relpath(d)} ({len(br)} branches)"}

    profiles = ["-".join(map(str, p)) for p in pf.enumerate_profiles(a.N, 2)]
    missing = [p for p in profiles if p not in rows]
    unexpected = sorted(p for p in rows if p not in set(profiles))
    if run_cutoff is None:
        cutoff_ok = a.allow_missing_run_record
    else:
        cutoff_ok = run_cutoff <= cutoff + 1e-12
    above = sorted(p for p in profiles if p in rows and rows[p]["status"] in CLOSED
                   and rows[p]["bound"] > cutoff + 1e-12)
    openp = [p for p in profiles if p in rows and not (rows[p]["status"] in CLOSED or
                                                        (rows[p]["status"] == "OPTIMAL" and rows[p]["bound"] <= cutoff))]
    bound = max((rows[p]["bound"] for p in profiles if p in rows), default=float("nan"))
    ok = (not missing and not openp and not unexpected and not above and not split_issues and cutoff_ok
          and bound <= cutoff + 1e-12 and not conflicts)
    nodes = sum(rows[p]["nodes"] for p in profiles if p in rows)
    cpu = sum(rows[p]["runtime_sec"] * (rows[p].get("threads") or 1) for p in profiles if p in rows)
    slow = sorted((p for p in profiles if p in rows), key=lambda p: -rows[p]["runtime_sec"])[:8]
    if a.prev_cert:
        pc = json.load(open(a.prev_cert))
        if "certified" in pc:
            prev_ub, prev_ok = pc["cutoff_upper_bound"], bool(pc["certified"]) and pc.get("N") == a.N - 1
        else:   # M0 anchor result (N = 4): closed run, upper bound = final dual bound rounded up at the 7th decimal
            prev_ub = math.ceil(pc["bound"] * 1e7 - 1e-6) / 1e7
            prev_ok = pc.get("status") == "OPTIMAL" and pc.get("N", a.N - 1) == a.N - 1
        prev_src = os.path.basename(a.prev_cert)
    else:
        prev_ub, prev_ok, prev_src = pf.CERT_UB.get(a.N - 1, float("nan")), None, "table CERT_UB (not re-audited)"
    margin = lb - (prev_ub + 0.5 / a.N)
    s2p_ok = margin > 0 and prev_ok is not False
    ok = ok and s2p_ok
    reasons = {}
    for p in profiles:
        if p in rows and rows[p]["status"] in CLOSED:
            r = closure_reason(a.N, tuple(int(x) for x in p.split("-")), cutoff, spine, energy=energy)
            if r:
                reasons[p] = r
    # records written without a solver run (vars = 0) must match the recomputed analytic reason, and conversely
    analytic_mismatch = []
    for p, (st, nv, cl) in csv_meta.items():
        if p not in set(profiles) or nv == "":
            continue
        r = closure_reason(a.N, tuple(int(x) for x in p.split("-")), cutoff, spine, energy=energy)
        want = None if r is None else ("INFEASIBLE" if r[0] == "energy" else "CUTOFF")
        if (nv == "0") != (r is not None) or (r is not None and st != want) or \
                (cl.startswith("analytic") and r is None):
            analytic_mismatch.append((p, st, nv, cl, None if r is None else r[0]))
    ok = ok and not analytic_mismatch
    n_energy = sum(1 for r in reasons.values() if r[0] == "energy")
    n_empty = sum(1 for r in reasons.values() if r[0] == "no candidate layer")
    worst_dE = max((r[1] for r in reasons.values() if r[0] == "energy"), default=None)
    res = {"N": a.N, "certified": ok, "lower_bound_feasible": lb, "cutoff_upper_bound": cutoff,
           "relative_gap": (cutoff - lb) / lb, "profiles_total": len(profiles), "profiles_closed": len(profiles) - len(missing) - len(openp),
           "missing": missing, "open": openp, "unexpected": unexpected, "duplicates_last_record_used": duplicates,
           "run_cutoff": run_cutoff, "run_cutoff_within_certificate": cutoff_ok,
           "closed_above_cutoff": above, "split_missing_or_open_layers": split_issues, "max_profile_bound": bound, "total_nodes": nodes,
           "total_thread_seconds": cpu, "S2plus_margin": margin, "S2plus_ok": s2p_ok,
           "S2plus_upper_bound_N_minus_1": prev_ub, "S2plus_upper_bound_source": prev_src,
           "S2plus_previous_certified": prev_ok, "family_SP": spine,
           "chain_checked": bool(a.prev_cert), "conflicting_repeated_records": conflicts,
           "analytic_records_checked": sum(1 for v in csv_meta.values() if v[1] != ""),
           "analytic_record_mismatches": analytic_mismatch,
           "analytic_closures": {"energy": n_energy, "no_candidate_layer": n_empty,
                                 "largest_energy_margin": (str(worst_dE) if worst_dE is not None else None),
                                 "profiles": {p: r[0] for p, r in reasons.items()}},
           "solver_closures": sum(1 for p in profiles if p in rows and rows[p]["status"] in CLOSED) - len(reasons),
           "status_counts": {s: sum(1 for p in profiles if p in rows and rows[p]["status"] == s)
                             for s in sorted({r["status"] for r in rows.values()})},
           "slowest": [{"profile": p, **rows[p]} for p in slow], "created": f"{datetime.now():%F %T}"}
    out = a.out or f"certificate_N{a.N}"
    with open(out + ".json", "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    lines = [f"# Layer-profile certificate, N = {a.N}", "",
             f"- certified: **{ok}**",
             f"- feasible value (lower bound): {lb:.7f}",
             f"- upper bound (cutoff): {cutoff:.8f}  (relative gap {res['relative_gap']:.2e})",
             f"- profiles: {res['profiles_closed']} / {len(profiles)} closed; status counts {res['status_counts']}",
             f"- missing {len(missing)}, unexpected {len(unexpected)}, repeated records {len(duplicates)} (last record used)",
             ("- cutoff of the run: record missing (accepted with --allow-missing-run-record)" if run_cutoff is None
              and cutoff_ok else f"- cutoff of the run: {run_cutoff} "
              f"({'<= certificate cutoff' if cutoff_ok else 'NOT <= certificate cutoff or record missing'})"),
             f"- closed records with bound above the cutoff: {len(above)}; split profiles with missing/open layers: "
             f"{len(split_issues)}",
             f"- S2+ (n_2 = 2) condition zbar({a.N - 1}) + 1/(2N) < lower bound: margin {margin:.5f} "
             f"({'holds' if s2p_ok else 'FAILS'}); zbar({a.N - 1}) = {prev_ub} from {prev_src}"
             + ("" if prev_ok is None else f", certified: {prev_ok}"),
             ("- S2+ precondition checked against the audited certificate of N - 1" if a.prev_cert else
              "- S2+ precondition checked against the table CERT_UB only (coverage audit; pass --prev-cert to check the chain)"),
             f"- closures without a solver run: {n_energy} by the energy test (largest Delta E {worst_dE}), "
             f"{n_empty} with no candidate layer; {res['solver_closures']} closed by the solver",
             f"- records without a solver run checked against the recomputed analytic reason: "
             f"{res['analytic_records_checked']} records, {len(analytic_mismatch)} mismatches; conflicting repeated "
             f"records: {len(conflicts)}",
             f"- total nodes {nodes:,}; total thread-seconds {cpu:,.0f}", "",
             "| slowest profile | status | nodes | seconds | source |", "|---|---|---|---|---|"]
    lines += [f"| {p} | {rows[p]['status']} | {rows[p]['nodes']:,} | {rows[p]['runtime_sec']:.1f} | {rows[p]['source']} |"
              for p in slow]
    if missing or openp or unexpected or duplicates or above or split_issues or not cutoff_ok or analytic_mismatch:
        lines += ["", f"missing: {missing}", f"open: {openp}", f"unexpected: {unexpected}", f"repeated: {duplicates}",
                  f"conflicting repeated records: {conflicts}", f"analytic record mismatches: {analytic_mismatch}",
                  f"closed above cutoff: {above}", f"split layers missing or open: {split_issues}",
                  f"run record: {'cutoff ' + str(run_cutoff) if run_cutoff is not None else 'missing'}"]
    with open(out + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
