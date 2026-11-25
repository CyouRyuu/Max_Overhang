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
"""
import os
import csv
import json
import glob
import argparse
from datetime import datetime

import certify_run as cr
import profile_run as pr

CLOSED = ("CUTOFF", "INFEASIBLE")


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
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    lb = cr.KNOWN_LB[a.N] if a.lb is None else a.lb
    cutoff = lb * (1 + a.gap) if a.cutoff is None else a.cutoff
    rows = {}
    with open(a.csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows[r["profile"]] = {"status": r["status"], "bound": float(r["bound"]), "nodes": int(float(r["nodes"])),
                                  "runtime_sec": float(r["runtime_sec"]), "source": os.path.basename(a.csv)}
    for d in a.race:
        for p in glob.glob(os.path.join(d, "P*_result.json")):
            r = json.load(open(p))
            prof = r["profile"]
            if r["status"] in CLOSED or (prof not in rows or rows[prof]["status"] not in CLOSED):
                if prof not in rows or rows[prof]["status"] not in CLOSED:
                    rows[prof] = {"status": r["status"], "bound": r["bound"], "nodes": r["nodes"],
                                  "runtime_sec": r["runtime_sec"], "threads": r.get("threads"),
                                  "seed": r.get("seed"), "source": os.path.relpath(p)}
    for k, d in enumerate(a.split):
        prof = a.split_profile[k] if k < len(a.split_profile) else os.path.basename(os.path.normpath(d))
        br = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(d, "L*_result.json")))]
        if br and all(b["status"] in CLOSED for b in br) and (prof not in rows or rows[prof]["status"] not in CLOSED):
            rows[prof] = {"status": "CUTOFF", "bound": cutoff, "nodes": sum(b["nodes"] for b in br),
                          "runtime_sec": sum(b["runtime_sec"] for b in br),
                          "source": f"pd_split {os.path.relpath(d)} ({len(br)} branches)"}

    profiles = ["-".join(map(str, p)) for p in pr.enumerate_profiles(a.N, 2)]
    missing = [p for p in profiles if p not in rows]
    openp = [p for p in profiles if p in rows and not (rows[p]["status"] in CLOSED or
                                                        (rows[p]["status"] == "OPTIMAL" and rows[p]["bound"] <= cutoff))]
    bound = max((rows[p]["bound"] for p in profiles if p in rows), default=float("nan"))
    ok = not missing and not openp and bound <= cutoff + 1e-12
    nodes = sum(rows[p]["nodes"] for p in profiles if p in rows)
    cpu = sum(rows[p]["runtime_sec"] * (rows[p].get("threads") or 1) for p in profiles if p in rows)
    slow = sorted((p for p in profiles if p in rows), key=lambda p: -rows[p]["runtime_sec"])[:8]
    margin = lb - (cr.CERT_UB.get(a.N - 1, float("nan")) + 0.5 / a.N)
    res = {"N": a.N, "certified": ok, "lower_bound_feasible": lb, "cutoff_upper_bound": cutoff,
           "relative_gap": (cutoff - lb) / lb, "profiles_total": len(profiles), "profiles_closed": len(profiles) - len(missing) - len(openp),
           "missing": missing, "open": openp, "max_profile_bound": bound, "total_nodes": nodes,
           "total_thread_seconds": cpu, "S2plus_margin": margin, "S2plus_uses_CERT_UB_N_minus_1": cr.CERT_UB.get(a.N - 1),
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
             f"- S2+ (n_2 = 2) margin: {margin:.5f} using certified UB z*({a.N - 1}) <= {cr.CERT_UB.get(a.N - 1)}",
             f"- total nodes {nodes:,}; total thread-seconds {cpu:,.0f}", "",
             "| slowest profile | status | nodes | seconds | source |", "|---|---|---|---|---|"]
    lines += [f"| {p} | {rows[p]['status']} | {rows[p]['nodes']:,} | {rows[p]['runtime_sec']:.1f} | {rows[p]['source']} |"
              for p in slow]
    if missing or openp:
        lines += ["", f"missing: {missing}", f"open: {openp}"]
    with open(out + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
