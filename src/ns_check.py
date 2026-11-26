#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ns_check.py -- certify that NO non-spinal stack comes within a margin of the optimum.

The SP certificates show z*(N) <= S*(N) + 1e-5.  They do not exclude a non-spinal stack with overhang close to
S*(N).  This script closes, for every n_2 = 2 profile and every candidate rightmost layer l*, the branch
    max X_{i*} + W   s.t.  profile model (certificate options + SP cuts),
                          X_first(l) <= X_{i*} for all l, X_1 <= X_{i*}   (i* = first block of layer l* is a max block),
                          sum_{t in J_l*} sum_{q>=1} T_{s_t, i_{t-1,q}} >= 1  (the support chain of i* is NOT a pure spine),
                          X_{i*} >= gamma - W
with Cutoff gamma = S*(N) (1 - rel).  A non-spinal stack has every max block on an impure chain, so it lies in
the branch of the layer of any of its max blocks.  All branches closed  =>  every stack with overhang > gamma
is spinal, i.e. every stack within a factor (1 - rel) of the optimum (in particular every optimal stack) is spinal.
If J_l* is empty the chain of i* is always a pure spine, so the branch is closed without solving.

--split-t0 partitions a branch (n, l*) by t0 = the LOWEST t in J_l* at which s_t touches a non-rightmost block
of layer t-1.  Sub-branch t0 fixes T_{s_t, i_{t-1,q}} = 0 (q >= 1) for t in J, t < t0, requires
sum_q T_{s_t0, i_{t0-1,q}} >= 1, and adds two valid rows:
    X_{s_t0} <= X_{s_{t0-1}} - 1 + 2w     (s_t0 touches a block at least 1 left of s_{t0-1})
    X_{s_{t0-1}} + W <= 2w S^_{t0-1}(N)   (s_1..s_{t0-1} is a unique-support chain: spinal bound, paper: Spinal Bound proposition)
t0 = l* is closed without solving (status TIE): then X_{s_l*} <= X_{s_{l*-1}}, so the lowest max layer
(principal block, Paterson-Zwick 2009, Def. 2.4) is below l*, and its chain s_1..s_{l*-1} is pure, i.e. the stack is spinal.

  python3 ns_check.py 13 --rel 1e-2 --workers 7 --outdir results_ns
Results: <outdir>/ns_N<N>_rel<rel>.csv (one row per (profile, l*)) and a summary line.
"""
import os
import sys
import csv
import json
import time
import argparse
import multiprocessing as mp

import gurobipy as gp
from gurobipy import GRB

import profile_run as pr
import pd_split as ps

CLOSED = ("INFEASIBLE", "CUTOFF", "EMPTY", "TIE")


from profiles import spinal_opt  # noqa: E402


def solve_branch(job):
    N, profile, lstar, opts, gamma, threads, tlim, params = job[:8]
    tz = job[8] if len(job) > 8 else None          # first impure layer (None = unsplit branch)
    t0 = time.time()
    J = [t for t in range(3, lstar + 1) if profile[t - 2] >= 2]
    ptxt = "-".join(map(str, profile))
    base = {"profile": ptxt, "lstar": lstar, "t0": ("" if tz is None else tz)}
    if not J:
        return {**base, "status": "EMPTY", "nodes": 0, "bound": float("-inf"),
                "objective": float("nan"), "runtime_sec": 0.0}
    if tz is not None and tz == lstar:
        return {**base, "status": "TIE", "nodes": 0, "bound": float("-inf"),
                "objective": float("nan"), "runtime_sec": 0.0}
    m = ps.build_branch(N, profile, opts, gamma, lstar, True, True, False)
    if m is None:
        return {**base, "status": "INFEASIBLE", "nodes": 0, "bound": float("-inf"),
                "objective": float("nan"), "runtime_sec": time.time() - t0}
    V = m.getVarByName
    _, blocks = pr.layers_of(profile)
    if tz is None:
        m.addConstr(gp.quicksum(V(f"T_{blocks[t][0]}_{j}") for t in J for j in blocks[t - 1][1:]) >= 1,
                    name="NonSpinal")
    else:
        W = 0.5; w = 0.5                              # eta = 0, as in pd_split.build_branch
        for t in J:
            if t < tz:
                for j in blocks[t - 1][1:]:
                    V(f"T_{blocks[t][0]}_{j}").UB = 0.0
        m.addConstr(gp.quicksum(V(f"T_{blocks[tz][0]}_{j}") for j in blocks[tz - 1][1:]) >= 1,
                    name="NonSpinal_t0")
        s0, sp = blocks[tz][0], blocks[tz - 1][0]
        m.addConstr(V(f"X_{s0}") <= V(f"X_{sp}") - 1.0 + 2.0 * w, name="ImpureStep")
        m.addConstr(V(f"X_{sp}") + W <= 2.0 * w * pr.spinal_bound(tz - 1, N) + 1e-9, name="LowerSpine")
        m.update()
    m.setParam("OutputFlag", 0); m.setParam("NonConvex", 2); m.setParam("Threads", threads)
    m.setParam("MIPGap", 1e-4); m.setParam("TimeLimit", float(tlim)); m.setParam("Cutoff", gamma)
    for kv in params:
        k, v = kv.split("=", 1)
        m.setParam(k.strip(), float(v) if any(c in v for c in ".eE") else int(v))
    m.optimize()
    st = m.Status
    if st == GRB.INF_OR_UNBD:
        m.setParam("DualReductions", 0); m.reset(); m.optimize(); st = m.Status
    name = {GRB.INFEASIBLE: "INFEASIBLE", GRB.CUTOFF: "CUTOFF", GRB.OPTIMAL: "OPTIMAL",
            GRB.TIME_LIMIT: "TIME_LIMIT"}.get(st, f"STATUS_{st}")
    witness = None
    if m.SolCount > 0:                               # feasible stack found (only when a branch is not closed)
        witness = pr.base.stack_from_values(N, {v.VarName: v.X for v in m.getVars()})
    res = {**base, "status": name, "nodes": int(m.NodeCount), "witness": witness,
           "bound": (m.ObjBound if name not in CLOSED else gamma),
           "objective": (m.ObjVal if m.SolCount > 0 else float("nan")), "runtime_sec": time.time() - t0}
    m.dispose()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("--rel", type=float, default=1e-2, help="gamma = S*(N) (1 - rel)")
    ap.add_argument("--gamma", type=float, default=None, help="absolute gamma (overrides --rel)")
    ap.add_argument("--strengthen", default="TB,S1,S2P,S4,C2")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--time-limit", type=float, default=7200.0)
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--outdir", default="results_ns")
    ap.add_argument("--split-t0", action="store_true",
                    help="split every branch by the lowest impure layer t0 (see module docstring)")
    ap.add_argument("--only", default=None,
                    help="re-run only these branches: comma-separated profile:lstar, e.g. 1-2-2-3-2:5,1-2-3-3:4")
    a = ap.parse_args()
    opts = sorted({s.strip().upper() for s in a.strengthen.split(",") if s.strip()} | {"EN", "SP"})
    # solver threshold in double precision; it does not exceed gamma_ns(N) = 0.99 S*(N) rounded upward at the
    # 7th decimal (checked by ns_audit.py), so every exclusion proved here holds above gamma_ns(N)
    gamma = a.gamma if a.gamma is not None else spinal_opt(a.N) * (1.0 - a.rel)
    os.makedirs(a.outdir, exist_ok=True)
    tag = f"ns_N{a.N}_" + (f"g{a.gamma:.7f}" if a.gamma is not None else f"rel{a.rel:g}")
    jobs = []
    for p in pr.enumerate_profiles(a.N, 2):
        ok, _, _ = pr.energy_profile_check(a.N, p)
        if not ok:
            continue
        H = len(p)
        for l in range(2, H + 1):
            jobs.append((a.N, p, l, opts, gamma, a.threads, a.time_limit, ["Heuristics=0.0", "IntFeasTol=1e-9"] + a.param))
    if a.only:
        keep = {(x.split(":")[0].strip(), int(x.split(":")[1])) for x in a.only.split(",") if x.strip()}
        jobs = [j for j in jobs if ("-".join(map(str, j[1])), j[2]) in keep]
        if len(jobs) != len(keep):
            raise SystemExit(f"--only: matched {len(jobs)} of {len(keep)} branches")
        tag += "_only"
    if a.split_t0:
        sub = []
        for j in jobs:
            p, l = j[1], j[2]
            Jl = [t for t in range(3, l + 1) if p[t - 2] >= 2]
            sub.extend([j + (t,) for t in Jl] if Jl else [j])
        jobs = sub
        tag += "_t0"
    t0 = time.time()
    rows = []
    path = os.path.join(a.outdir, tag + ".csv")
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["profile", "lstar", "t0", "status", "objective", "bound", "nodes", "runtime_sec"])
        wr.writeheader()
        it = (map(solve_branch, jobs) if a.workers <= 1 else
              mp.get_context("spawn").Pool(a.workers, maxtasksperchild=20).imap_unordered(solve_branch, jobs))
        witnesses = []
        for r in it:
            wit = r.pop("witness", None)
            if wit is not None:
                witnesses.append({"profile": r["profile"], "lstar": r["lstar"], "t0": r["t0"],
                                  "overhang": r["objective"], "stack": wit})
            rows.append(r); wr.writerow(r); f.flush()
    if witnesses:                                    # feasible non-spinal stacks above gamma (reverse controls)
        witnesses.sort(key=lambda x: -x["overhang"])
        with open(os.path.join(a.outdir, tag + "_witnesses.json"), "w") as fw:
            json.dump({"N": a.N, "gamma": gamma, "witnesses": witnesses}, fw, indent=1)
    bad = [r for r in rows if r["status"] not in CLOSED]
    best = max([r["objective"] for r in rows if r["objective"] == r["objective"]], default=float("nan"))
    cpu = sum(r["runtime_sec"] for r in rows)
    line = (f"N={a.N}" + (f" ({len(a.only.split(','))} selected branches, time limit {a.time_limit:.0f}s, {a.threads} thread(s))"
               if a.only else "")
            + (" (split by lowest impure layer t0)" if a.split_t0 else "")
            + f" gamma={gamma:.12f} (S*={spinal_opt(a.N):.7f}) branches={len(rows)} open={len(bad)} "
            f"nodes={sum(r['nodes'] for r in rows)} thread_s={cpu:.0f} wall={time.time()-t0:.0f}s "
            f"best_nonspinal_found={best:.7f} => " + ("NO non-spinal stack above gamma" if not bad else "NOT CERTIFIED"))
    print(line)
    with open(os.path.join(a.outdir, "ns_summary.txt"), "a") as f:
        f.write(line + "\n")


if __name__ == "__main__":
    main()
