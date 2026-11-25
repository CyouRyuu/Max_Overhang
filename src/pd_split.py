#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pd_split.py -- rightmost-layer split of ONE hard layer profile (independent cross-check in the paper's methodological-progression subsection).

For a profile n and cutoff gamma = LB (1 + gap), every stack of profile n with overhang > gamma has its
rightmost block at the first block i* = first(l*) of some candidate layer l* (layers 2..H whose first block
can exceed gamma).  One sub-subproblem per l*:
    max X_{i*} + W   s.t.  profile model (same options as profile_run.py),
                           X_{first(l)} <= X_{i*}  (l != l*),   X_1 <= X_{i*},
                           X_{i*} >= gamma - W                 (only stacks that beat the cutoff),
                           X_{first(l)} >= (gamma - W) - 2w (l* - l)  for l < l*   (support chain of i*),
solved with Cutoff gamma.  INFEASIBLE or CUTOFF closes the branch.  All branches closed => the profile cannot
beat gamma, which is exactly what profile_run.py needs for that profile.

Phase 1: every branch with 1 thread for --phase1-time seconds (the easy ones close here).
Phase 2: the open branches share --threads-budget threads (budget // #open threads each), no restart limit.
Results: <outdir>/L<l*>_result.json (+ Gurobi logs); finished branches are skipped on a re-run.

Example (N = 13, the profile containing the optimum, certificate options without SP):
  nohup nice -n 10 python3 pd_split.py 13 1-2-2-2-3-2-1 --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer \
      --param Heuristics=0.0 --threads-budget 7 --outdir results_pd13_split > pd_split.log 2>&1 &
"""
import os
import sys
import json
import time
import argparse
import multiprocessing as mp
from datetime import datetime

import gurobipy as gp
from gurobipy import GRB

import certify_run as cr
import profile_run as pr

CLOSED = ("INFEASIBLE", "CUTOFF")


def build_branch(N, profile, opts, cutoff, lstar, rank_box, c2_layer, tail):
    W = 0.5; w = W                      # eta = 0
    model, info = pr.build_profile_model(N, profile, "M1u", set(opts), 0.0, cutoff, tail, lambda s: None,
                                         rank_box=rank_box, c2_layer=c2_layer)
    V = model.getVarByName
    _, blocks = pr.layers_of(profile)
    H = len(profile)
    for g in list(model.getGenConstrs()):
        model.remove(g)
    istar = blocks[lstar][0]
    xs = V(f"X_{istar}")
    for l in range(2, H + 1):
        if l != lstar:
            model.addConstr(V(f"X_{blocks[l][0]}") <= xs, name=f"RM_dom_{l}")
    model.addConstr(V("X_1") <= xs, name="RM_dom_1")
    model.addConstr(xs >= cutoff - W, name="RM_beats_cutoff")
    xs.LB = max(xs.LB, cutoff - W)
    model.update()
    for l in range(1, lstar):
        v = V(f"X_{blocks[l][0]}")
        lb = (cutoff - W) - 2 * w * (lstar - l)
        if lb > v.UB + 1e-9:
            return None                  # pruned by bounds alone
        v.LB = max(v.LB, lb)
    model.update()
    if rank_box:
        pr.tighten_contact_boxes(model, profile, blocks, V, w)
    model.setObjective(xs + W, GRB.MAXIMIZE)
    model.update()
    return model


def solve_branch(job):
    (N, profile, opts, cutoff, lstar, rank_box, c2_layer, tail, threads, tlim, seed, params, outdir, phase) = job
    t0 = time.time()
    model = build_branch(N, profile, opts, cutoff, lstar, rank_box, c2_layer, tail)
    if model is None:
        return {"layer": lstar, "status": "INFEASIBLE", "note": "bounds", "nodes": 0, "bound": float("-inf"),
                "runtime_sec": time.time() - t0, "phase": phase, "threads": threads}
    model.setParam("OutputFlag", 1); model.setParam("LogToConsole", 0)
    model.setParam("LogFile", os.path.join(outdir, f"L{lstar}_phase{phase}.log"))
    model.setParam("NonConvex", 2); model.setParam("Threads", threads); model.setParam("Seed", seed)
    model.setParam("MIPGap", 5e-4); model.setParam("TimeLimit", float(tlim)); model.setParam("Cutoff", cutoff)
    for kv in params:
        k, v = kv.split("=", 1); v = v.strip()
        try:
            val = int(v)
        except ValueError:
            try:
                val = float(v)
            except ValueError:
                val = v
        model.setParam(k.strip(), val)
    model.optimize()
    st = model.Status
    if st == GRB.INF_OR_UNBD:
        model.setParam("DualReductions", 0); model.reset(); model.optimize(); st = model.Status
    name = {GRB.INFEASIBLE: "INFEASIBLE", GRB.CUTOFF: "CUTOFF", GRB.OPTIMAL: "OPTIMAL",
            GRB.TIME_LIMIT: "TIME_LIMIT"}.get(st, f"STATUS_{st}")
    res = {"layer": lstar, "status": name, "nodes": int(model.NodeCount),
           "bound": (model.ObjBound if name not in CLOSED else cutoff), "objective":
           (model.ObjVal if model.SolCount > 0 else None), "runtime_sec": time.time() - t0,
           "phase": phase, "threads": threads}
    if model.SolCount > 0:
        res["stack_X"] = {v.VarName: v.X for v in model.getVars() if v.VarName.startswith("X_")}
    model.dispose()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("profile", help="e.g. 1-2-2-2-3-2-1")
    ap.add_argument("--strengthen", default="TB,S1,S2P,S4,C2")
    ap.add_argument("--rank-box", action="store_true")
    ap.add_argument("--c2-layer", action="store_true")
    ap.add_argument("--tail", action="store_true")
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--gap", type=float, default=5e-4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads-budget", type=int, default=4)
    ap.add_argument("--phase1-time", type=float, default=1800.0)
    ap.add_argument("--time-limit", type=float, default=1209600.0, help="phase-2 limit per branch (s)")
    ap.add_argument("--nodefile-start", type=float, default=1.5)
    ap.add_argument("--soft-mem-limit", type=float, default=6.0)
    ap.add_argument("--cutoff", type=float, default=None,
                    help="objective cutoff (default: KNOWN_LB[N] * (1 + gap))")
    ap.add_argument("--outdir", default="results_split")
    a = ap.parse_args()

    profile = tuple(int(x) for x in a.profile.split("-"))
    assert sum(profile) == a.N and profile[:2] == (1, 2), "profile must be an n_2 = 2 profile of N blocks"
    opts = sorted({s.strip().upper() for s in a.strengthen.split(",") if s.strip()})
    cutoff = a.cutoff if a.cutoff is not None else cr.KNOWN_LB[a.N] * (1.0 + a.gap)
    os.makedirs(a.outdir, exist_ok=True)
    nodedir = os.path.join(a.outdir, "nodefiles"); os.makedirs(nodedir, exist_ok=True)
    params = [f"NodefileStart={a.nodefile_start}", f"SoftMemLimit={a.soft_mem_limit}",
              f"NodefileDir={nodedir}"] + list(a.param)

    # candidate layers (same filter as profile_run.py) and the energy layer bound (paper: energy corollary)
    _, blocks = pr.layers_of(profile)
    m0, info = pr.build_profile_model(a.N, profile, "M1u", set(opts), 0.0, cutoff, a.tail, lambda s: None,
                                      rank_box=a.rank_box, c2_layer=a.c2_layer)
    cand_layers = [l for l in range(2, len(profile) + 1) if blocks[l][0] in info["candidates"]]
    m0.dispose()
    ok, dE, ubs = pr.energy_profile_check(a.N, profile)
    print(f"[{datetime.now():%F %T}] N={a.N} profile {a.profile}  cutoff {cutoff:.8f}  candidate layers {cand_layers}"
          f"  energy dE={dE:.3f}  layer bounds " + ", ".join(f"{l}:{ubs.get(l, float('nan')):.4f}" for l in cand_layers),
          flush=True)

    def done(l):
        p = os.path.join(a.outdir, f"L{l}_result.json")
        if os.path.exists(p):
            r = json.load(open(p))
            return r["status"] in CLOSED or r["status"] == "OPTIMAL"
        return False

    def record(r):
        with open(os.path.join(a.outdir, f"L{r['layer']}_result.json"), "w", encoding="utf-8") as f:
            json.dump(r, f, indent=2)
        print(f"[{datetime.now():%F %T}]   layer {r['layer']} (first block {blocks[r['layer']][0]}): {r['status']}"
              f"  bound={r['bound']}  nodes={r['nodes']}  t={r['runtime_sec']:.1f}s  phase {r['phase']}"
              f"  threads {r['threads']}" + (f"  NEW SOLUTION {r['objective']:.7f}" if r.get("objective") else ""),
              flush=True)

    todo = [l for l in cand_layers if not done(l) and not (l in ubs and ubs[l] <= cutoff)]
    for l in cand_layers:
        if l in ubs and ubs[l] <= cutoff and not done(l):
            record({"layer": l, "status": "INFEASIBLE", "note": "energy bound", "nodes": 0, "bound": ubs[l],
                    "runtime_sec": 0.0, "phase": 0, "threads": 0})
    ctx = mp.get_context("spawn")
    for phase in (1, 2):
        todo = [l for l in todo if not done(l)]
        if not todo:
            break
        thr = 1 if phase == 1 else max(1, a.threads_budget // len(todo))
        wk = min(len(todo), a.threads_budget) if phase == 1 else len(todo)
        tl = a.phase1_time if phase == 1 else a.time_limit
        print(f"[{datetime.now():%F %T}] phase {phase}: layers {todo}, {wk} workers x {thr} threads, limit {tl:.0f}s",
              flush=True)
        jobs = [(a.N, profile, opts, cutoff, l, a.rank_box, a.c2_layer, a.tail, thr, tl, a.seed, params,
                 a.outdir, phase) for l in todo]
        with ctx.Pool(wk, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(solve_branch, jobs):
                if r["status"] in CLOSED or r["status"] == "OPTIMAL" or phase == 2:
                    record(r)
                else:
                    print(f"[{datetime.now():%F %T}]   layer {r['layer']}: still open after phase 1"
                          f" (bound {r['bound']:.7f}, {r['nodes']} nodes) -> phase 2", flush=True)
    closed = all(done(l) and json.load(open(os.path.join(a.outdir, f"L{l}_result.json")))["status"] in CLOSED
                 for l in cand_layers)
    print(f"[{datetime.now():%F %T}] profile {a.profile}: " +
          ("ALL BRANCHES CLOSED -> no stack with this profile beats the cutoff" if closed else "not closed"),
          flush=True)


if __name__ == "__main__":
    main()
