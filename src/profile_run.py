#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
profile_run.py -- layer-profile decomposition (paper: Layer-Profile Decomposition theorem).

Every stack has a layer profile (n_1, ..., n_H) with n_1 = 1.  The layer-width doubling theorem gives
n_{l+1} <= 2 n_l and n_l <= 2l - 2, and the base-singleton reduction (S2+) gives n_2 = 2 at every optimum
for N = 5..16 (it uses the certified bound for N - 1, see certify_run.CERT_UB).  Blocks are numbered
layer by layer, so a profile fixes every layer index L_i.  One subproblem is solved per profile:

  * the monolithic model (M1u or Mrel) with --strengthen rows (default S1,S2P,C2 and TB),
  * L_i and LayerJump_i fixed by the profile (contacts between non-adjacent layers drop out),
  * objective restricted to the first block of each layer l >= 2 (the rightmost block of its layer),
  * planar non-crossing rows  T[i1,j2] + T[i2,j1] <= 1  for adjacent layers, (i2-i1)+(j2-j1) >= 3,
  * every singleton layer {u}: the blocks of the next layer rest on u only (T = 1); with one block v
    above, F_vu = M (blocks above the layer of u); with two blocks v1 < v2 above, the generalised C2
    rows (difference RLT on F_v1u * (P_v1u - P_v2u), or the Mrel form) and the sign rows
    X_v1 >= X_u >= X_v2,
  * optional --tail: harmonic canonical form of a single-column top tail (S3; needs a singleton layer
    directly below the tail, tail length m >= 2).

Every subproblem gets Cutoff = LB (1 + gap) with LB a known feasible value.  Status CUTOFF (or
INFEASIBLE) proves that no stack with that profile beats the cutoff.  The certified bound is the
largest subproblem bound; the result is written in the certify_run.py JSON format.

Usage
  python3 profile_run.py 11 --model M1u --strengthen TB,S1,S2P,C2 --seed 42 --workers 1 --outdir results_profile
  python3 profile_run.py 13 --model M1u --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --workers 8 \
          --time-limit 86400 --keep-logs --outdir results_profile13      (add --resume to continue a run)
"""

import os
import sys
import csv
import json
import time
import socket
import argparse
import multiprocessing as mp
from datetime import datetime

import gurobipy as gp
from gurobipy import GRB

import certify_run as cr
from profiles import enumerate_profiles, layers_of, energy_profile_check, spinal_bound  # noqa: F401
import certify_m0 as base


def add_rank_box(model, N, profile, blocks, V, w, W):
    """Layer-rank box (paper: profile-specific bounds, part (i)).  For the block of rank r (0 = rightmost)
    in layer l:  X_1 - 2w(l-1) + (n_l - 1 - r) <= X_i <= X_1 + 2w(l-1) - r  (support chains to block 1 and
    same-layer spacing >= 1).  Also tightens the contact-point and torque boxes of every adjacent pair."""
    X1 = V("X_1"); x1lb, x1ub = max(X1.LB, -W), X1.UB
    H = len(profile); n = 0
    for l in range(2, H + 1):
        nl = profile[l - 1]
        for r, i in enumerate(blocks[l]):
            up = 2 * w * (l - 1) - r; dn = 2 * w * (l - 1) - (nl - 1 - r)
            xi = V(f"X_{i}")
            xi.LB = max(xi.LB, x1lb - dn); xi.UB = min(xi.UB, x1ub + up)
            model.addConstr(xi - X1 <= up, name=f"RB_up_{i}"); model.addConstr(X1 - xi <= dn, name=f"RB_dn_{i}")
            n += 2
    model.update()
    tighten_contact_boxes(model, profile, blocks, V, w)
    return n


def tighten_contact_boxes(model, profile, blocks, V, w):
    """Recompute the P and Torque boxes of every adjacent pair from the current X / F / T bounds
    (keep 0 in the box while T may be 0)."""
    H = len(profile)
    for l in range(1, H):
        for i in blocks[l + 1]:
            for j in blocks[l]:
                T, F, P, tq = V(f"T_{i}_{j}"), V(f"F_{i}_{j}"), V(f"P_{i}_{j}"), V(f"Torque_{i}_{j}")
                if P is None:
                    continue
                xi, xj = V(f"X_{i}"), V(f"X_{j}")
                plo = max(xi.LB, xj.LB) - w; phi = min(xi.UB, xj.UB) + w
                if T.LB < 0.5:          # T may be 0, then P = 0
                    plo, phi = min(plo, 0.0), max(phi, 0.0)
                plo = max(P.LB, plo); phi = min(P.UB, phi)
                if plo > phi:
                    continue            # contact impossible; leave it to presolve
                P.LB, P.UB = plo, phi
                cands = [a * b for a in (F.LB, F.UB) for b in (plo, phi)] + ([0.0] if T.LB < 0.5 else [])
                tq.LB = max(tq.LB, min(cands)); tq.UB = min(tq.UB, max(cands))
    model.update()


def add_c2_layer_cuts(model, profile, blocks, V, rel, w):
    """Family C2-Layer (paper: two-block layer cut-set RLT cuts).  For every non-top layer l with
    n_l = 2 (blocks u1 < u2): the blocks above weigh M, F1 = sum_k F_{k,u1}, Delta = X_u1 - X_u2 in
    [1, Dbar] with Dbar = min(4w (l - s0), box), s0 = highest singleton layer below l, and
      sum_{b above} X_b - M X_u2 - sum_k [mu_top(k,u1) + mu_top(k,u2)] = F1 * Delta   (exact identity).
    M1u has no mu_top variables; |sum mu_top| <= M w gives the same 4 McCormick rows with slack M w."""
    H = len(profile); n = 0
    for l in range(2, H):
        if profile[l - 1] != 2:
            continue
        u1, u2 = blocks[l]
        upper = [b for t in range(l + 1, H + 1) for b in blocks[t]]
        M = float(len(upper)); nxt = blocks[l + 1]
        s0 = max(t for t in range(1, l) if profile[t - 1] == 1)
        dub = min(4.0 * w * (l - s0), V(f"X_{u1}").UB - V(f"X_{u2}").LB)
        if dub < 1.0 - 1e-9:
            raise RuntimeError(f"profile {profile}: layer {l} cannot hold two blocks")
        dub = max(dub, 1.0)
        Delta = V(f"X_{u1}") - V(f"X_{u2}")
        model.addConstr(Delta <= dub, name=f"C2L_DeltaUB_{l}")
        F1 = gp.quicksum(V(f"F_{k}_{u1}") for k in nxt)
        Wexp = gp.quicksum(V(f"X_{b}") for b in upper) - M * V(f"X_{u2}")
        if rel:
            Wexp = Wexp - gp.quicksum(V(f"MuTop_{k}_{u1}") + V(f"MuTop_{k}_{u2}") for k in nxt)
            sl = 0.0
        else:
            sl = M * w
        model.addConstr(Wexp >= F1 - sl, name=f"C2L_a_{l}")
        model.addConstr(Wexp >= M * Delta + dub * F1 - M * dub - sl, name=f"C2L_b_{l}")
        model.addConstr(Wexp <= M * Delta + F1 - M + sl, name=f"C2L_c_{l}")
        model.addConstr(Wexp <= dub * F1 + sl, name=f"C2L_d_{l}")
        n += 1
    model.update()
    return n


def add_spinal_cuts(model, N, profile, blocks, V, w, W, cands, cutoff):
    """SP: for every candidate rightmost block s_l (l >= 2):
         X_{s_l} + W <= 2w S^_l(N) + M_l * sum_{t in J_l} sum_{q >= 1} T_{s_t, i_{t-1,q}},
       J_l = {t in 3..l : n_{t-1} >= 2}.  If the sum is 0, every s_t (t <= l) rests only on s_{t-1}
       (n_1 = 1 covers t = 2), so the stack is spinal up to s_l and the spinal bound (paper: Spinal Bound proposition) applies.  Valid for every
       feasible stack (no cutoff needed).  Returns (#cuts, list of candidates that stay)."""
    n = 0
    keep = []
    for i in cands:
        l = next(t for t in blocks if blocks[t][0] == i)
        ub = 2.0 * w * spinal_bound(l, N) + 1e-9
        xi = V(f"X_{i}")
        J = [t for t in range(3, l + 1) if profile[t - 2] >= 2]
        if not J:
            xi.UB = min(xi.UB, ub - W)
            if cutoff is None or ub > cutoff:
                keep.append(i)
            continue
        keep.append(i)
        M = xi.UB + W - ub
        if M <= 0:
            continue
        model.addConstr(xi + W <= ub + M * gp.quicksum(V(f"T_{blocks[t][0]}_{j}")
                                                       for t in J for j in blocks[t - 1][1:]),
                        name=f"SP_{l}")
        n += 1
    model.update()
    return n, keep


def build_profile_model(N, profile, model_name, opts, eta, cutoff, tail, log, rank_box=False, c2_layer=False):
    W = 0.5; EPS = W * float(eta); w = W - EPS
    model = cr.build_model(model_name, N, eta, log, tight=("TB" in opts))
    cr.apply_strengthening(model, model_name, N, eta, opts, log)
    model.update()
    V = model.getVarByName
    rel = (model_name == "Mrel")
    lay, blocks = layers_of(profile)
    H = len(profile)
    for i in range(1, N + 1):
        v = V(f"L_{i}"); v.LB = lay[i]; v.UB = lay[i]
    for i in range(1, N):
        v = V(f"LayerJump_{i}"); jb = 1.0 if lay[i + 1] > lay[i] else 0.0; v.LB = jb; v.UB = jb
    nfix = 0
    for i in range(2, N + 1):
        for j in range(1, i):
            if lay[i] != lay[j] + 1:
                V(f"T_{i}_{j}").UB = 0.0; V(f"F_{i}_{j}").UB = 0.0; nfix += 1
                for nm in ((f"Torque_{i}_{j}", f"P_{i}_{j}") if not rel else (f"D_{i}_{j}", f"MuBot_{i}_{j}", f"MuTop_{i}_{j}")):
                    v = V(nm)
                    if v is not None:
                        v.LB = 0.0; v.UB = 0.0
    model.update()
    n_rb = add_rank_box(model, N, profile, blocks, V, w, W) if (rank_box and not rel) else 0
    # objective: rightmost block of some layer l >= 2; with a cutoff, drop layers whose rightmost block
    # cannot beat it (X_first(l) <= X_1 + 2w(l-1), so overhang <= X_1^ub + W + 2w(l-1))
    removed = 0
    for g in model.getGenConstrs():
        if g.GenConstrName == "LinkMaxEdge":
            model.remove(g); removed += 1
    assert removed == 1, "LinkMaxEdge not found exactly once"
    mo = V("max_overhang")
    x1ub = V("X_1").UB
    e_ok, e_dE, e_ubs = True, None, {}
    if "EN" in opts:
        e_ok, e_dE, e_ubs = energy_profile_check(N, profile, eta)
        if "S1" not in opts:
            e_ubs = {}            # the energy layer bounds need sum_b X_b = -N eps
    cands = [blocks[l][0] for l in range(2, H + 1)
             if (cutoff is None or x1ub + W + 2 * w * (l - 1) > cutoff)
             and (cutoff is None or l not in e_ubs or e_ubs[l] > cutoff)]
    n_sp = 0
    if "SP" in opts and cands and e_ok:
        n_sp, cands = add_spinal_cuts(model, N, profile, blocks, V, w, W, cands, cutoff)
    if not cands or not e_ok:
        model.update()
        return model, {"pairs_fixed0": nfix, "noncross": 0, "c2_general": 0, "tail_m": 0,
                       "candidates": [], "rank_box": n_rb, "c2_layer": 0,
                       "energy_infeasible": not e_ok, "energy_dE": e_dE}
    model.addGenConstrMax(mo, [V(f"RE_{i}") for i in cands], name="LinkMaxEdgeProfile")
    # planar non-crossing rows between adjacent layers
    n_nc = 0
    thr = 3 if eta == 0 else 2
    for l in range(1, H):
        J, I = blocks[l], blocks[l + 1]
        for a in range(len(I)):
            for b_ in range(a + 1, len(I)):
                for c in range(len(J)):
                    for d in range(c + 1, len(J)):
                        i1, i2, j1, j2 = I[a], I[b_], J[c], J[d]
                        if (i2 - i1) + (j2 - j1) >= thr:
                            model.addConstr(V(f"T_{i1}_{j2}") + V(f"T_{i2}_{j1}") <= 1, name=f"NC_{i1}_{i2}_{j1}_{j2}")
                            n_nc += 1
    # singleton layers
    n_c2 = 0
    for l in range(1, H):
        if len(blocks[l]) != 1:
            continue
        u = blocks[l][0]; up = blocks[l + 1]; M = sum(profile[l:])
        for v in up:
            V(f"T_{v}_{u}").LB = 1.0
        if len(up) == 1:
            f = V(f"F_{up[0]}_{u}"); f.LB = float(M); f.UB = float(M)
        elif len(up) == 2 and not (l == 1 and "C2" in opts):
            v1, v2 = up
            f = V(f"F_{v1}_{u}"); FL, FU = 1.0, float(M - 1)
            if rel:
                Q = (V(f"MuTop_{v1}_{u}") - V(f"MuBot_{v1}_{u}")) + (V(f"MuTop_{v2}_{u}") - V(f"MuBot_{v2}_{u}")) \
                    - M * V(f"D_{v2}_{u}")
                d = V(f"D_{v1}_{u}") - V(f"D_{v2}_{u}"); dL, dU = 1.0, 4 * w
                model.addConstr(V(f"D_{v1}_{u}") >= 0, name=f"SgnD_{v1}_{u}")
                model.addConstr(V(f"D_{v2}_{u}") <= 0, name=f"SgnD_{v2}_{u}")
            else:
                tout = (V("Torque_1_0") if u == 1 else
                        gp.quicksum(V(f"Torque_{u}_{j}") for j in blocks[l - 1]))
                Q = tout - V(f"X_{u}") - M * V(f"P_{v2}_{u}")
                d = V(f"P_{v1}_{u}") - V(f"P_{v2}_{u}"); dL, dU = 1.0 - 2 * w, 2 * w
            model.addConstr(V(f"X_{v1}") >= V(f"X_{u}"), name=f"Sgn_{v1}_{u}")
            model.addConstr(V(f"X_{v2}") <= V(f"X_{u}"), name=f"Sgn_{v2}_{u}")
            model.addConstr(Q >= FL * d + dL * f - FL * dL, name=f"C2g_a_{u}")
            model.addConstr(Q >= FU * d + dU * f - FU * dU, name=f"C2g_b_{u}")
            model.addConstr(Q <= FU * d + dL * f - FU * dL, name=f"C2g_c_{u}")
            model.addConstr(Q <= FL * d + dU * f - FL * dU, name=f"C2g_d_{u}")
            n_c2 += 1
    # S3: harmonic single-column tail
    tail_m = 0
    if tail:
        k = 0
        while k < H and profile[H - 1 - k] == 1:
            k += 1
        if k >= 3:
            s = blocks[H - k + 1][0]; m = N - s; tail_m = m
            for r in range(s + 2, N + 1):
                model.addConstr(V(f"X_{r}") - V(f"X_{r - 1}") == w / (N - r + 1), name=f"S3_dx_{r}")
            if rel:   # mub_{s+1,s} = F (P - X_{s+1}) = m * (m-1) w / m
                model.addConstr(V(f"MuBot_{s + 1}_{s}") == (m - 1) * w, name="S3_base")
            else:
                model.addConstr(V(f"P_{s + 1}_{s}") == V(f"X_{s + 1}") + (m - 1) * w / m, name="S3_base")
    n_c2l = add_c2_layer_cuts(model, profile, blocks, V, rel, w) if c2_layer else 0
    if rank_box and not rel:      # singleton layers fixed T = 1 after the first pass
        tighten_contact_boxes(model, profile, blocks, V, w)
    model.update()
    return model, {"pairs_fixed0": nfix, "noncross": n_nc, "c2_general": n_c2, "tail_m": tail_m,
                   "candidates": cands, "rank_box": n_rb, "c2_layer": n_c2l, "spinal_cuts": n_sp}


def solve_one(job):
    (N, profile, model_name, opts, eta, cutoff, gap, seed, threads, tlim, tail, logdir, rank_box, c2_layer,
     params) = job
    t0 = time.time()
    lines = []
    model, info = build_profile_model(N, profile, model_name, set(opts), eta, cutoff, tail, lines.append,
                                      rank_box=rank_box, c2_layer=c2_layer)
    ptxt = "-".join(map(str, profile))
    if not info["candidates"]:      # pruned without solving
        model.dispose()
        if info.get("energy_infeasible"):   # energy test: no feasible stack has this profile
            return {"profile": ptxt, "H": len(profile), "status": "INFEASIBLE", "objective": float("nan"),
                    "bound": float("-inf"), "nodes": 0, "runtime_sec": time.time() - t0, "stack": None,
                    **info, "vars": 0, "constrs": 0}
        return {"profile": ptxt, "H": len(profile), "status": "CUTOFF", "objective": float("nan"),
                "bound": cutoff, "nodes": 0, "runtime_sec": time.time() - t0, "stack": None, **info,
                "vars": 0, "constrs": 0}
    model.setParam("OutputFlag", 1); model.setParam("LogToConsole", 0)
    if logdir:
        model.setParam("LogFile", os.path.join(logdir, f"P{ptxt}.log"))
    model.setParam("NonConvex", 2); model.setParam("Threads", threads)
    model.setParam("MIPGap", gap); model.setParam("Seed", seed)
    model.setParam("TimeLimit", float(tlim)); model.setParam("Cutoff", cutoff)
    model.setParam("NodefileStart", 1.0); model.setParam("SoftMemLimit", 4.0)
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
    nodes_extra = 0
    if st == GRB.INF_OR_UNBD:       # confirm infeasibility explicitly (presolve may report INF_OR_UNBD)
        nodes_extra = int(model.NodeCount)
        model.setParam("DualReductions", 0); model.reset(); model.optimize(); st = model.Status
    if st == GRB.CUTOFF:
        status = "CUTOFF"; bound = cutoff
    elif st == GRB.INFEASIBLE:
        status = "INFEASIBLE"; bound = float("-inf")
    elif st == GRB.OPTIMAL:
        status = "OPTIMAL"; bound = model.ObjBound
    else:
        status = base.STATUS_NAMES.get(st, f"STATUS_{st}")
        bound = max(model.ObjBound, cutoff) if model.SolCount == 0 else model.ObjBound
    obj = model.ObjVal if model.SolCount > 0 else float("nan")
    stack = None
    if model.SolCount > 0:
        val = {v.VarName: v.X for v in model.getVars()}
        stack = base.stack_from_values(N, val)
    res = {"profile": ptxt, "H": len(profile), "status": status, "objective": obj, "bound": bound,
           "nodes": int(model.NodeCount) + nodes_extra, "runtime_sec": time.time() - t0, "stack": stack, **info,
           "vars": model.NumVars, "constrs": model.NumConstrs}
    model.dispose()
    return res


def main():
    ap = argparse.ArgumentParser(description="Layer-profile decomposition (unconditional, n_2 = 2 via S2+)")
    ap.add_argument("N", type=int)
    ap.add_argument("--model", choices=["M1u", "Mrel"], default="M1u")
    ap.add_argument("--strengthen", default="TB,S1,S2P,C2")
    ap.add_argument("--tail", action="store_true", help="add S3 harmonic tail canonicalisation")
    ap.add_argument("--gap", type=float, default=5e-4)
    ap.add_argument("--lb", type=float, default=None, help="known feasible value (default certify_run.KNOWN_LB)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=1, help="threads per subproblem")
    ap.add_argument("--workers", type=int, default=1, help="subproblems solved in parallel")
    ap.add_argument("--time-limit", type=float, default=3600.0, help="seconds per subproblem")
    ap.add_argument("--outdir", default="results_profile")
    ap.add_argument("--keep-logs", action="store_true", help="write one Gurobi log per profile")
    ap.add_argument("--rank-box", action="store_true", help="layer-rank coordinate box (M1u)")
    ap.add_argument("--c2-layer", action="store_true", help="Family C2-Layer rows on two-block layers")
    ap.add_argument("--energy", action="store_true",
                    help="second-moment energy test: drop infeasible profiles and hopeless layers (paper: energy filter)")
    ap.add_argument("--spine", action="store_true",
                    help="SP: spinal cut X_s + W <= 2w S*_l(N) unless the right spine has a non-rightmost support")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="absolute cutoff (overrides lb*(1+gap)); e.g. spinal optimum + 1e-5")
    ap.add_argument("--param", action="append", default=[], help="Gurobi parameter NAME=VALUE (repeatable)")
    ap.add_argument("--ptag", default="", help="label for --param settings, appended to the tag")
    ap.add_argument("--nodefile-start", type=float, default=1.0, help="GB per subproblem before node files")
    ap.add_argument("--soft-mem-limit", type=float, default=4.0, help="GB per subproblem")
    ap.add_argument("--only", default=None,
                    help="comma-separated profiles to solve, e.g. 1-2-2-2-3-2-1 (default: all n_2 = 2 profiles)")
    ap.add_argument("--resume", action="store_true",
                    help="keep closed profiles (CUTOFF/INFEASIBLE/OPTIMAL) from an earlier run with the same tag")
    a = ap.parse_args()

    opts = {s.strip().upper() for s in a.strengthen.split(",") if s.strip()}
    if "S2P" not in opts:
        sys.exit("ERROR: the decomposition enumerates n_2 = 2 profiles only, so S2P is required")
    if a.model != "Mrel" and "C1" in opts:
        sys.exit("ERROR: C1 is for Mrel only")
    eta = 0.0
    lb = a.lb if a.lb is not None else cr.KNOWN_LB[a.N]
    cutoff = lb * (1.0 + a.gap) if a.cutoff is None else a.cutoff
    order = ["TB", "S1", "S2", "S2P", "S4", "C1", "C2"]
    if a.energy:
        opts.add("EN")
    if a.spine:
        opts.add("SP")
    stag = "".join("-" + o for o in order if o in opts) + ("-RB" if a.rank_box else "") + \
        ("-EN" if a.energy else "") + ("-SP" if a.spine else "") + \
        (f"-CUT{a.cutoff:.7f}" if a.cutoff is not None else "") + \
        ("-C2L" if a.c2_layer else "") + ("-S3" if a.tail else "") + (("-" + a.ptag) if a.ptag else "") + \
        ("-ONLY" if a.only else "")
    tag = f"{a.model}{stag}-PDcert_N{a.N}_s{a.seed}"
    os.makedirs(a.outdir, exist_ok=True)
    path = os.path.join(a.outdir, tag)
    logdir = None
    if a.keep_logs:
        logdir = path + "_logs"; os.makedirs(logdir, exist_ok=True)
    csv_path = f"{path}_profiles.csv"
    prev = {}
    if a.resume and os.path.exists(csv_path):
        with open(csv_path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                for k in ("objective", "bound", "runtime_sec"):
                    r[k] = float(r[k])
                for k in ("H", "nodes"):
                    r[k] = int(float(r[k]))
                r["stack"] = None
                prev[r["profile"]] = r          # the last row of a profile wins
        prev = {k: r for k, r in prev.items() if r["status"] in ("CUTOFF", "INFEASIBLE", "OPTIMAL")}
    log = base.TeeLog(f"{path}.out")

    profiles = enumerate_profiles(a.N, 2)
    if a.only:
        keep = {q.strip() for q in a.only.split(",") if q.strip()}
        profiles = [p for p in profiles if "-".join(map(str, p)) in keep]
        if len(profiles) != len(keep):
            sys.exit(f"ERROR: --only contains profiles that are not n_2 = 2 profiles of {a.N} blocks")
    # longest-first heuristic: many candidate contact pairs first
    profiles.sort(key=lambda p: -sum(p[l] * p[l + 1] for l in range(len(p) - 1)))
    ver = ".".join(map(str, gp.gurobi.version()))
    log("=" * 80)
    log(f"META: ProfileDecomp | model={a.model} strengthen={sorted(opts)} tail={a.tail} | N={a.N} | "
        f"profiles={len(profiles)} | LB={lb:.7f} | Cutoff={cutoff:.7f} | gap={a.gap} | seed={a.seed} | "
        f"workers={a.workers} x threads={a.threads} | TL/sub={a.time_limit:.0f}s")
    log(f"  Host={socket.gethostname()} | Gurobi={ver} | Start={datetime.now():%F %T} | PID={os.getpid()}")
    log("=" * 80)

    if prev:
        log(f"RESUME: {len(prev)} closed profiles kept from {csv_path}")
    nodedir = path + "_nodefiles"; os.makedirs(nodedir, exist_ok=True)
    solver_params = [f"NodefileStart={a.nodefile_start}", f"SoftMemLimit={a.soft_mem_limit}",
                     f"NodefileDir={nodedir}"] + list(a.param)
    jobs = [(a.N, p, a.model, sorted(opts), eta, cutoff, a.gap, a.seed, a.threads, a.time_limit, a.tail, logdir,
             a.rank_box, a.c2_layer, solver_params)
            for p in profiles if "-".join(map(str, p)) not in prev]
    t0 = time.time()
    rows = list(prev.values())
    fields = ["profile", "H", "status", "objective", "bound", "nodes", "runtime_sec",
              "pairs_fixed0", "noncross", "c2_general", "tail_m", "rank_box", "c2_layer", "energy_dE",
              "vars", "constrs"]
    new_file = not (a.resume and os.path.exists(csv_path))
    with open(csv_path, "w" if new_file else "a", newline="", encoding="utf-8") as fcsv:
        wr = csv.DictWriter(fcsv, fieldnames=fields, extrasaction="ignore")
        if new_file:
            wr.writeheader()
        if a.workers <= 1:
            it = map(solve_one, jobs)
            pool = None
        else:
            pool = mp.get_context("spawn").Pool(a.workers, maxtasksperchild=1)
            it = pool.imap_unordered(solve_one, jobs)
        for k, r in enumerate(it, 1):
            rows.append(r); wr.writerow(r); fcsv.flush()
            log(f"[{datetime.now():%F %T}] {k:3d}/{len(jobs)} P={r['profile']:<22s} {r['status']:10s} "
                f"bound={r['bound']:.7f} nodes={r['nodes']:>10d} t={r['runtime_sec']:9.1f}s")
        if pool:
            pool.close(); pool.join()
    wall = time.time() - t0

    best = max([r["objective"] for r in rows if r["objective"] == r["objective"]] + [lb])
    bound = max(r["bound"] for r in rows)
    closed = all(r["status"] in ("CUTOFF", "INFEASIBLE", "OPTIMAL") for r in rows)
    status = "OPTIMAL" if closed else "TIME_LIMIT"
    gap_rel = (bound - best) / best
    nodes = sum(r["nodes"] for r in rows)
    cpu = sum(r["runtime_sec"] for r in rows)
    slow = sorted(rows, key=lambda r: -r["runtime_sec"])[:5]
    best_row = max((r for r in rows if r["objective"] == r["objective"]), key=lambda r: r["objective"], default=None)
    base.atomic_write_json(f"{path}_result.json", {
        "method": f"ProfileDecomp {a.model} {sorted(opts)} rank_box={a.rank_box} c2_layer={a.c2_layer} "
                  f"tail={a.tail} params={a.param} (unconditional; n_2 = 2 by S2+)",
        "N": a.N, "eta": eta, "seed": a.seed, "threads": a.threads * a.workers, "mip_gap_target": a.gap,
        "time_limit_sec": a.time_limit, "gurobi_version": ver, "host": socket.gethostname(),
        "status": status, "objective": best, "objective_source": "subproblem" if best_row and best_row["objective"] >= lb else "known LB",
        "bound": bound, "gap_pct": 100.0 * gap_rel, "runtime_sec": wall, "cpu_sec": cpu, "nodes": nodes,
        "profiles": len(rows), "cutoff": cutoff, "milestones": {},
        "slowest": [{k: r[k] for k in ("profile", "status", "nodes", "runtime_sec")} for r in slow],
        "stack": best_row["stack"] if best_row else None,
    })
    log("")
    log("=" * 80)
    log("SUMMARY_REPORT:")
    log(f"  Status           : {status}")
    log(f"  ObjVal           : {best:.7f}")
    log(f"  BestBound        : {bound:.7f}")
    log(f"  Gap_rel          : {gap_rel:.3e}")
    log(f"  Runtime_Sec      : {wall:.1f}  (sum over subproblems {cpu:.1f})")
    log(f"  NodeCount        : {nodes}")
    log(f"  Slowest profiles : " + "; ".join(f"{r['profile']} {r['runtime_sec']:.1f}s/{r['nodes']}n" for r in slow))
    log("=" * 80)
    log.close()
    try:
        os.rmdir(nodedir)          # empty after every subproblem finished
    except OSError:
        pass


if __name__ == "__main__":
    main()
