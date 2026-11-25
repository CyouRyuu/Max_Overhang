#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
certify_m0.py
Long-horizon, multi-threaded certification runs of the UNCONDITIONAL baseline model M0
for the single-root maximum overhang problem.

Purpose
-------
The ablation study (run_ablation.py) compares M0..M4 under a fixed, single-thread, 3600 s budget.
This script is a separate "certification" experiment: it runs only M0 (no spinal / two-support /
empirical-envelope assumptions) with a tight gap, several threads and a multi-day time limit,
so that any OPTIMAL status is an unconditional global-optimality certificate.

Formulation
-----------
build_m0() is a verbatim copy of sections 1-5 of solve_overhang_tower_baseline() in
maxOverhang_baseline.py. Do not edit one without the other.

Extras over maxOverhang_baseline.py
-----------------------------------
* tight default gap (MIPGap = 1e-6) and records the wall time at which the gap first falls
  below 5e-4, 1e-4, 1e-5 and 1e-6 (milestones), so a run that stops early is still usable;
* every new incumbent is saved atomically to <tag>_incumbent.json (all variable values),
  which can be passed back with --start-json to warm-start a restarted run;
* bound/incumbent trace written to <tag>_trace.csv every --log-interval seconds;
* memory guards: node files go to disk after --nodefile-start GB, solve stops gracefully at
  --soft-mem-limit GB (Gurobi SoftMemLimit);
* SIGTERM / SIGINT terminate the solve gracefully and still write the final summary;
* final machine-readable result in <tag>_result.json plus the SUMMARY_REPORT block used by
  parse_results.py.
"""

import os
import sys
import csv
import json
import time
import math
import signal
import socket
import argparse
import platform
from datetime import datetime

import gurobipy as gp
from gurobipy import GRB


STOP_REQUESTED = False
MILESTONES = [5e-4, 1e-4, 1e-5, 1e-6]


def _request_stop(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True


class TeeLog:
    """Write messages to stdout and to a log file, flushing immediately."""

    def __init__(self, path):
        self.f = open(path, "a", buffering=1, encoding="utf-8")

    def __call__(self, msg=""):
        line = f"{msg}"
        print(line, flush=True)
        self.f.write(line + "\n")

    def close(self):
        self.f.close()


def build_m0(model, N, eta=0.0):
    """Verbatim copy of the M0 formulation (sections 1-5 of maxOverhang_baseline.py)."""
    W = 0.5
    MASS = 1.0
    EPSILON = W * float(eta)
    MAX_X = float(N)
    MAX_F = float(N) * MASS
    MAX_TORQUE = MAX_F * MAX_X
    BIG_M = MAX_X * 2.0 + 1.0

    # 1. Decision Variables
    X = {}
    X[0] = model.addVar(lb=0.0, ub=0.0, vtype=GRB.CONTINUOUS, name="X_0")
    X[1] = model.addVar(lb=-W + EPSILON, ub=W - 2.0 * EPSILON, vtype=GRB.CONTINUOUS, name="X_1")
    for i in range(2, N + 1):
        X[i] = model.addVar(lb=-MAX_X, ub=MAX_X - W, vtype=GRB.CONTINUOUS, name=f"X_{i}")

    L = {}
    L[0] = model.addVar(lb=0, ub=0, vtype=GRB.INTEGER, name="L_0")
    L[1] = model.addVar(lb=1, ub=1, vtype=GRB.INTEGER, name="L_1")
    for i in range(2, N + 1):
        L[i] = model.addVar(lb=2, ub=N, vtype=GRB.INTEGER, name=f"L_{i}")

    T = {}; F = {}; P = {}; Torque = {}
    for i in range(1, N + 1):
        for j in range(0, N + 1):
            if i == j:
                T[i, j] = model.addVar(lb=0.0, ub=0.0, vtype=GRB.BINARY, name=f"T_{i}_{j}")
                F[i, j] = model.addVar(lb=0.0, ub=0.0, vtype=GRB.CONTINUOUS, name=f"F_{i}_{j}")
                P[i, j] = model.addVar(lb=0.0, ub=0.0, vtype=GRB.CONTINUOUS, name=f"P_{i}_{j}")
                Torque[i, j] = model.addVar(lb=0.0, ub=0.0, vtype=GRB.CONTINUOUS, name=f"Torque_{i}_{j}")
            elif j == 0:
                t_val = 1.0 if i == 1 else 0.0
                T[i, 0] = model.addVar(lb=t_val, ub=t_val, vtype=GRB.BINARY, name=f"T_{i}_0")
                F[i, 0] = model.addVar(lb=0.0, ub=MAX_F * t_val, vtype=GRB.CONTINUOUS, name=f"F_{i}_0")
                P[i, 0] = model.addVar(lb=-MAX_X * t_val, ub=MAX_X * t_val, vtype=GRB.CONTINUOUS, name=f"P_{i}_0")
                Torque[i, 0] = model.addVar(lb=-MAX_TORQUE * t_val, ub=MAX_TORQUE * t_val,
                                            vtype=GRB.CONTINUOUS, name=f"Torque_{i}_0")
            else:
                T[i, j] = model.addVar(vtype=GRB.BINARY, name=f"T_{i}_{j}")
                F[i, j] = model.addVar(lb=0.0, ub=MAX_F, vtype=GRB.CONTINUOUS, name=f"F_{i}_{j}")
                P[i, j] = model.addVar(lb=-MAX_X, ub=MAX_X, vtype=GRB.CONTINUOUS, name=f"P_{i}_{j}")
                Torque[i, j] = model.addVar(lb=-MAX_TORQUE, ub=MAX_TORQUE,
                                            vtype=GRB.CONTINUOUS, name=f"Torque_{i}_{j}")

    # 2. Geometric Contact and DAG Topological Constraints
    for i in range(1, N + 1):
        model.addConstr(gp.quicksum(T[i, j] for j in range(0, N + 1) if i != j) >= 1, name=f"DAGSup_{i}")
        for j in range(0, N + 1):
            if i == j:
                continue
            model.addConstr(F[i, j] <= MAX_F * T[i, j])
            model.addConstr(P[i, j] >= -MAX_X * T[i, j])
            model.addConstr(P[i, j] <= MAX_X * T[i, j])
            model.addConstr(Torque[i, j] >= -MAX_TORQUE * T[i, j])
            model.addConstr(Torque[i, j] <= MAX_TORQUE * T[i, j])

            model.addConstr(P[i, j] >= X[i] - W + EPSILON - BIG_M * (1 - T[i, j]))
            model.addConstr(P[i, j] <= X[i] + W - EPSILON + BIG_M * (1 - T[i, j]))
            if j == 0:
                model.addConstr(P[i, 0] <= -EPSILON + BIG_M * (1 - T[i, 0]))
            else:
                model.addConstr(P[i, j] >= X[j] - W + EPSILON - BIG_M * (1 - T[i, j]))
                model.addConstr(P[i, j] <= X[j] + W - EPSILON + BIG_M * (1 - T[i, j]))

            model.addConstr(L[i] - L[j] <= 1 + N * (1 - T[i, j]))
            model.addConstr(L[i] - L[j] >= 1 - N * (1 - T[i, j]))
            model.addConstr(Torque[i, j] == F[i, j] * P[i, j], name=f"Bilinear_{i}_{j}")

    # 3. Same-Layer Non-Overlap and Layer Ordering Constraints
    for i in range(1, N):
        model.addConstr(L[i] <= L[i + 1], name=f"SymLayer_{i}")
        b = model.addVar(vtype=GRB.BINARY, name=f"LayerJump_{i}")
        model.addConstr(L[i + 1] - L[i] >= b)
        model.addConstr(L[i + 1] - L[i] <= N * b)
        model.addConstr(X[i] - X[i + 1] >= 2 * W - BIG_M * b, name=f"NonOverlap_{i}")

    # 4. Static Force and Moment Equilibrium Equations
    model.addConstr(gp.quicksum(X[i] for i in range(1, N + 1)) <= -float(N) * EPSILON, name="Global_CoM")
    for i in range(1, N + 1):
        f_in = gp.quicksum(F[k, i] for k in range(1, N + 1) if k != i) + MASS
        f_out = gp.quicksum(F[i, j] for j in range(0, N + 1) if j != i)
        model.addConstr(f_in == f_out, name=f"ForceBal_{i}")

        t_in = gp.quicksum(Torque[k, i] for k in range(1, N + 1) if k != i) + MASS * X[i]
        t_out = gp.quicksum(Torque[i, j] for j in range(0, N + 1) if j != i)
        model.addConstr(t_in == t_out, name=f"TorqueBal_{i}")

    # 5. Objective Function
    right_edges = [model.addVar(lb=-MAX_X, ub=MAX_X, vtype=GRB.CONTINUOUS, name=f"RE_{i}") for i in range(1, N + 1)]
    for i in range(1, N + 1):
        model.addConstr(right_edges[i - 1] == X[i] + W)

    max_overhang = model.addVar(lb=0.5 - 2.0 * EPSILON, ub=MAX_X, vtype=GRB.CONTINUOUS, name="max_overhang")
    model.addGenConstrMax(max_overhang, right_edges, name="LinkMaxEdge")
    model.setObjective(max_overhang, GRB.MAXIMIZE)
    model.update()

    return {"X": X, "L": L, "T": T, "F": F, "P": P, "max_overhang": max_overhang}


def atomic_write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)


def stack_from_values(N, val):
    """Build a human-readable stack description from a {var_name: value} dict."""
    rows = []
    for i in range(1, N + 1):
        layer = int(round(val[f"L_{i}"]))
        pos = val[f"X_{i}"]
        supports = [j for j in range(0, N + 1) if j != i and val.get(f"T_{i}_{j}", 0.0) > 0.5]
        rows.append({"block": i, "layer": layer, "center_x": pos, "supported_by": supports})
    rows.sort(key=lambda r: (-r["layer"], r["center_x"]))
    return rows


class CertifyCallback:
    def __init__(self, N, all_vars, tag_path, t0, log_interval, log):
        self.N = N
        self.all_vars = all_vars
        self.names = [v.VarName for v in all_vars]
        self.tag_path = tag_path
        self.t0 = t0
        self.log_interval = log_interval
        self.log = log
        self.last_log = -1e9
        self.milestones = {}
        self.terminated = False
        self.first_sol_time = None
        self.trace_path = f"{tag_path}_trace.csv"
        new = not os.path.exists(self.trace_path)
        self.trace_f = open(self.trace_path, "a", newline="", encoding="utf-8")
        self.trace = csv.writer(self.trace_f)
        if new:
            self.trace.writerow(["wall_sec", "incumbent", "bound", "gap_rel", "nodes"])

    def _gap(self, best, bnd):
        if best <= -1e99 or abs(best) < 1e-12:
            return math.inf
        return abs(bnd - best) / abs(best)

    def __call__(self, model, where):
        if STOP_REQUESTED and not self.terminated:
            self.terminated = True
            self.log(f"[{datetime.now():%F %T}] Stop signal received -> terminating solve gracefully.")
            model.terminate()
            return

        if where == GRB.Callback.MIP:
            now = time.time() - self.t0
            best = model.cbGet(GRB.Callback.MIP_OBJBST)
            bnd = model.cbGet(GRB.Callback.MIP_OBJBND)
            nodes = model.cbGet(GRB.Callback.MIP_NODCNT)
            gap = self._gap(best, bnd)
            for m in MILESTONES:
                if m not in self.milestones and gap <= m:
                    self.milestones[m] = {"wall_sec": round(now, 2), "incumbent": best,
                                          "bound": bnd, "nodes": int(nodes)}
                    self.log(f"[{datetime.now():%F %T}] MILESTONE gap<={m:g} at {now:.1f}s "
                             f"(inc={best:.7f}, bnd={bnd:.7f}, nodes={int(nodes)})")
            if now - self.last_log >= self.log_interval:
                self.last_log = now
                self.trace.writerow([f"{now:.1f}", f"{best:.9f}", f"{bnd:.9f}",
                                     f"{gap:.3e}" if math.isfinite(gap) else "inf", int(nodes)])
                self.trace_f.flush()

        elif where == GRB.Callback.MIPSOL:
            now = time.time() - self.t0
            if self.first_sol_time is None:
                self.first_sol_time = now
            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
            vals = model.cbGetSolution(self.all_vars)
            val = dict(zip(self.names, vals))
            atomic_write_json(f"{self.tag_path}_incumbent.json", {
                "N": self.N, "objective": obj, "found_at_sec": round(now, 2),
                "stack": stack_from_values(self.N, val), "values": val,
            })
            self.log(f"[{datetime.now():%F %T}] New incumbent {obj:.9f} at {now:.1f}s (saved)")

    def close(self):
        self.trace_f.close()


STATUS_NAMES = {
    GRB.OPTIMAL: "OPTIMAL", GRB.TIME_LIMIT: "TIMEOUT", GRB.INTERRUPTED: "INTERRUPTED",
    GRB.INFEASIBLE: "INFEASIBLE", GRB.MEM_LIMIT: "MEM_LIMIT", GRB.NODE_LIMIT: "NODE_LIMIT",
}


def main():
    ap = argparse.ArgumentParser(description="Long-horizon unconditional certification of M0")
    ap.add_argument("N", type=int)
    ap.add_argument("--time-limit", type=float, default=259200, help="seconds (default 72 h)")
    ap.add_argument("--threads", "-t", type=int, default=4)
    ap.add_argument("--gap", type=float, default=1e-6, help="relative MIPGap (default 1e-6)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--eta", type=float, default=0.0)
    ap.add_argument("--soft-mem-limit", type=float, default=6.0, help="GB, graceful stop (default 6)")
    ap.add_argument("--nodefile-start", type=float, default=2.0, help="GB before nodes go to disk (default 2)")
    ap.add_argument("--log-interval", type=float, default=60.0, help="trace interval in seconds")
    ap.add_argument("--start-json", type=str, default=None, help="incumbent JSON from a previous run (MIP start)")
    ap.add_argument("--outdir", type=str, default="results_certify_M0")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    eta_sfx = f"_eta{args.eta:.2f}" if args.eta > 0 else ""
    tag = f"M0cert_N{args.N}{eta_sfx}_s{args.seed}"
    tag_path = os.path.join(args.outdir, tag)
    nodedir = os.path.join(args.outdir, f"nodefiles_{tag}")
    os.makedirs(nodedir, exist_ok=True)

    log = TeeLog(f"{tag_path}.out")
    signal.signal(signal.SIGTERM, _request_stop)
    signal.signal(signal.SIGINT, _request_stop)

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 1)
    env.setParam("LogToConsole", 0)
    env.setParam("LogFile", f"{tag_path}_gurobi.log")
    env.start()
    model = gp.Model(f"M0_certify_N{args.N}", env=env)
    v = build_m0(model, args.N, args.eta)

    model.setParam("NonConvex", 2)
    model.setParam("Threads", args.threads)
    model.setParam("MIPGap", args.gap)
    model.setParam("Seed", args.seed)
    model.setParam("TimeLimit", float(args.time_limit))
    model.setParam("NodefileStart", args.nodefile_start)
    model.setParam("NodefileDir", nodedir)
    model.setParam("SoftMemLimit", args.soft_mem_limit)

    if args.start_json:
        with open(args.start_json, encoding="utf-8") as f:
            start_vals = json.load(f)["values"]
        n_set = 0
        for var in model.getVars():
            if var.VarName in start_vals:
                var.Start = start_vals[var.VarName]
                n_set += 1
        log(f"INFO: MIP start loaded from {args.start_json} ({n_set} values)")

    ver = ".".join(map(str, gp.gurobi.version()))
    log("=" * 80)
    log(f"META: Method=Baseline_Certify | N={args.N} | Eta={args.eta:.4f} | Seed={args.seed} | "
        f"TimeLimit={args.time_limit:.0f}s | Threads={args.threads} | MIPGap={args.gap}")
    log(f"  Host={socket.gethostname()} | Gurobi={ver} | Python={platform.python_version()} | "
        f"Start={datetime.now():%F %T} | PID={os.getpid()}")
    log(f"  SoftMemLimit={args.soft_mem_limit} GB | NodefileStart={args.nodefile_start} GB | "
        f"Vars={model.NumVars} | LinConstrs={model.NumConstrs} | QConstrs={model.NumQConstrs}")
    log("=" * 80)

    t0 = time.time()
    cb = CertifyCallback(args.N, model.getVars(), tag_path, t0, args.log_interval, log)
    try:
        model.optimize(cb)
    except gp.GurobiError as e:
        log(f"[{datetime.now():%F %T}] GUROBI ERROR: {e}")
        cb.close()
        log.close()
        _cleanup_nodedir(nodedir)
        sys.exit(2)
    finally:
        if not cb.trace_f.closed:
            cb.close()
    runtime = time.time() - t0

    status = STATUS_NAMES.get(model.Status, f"STATUS_{model.Status}")
    has_sol = model.SolCount > 0
    obj = model.ObjVal if has_sol else float("nan")
    bnd = model.ObjBound
    gap_pct = model.MIPGap * 100.0 if has_sol else float("nan")
    nodes = int(model.NodeCount)

    stack = None
    if has_sol:
        val = {var.VarName: var.X for var in model.getVars()}
        stack = stack_from_values(args.N, val)
        model.write(f"{tag_path}.sol")

    atomic_write_json(f"{tag_path}_result.json", {
        "method": "Baseline_Certify (M0, unconditional)", "N": args.N, "eta": args.eta,
        "seed": args.seed, "threads": args.threads, "mip_gap_target": args.gap,
        "time_limit_sec": args.time_limit, "gurobi_version": ver, "host": socket.gethostname(),
        "status": status, "objective": obj, "bound": bnd, "gap_pct": gap_pct,
        "runtime_sec": runtime, "nodes": nodes, "sol_count": model.SolCount,
        "first_sol_sec": cb.first_sol_time,
        "milestones": {f"{k:g}": val for k, val in sorted(cb.milestones.items(), reverse=True)},
        "tolerances": {"FeasibilityTol": model.Params.FeasibilityTol,
                       "IntFeasTol": model.Params.IntFeasTol, "OptimalityTol": model.Params.OptimalityTol},
        "stack": stack,
    })

    log("")
    log("=" * 80)
    log("SUMMARY_REPORT:")
    log(f"  Method           : Baseline_Certify")
    log(f"  N                : {args.N}")
    log(f"  Eta              : {args.eta:.4f}")
    log(f"  Status           : {status}")
    log(f"  ObjVal           : {obj:.7f}")
    log(f"  BestBound        : {bnd:.7f}")
    log(f"  MIPGap_Pct       : {gap_pct:.6f}%")
    log(f"  Runtime_Sec      : {runtime:.4f}")
    log(f"  FirstSol_Sec     : {cb.first_sol_time if cb.first_sol_time is not None else float('nan'):.4f}")
    log(f"  NodeCount        : {nodes}")
    log(f"  SolCount         : {model.SolCount}")
    for k, m in sorted(cb.milestones.items(), reverse=True):
        log(f"  Milestone gap<={k:<7g}: {m['wall_sec']:.1f}s")
    log("=" * 80)
    if stack:
        log("\nSTACK_CONFIGURATION:")
        for r in stack:
            log(f"  Layer {r['layer']:2d} | Block {r['block']:2d} | Center X: {r['center_x']:10.7f} | "
                f"Supported by: {r['supported_by']}")
    log.close()
    _cleanup_nodedir(nodedir)


def _cleanup_nodedir(nodedir):
    """Node files are no longer needed once the solve has ended."""
    try:
        for fn in os.listdir(nodedir):
            os.remove(os.path.join(nodedir, fn))
        os.rmdir(nodedir)
    except OSError:
        pass


if __name__ == "__main__":
    main()
