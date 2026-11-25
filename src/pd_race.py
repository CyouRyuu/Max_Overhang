#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pd_race.py -- solve the profiles that are still open in a profile_run.py CSV, in a SEPARATE process
(e.g. with more threads / another seed), without touching the running profile_run.py.

It reads <csv>, enumerates all n_2 = 2 profiles for N, and solves every profile that has no closed row
(CUTOFF / INFEASIBLE / OPTIMAL) with exactly the given model options.  Results go to <outdir> (one JSON and
one Gurobi log per profile).  Whichever run closes a profile first, the certificate is the same:
the profile cannot beat the cutoff.

Example (N = 13, layer-profile decomposition without SP, Table of methodological progression):
  nohup nice -n 10 python3 pd_race.py 13 results_pd13/<tag>_profiles.csv \
      --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --param Heuristics=0.0 \
      --threads 3 --seed 1 --outdir results_pd13_race > pd_race.log 2>&1 &
"""
import os
import sys
import csv
import json
import argparse

import certify_run as cr
import profile_run as pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("csv")
    ap.add_argument("--model", default="M1u")
    ap.add_argument("--strengthen", default="TB,S1,S2P,S4,C2")
    ap.add_argument("--rank-box", action="store_true")
    ap.add_argument("--c2-layer", action="store_true")
    ap.add_argument("--tail", action="store_true")
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--gap", type=float, default=5e-4)
    ap.add_argument("--cutoff", type=float, default=None, help="absolute cutoff (overrides LB*(1+gap))")
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--time-limit", type=float, default=604800)
    ap.add_argument("--nodefile-start", type=float, default=1.5)
    ap.add_argument("--soft-mem-limit", type=float, default=8.0)
    ap.add_argument("--outdir", default="results_race")
    a = ap.parse_args()

    closed = set()
    with open(a.csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["status"] in ("CUTOFF", "INFEASIBLE", "OPTIMAL"):
                closed.add(r["profile"])
    todo = [p for p in pr.enumerate_profiles(a.N, 2) if "-".join(map(str, p)) not in closed]
    print(f"{len(closed)} closed profiles in {a.csv}; open: {['-'.join(map(str, p)) for p in todo]}", flush=True)
    if not todo:
        return
    os.makedirs(a.outdir, exist_ok=True)
    nodedir = os.path.join(a.outdir, "nodefiles"); os.makedirs(nodedir, exist_ok=True)
    opts = sorted({s.strip().upper() for s in a.strengthen.split(",") if s.strip()})
    cutoff = cr.KNOWN_LB[a.N] * (1.0 + a.gap) if a.cutoff is None else a.cutoff
    params = [f"NodefileStart={a.nodefile_start}", f"SoftMemLimit={a.soft_mem_limit}",
              f"NodefileDir={nodedir}"] + list(a.param)
    for p in todo:
        ptxt = "-".join(map(str, p))
        print(f"solving {ptxt} with {a.threads} threads, seed {a.seed}, cutoff {cutoff:.8f} ...", flush=True)
        job = (a.N, p, a.model, opts, 0.0, cutoff, a.gap, a.seed, a.threads, a.time_limit, a.tail, a.outdir,
               a.rank_box, a.c2_layer, params)
        r = pr.solve_one(job)
        r.pop("stack", None)
        with open(os.path.join(a.outdir, f"P{ptxt}_s{a.seed}_result.json"), "w", encoding="utf-8") as f:
            json.dump({**r, "seed": a.seed, "threads": a.threads, "cutoff": cutoff, "options": vars(a)}, f, indent=2)
        print(f"  {ptxt}: {r['status']}  bound={r['bound']:.7f}  nodes={r['nodes']}  t={r['runtime_sec']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
