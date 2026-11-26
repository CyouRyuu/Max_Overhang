"""Write results/INDEX.csv: one row per archived run, with its solver version, model options and parameters,
cutoff or threshold, threads, time limit and files.  Needs no solver.

  python3 src/make_index.py            (run from the repository root)

Profile runs (profile_run.py, certify_m0.py) record these fields in their *_result.json.  The non-spinal runs
(ns_check.py) record only branches; their settings are those of scripts/nonspinal*.sh, tolerance_check_ns16.sh
and reverse_controls.sh (single-thread branches, Heuristics = 0, IntFeasTol = 1e-9, MIPGap = 1e-4) and their
threshold is read from ns_summary.txt.
"""
import csv
import glob
import json
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RES = os.path.join(ROOT, "results")
NS_SETTINGS = ("ns_check.py: threads 1 per branch, Heuristics=0, IntFeasTol=1e-9, MIPGap=1e-4, NonConvex=2; "
               "time limit 7200 s per branch (86400 s for the N = 16 t0 split); threshold shown to 7 decimals")


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, "**", "*_result.json"), recursive=True)):
        d = json.load(open(f))
        if "method" not in d:
            continue
        rel = os.path.relpath(f, ROOT)
        rows.append({"record": rel, "kind": rel.split(os.sep)[1], "N": d.get("N"),
                     "gurobi": d.get("gurobi_version"), "settings": d.get("method"),
                     "cutoff_or_threshold": d.get("cutoff"), "threads": d.get("threads"),
                     "time_limit_sec": d.get("time_limit_sec"), "status": d.get("status"),
                     "bound": d.get("bound"), "nodes": d.get("nodes"),
                     "files": ";".join(sorted(os.path.basename(g) for g in glob.glob(f.replace("_result.json", "*"))))})
    for f in sorted(glob.glob(os.path.join(RES, "**", "ns_summary.txt"), recursive=True)):
        rel = os.path.relpath(os.path.dirname(f), ROOT)
        for line in open(f):
            m = re.match(r"N=(\d+) gamma=([0-9.]+).*branches=(\d+) open=(\d+) nodes=(\d+)", line)
            split = re.match(r"N=(\d+) .*t0.*?\) gamma=([0-9.]+).*branches=(\d+) open=(\d+) nodes=(\d+)", line)
            tol = "tolerance_check" in rel
            settings = NS_SETTINGS
            if split:   # t0 sub-branches (Appendix B): the split of the main run, or its tolerance check
                m = split
                settings = NS_SETTINGS.replace("time limit 7200 s per branch (86400 s for the N = 16 t0 split)",
                                               "t0 sub-branches only, time limit "
                                               + ("172800 s; FeasibilityTol=1e-9, OptimalityTol=1e-9, "
                                                  "NumericFocus=3" if tol else "86400 s"))
                settings = settings.replace("; threshold shown to 7 decimals", "") if tol else settings
            elif tol:
                m = None
            if m:
                rows.append({"record": rel, "kind": rel.split(os.sep)[1], "N": int(m.group(1)), "gurobi": "13.0.1",
                             "settings": settings, "cutoff_or_threshold": m.group(2), "threads": 1,
                             "time_limit_sec": "", "status": f"branches {m.group(3)}, open {m.group(4)}",
                             "bound": "", "nodes": m.group(5),
                             "files": ";".join(sorted(os.path.basename(g) for g in glob.glob(
                                 os.path.join(os.path.dirname(f), f"ns_N{m.group(1)}_*"))
                                 if ("only_t0" in g) == bool(split)))})
    keys = ["record", "kind", "N", "gurobi", "settings", "cutoff_or_threshold", "threads", "time_limit_sec",
            "status", "bound", "nodes", "files"]
    with open(os.path.join(RES, "INDEX.csv"), "w", newline="", encoding="utf-8") as fo:
        w = csv.DictWriter(fo, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"results/INDEX.csv: {len(rows)} runs")


if __name__ == "__main__":
    main()
