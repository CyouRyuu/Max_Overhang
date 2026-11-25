#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_certify.py
Collect *_result.json files written by certify_run.py / certify_m0.py into one CSV,
and print a compact table. For runs still in progress (no result.json yet), report the
last row of the *_trace.csv file instead (status RUNNING).

Usage:  python3 summarize_certify.py <dir> [<dir> ...] [--csv out.csv]
"""

import os
import sys
import csv
import glob
import json
import argparse
import re


def last_trace_row(path):
    last = None
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            last = row
    return last


def collect(dirs):
    rows = []
    for d in dirs:
        done = set()
        for p in sorted(glob.glob(os.path.join(d, "*_result.json"))):
            with open(p, encoding="utf-8") as f:
                r = json.load(f)
            tag = os.path.basename(p)[: -len("_result.json")]
            done.add(tag)
            ms = r.get("milestones") or {}
            rows.append({
                "dir": d, "tag": tag, "model": tag.split("cert_")[0], "N": r["N"], "seed": r.get("seed"), "threads": r.get("threads"),
                "status": r["status"], "objective": r["objective"], "bound": r["bound"],
                "gap_rel": (r["gap_pct"] / 100.0) if r["gap_pct"] == r["gap_pct"] else float("nan"),
                "runtime_h": r["runtime_sec"] / 3600.0, "nodes": r["nodes"],
                "t_gap_5e-4_h": ms.get("0.0005", {}).get("wall_sec", float("nan")) / 3600.0 if ms.get("0.0005") else float("nan"),
                "t_gap_1e-4_h": ms.get("0.0001", {}).get("wall_sec", float("nan")) / 3600.0 if ms.get("0.0001") else float("nan"),
                "t_gap_1e-5_h": ms.get("1e-05", {}).get("wall_sec", float("nan")) / 3600.0 if ms.get("1e-05") else float("nan"),
            })
        for p in sorted(glob.glob(os.path.join(d, "*_trace.csv"))):
            tag = os.path.basename(p)[: -len("_trace.csv")]
            if tag in done:
                continue
            r = last_trace_row(p)
            if not r:
                continue
            N = int(tag.split("_N")[1].split("_")[0])
            inc, bnd = float(r["incumbent"]), float(r["bound"])
            m_seed = re.search(r"_s(\d+)$", tag)
            rows.append({
                "dir": d, "tag": tag, "model": tag.split("cert_")[0], "N": N,
                "seed": int(m_seed.group(1)) if m_seed else None, "threads": "",
                "status": "RUNNING", "objective": inc, "bound": bnd,
                "gap_rel": float(r["gap_rel"]) if r["gap_rel"] != "inf" else float("inf"),
                "runtime_h": float(r["wall_sec"]) / 3600.0, "nodes": int(r["nodes"]),
                "t_gap_5e-4_h": float("nan"), "t_gap_1e-4_h": float("nan"), "t_gap_1e-5_h": float("nan"),
            })
    rows.sort(key=lambda x: (x["N"], x["model"], x["seed"] or 0))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    rows = collect(a.dirs)
    if not rows:
        print("no results found")
        return
    hdr = f"{'model':5s} {'N':>3s} {'seed':>4s} {'status':11s} {'objective':>11s} {'bound':>11s} {'gap':>9s} {'time_h':>8s} {'nodes':>14s} {'t5e-4_h':>8s} {'t1e-4_h':>8s}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['model']:5s} {r['N']:3d} {str(r['seed']):>4s} {r['status']:11s} {r['objective']:11.7f} {r['bound']:11.7f} "
              f"{r['gap_rel']:9.2e} {r['runtime_h']:8.3f} {r['nodes']:14,d} "
              f"{r['t_gap_5e-4_h']:8.3f} {r['t_gap_1e-4_h']:8.3f}")
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nCSV written to {a.csv}")


if __name__ == "__main__":
    main()
