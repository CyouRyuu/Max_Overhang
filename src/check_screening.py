"""Exact re-check of the analytic screening steps at eta = 0, without an optimization solver.

The runs decide in double precision (profiles.energy_profile_check, pd_certificate.closure_reason) which profiles
need no model and which candidate layers are dropped from a profile subproblem.  This script repeats every such
decision in exact rational arithmetic, with the cutoff equal to the exact value of the double passed to the solver,
and writes one record per decision:

  certificates (cutoff gamma_N, N = 5..16):
    * energy test: Delta E(n) < 0, decided exactly (Delta E is a multiple of 1/4 at eta = 0);
    * every layer l >= 2 of every remaining profile that the runs drop from the candidate set, with the criterion
      that drops it, each checked exactly:
        base     x1ub + W + 2w(l - 1) <= gamma                               (rational)
        energy   z^E(n, l) <= gamma, the energy layer bound, decided without square roots:
                 case 1 (m = 1 or 4 A_m dE <= m^2):  c = gamma - w - (m - 1)/2 >= 0 and c^2 >= dE / A_m,
                 case 2:                             c = gamma - w >= 0 and c^2 >= (1/A_m + (m-1)/m)(dE + m(m-1)/4)
        sp       2w S^_l(N) <= gamma for a pure chain (J_l empty), with the rigorous interval enclosure of S^_l(N)
                 (spine_faces.shat_enclosure, outward rounding);
  non-spinal branches (threshold gamma_ns(N) = 0.99 S*(N) as passed to the solver, N = 5..16):
    * every branch (n, l*) with J_{l*} nonempty is tested for a conflict between the lower bounds
      X_{s_l} >= gamma_ns - W - 2w(l* - l) of the branch and the upper bounds x1ub + 2w(l - 1) (layer-rank box) and,
      if J_l is empty, 2w S^_l(N) - W (pure chain); a conflict proves the branch empty without a solver.  The
      archived records with status INFEASIBLE and no node are then split into branches closed by such an exact
      conflict and branches closed by the solver in presolve.

Outputs (default results/screening/): screening_N<N>.csv (one row per decision) and screening_summary.txt.
Every double-precision decision of the runs must agree with the exact one; any disagreement fails the check.

  python3 src/check_screening.py [--out results/screening]
"""
import argparse
import csv
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profiles as pf  # noqa: E402
from pd_certificate import closure_reason, split_candidate_layers  # noqa: E402
from spine_faces import shat_enclosure  # noqa: E402
from mpmath import iv, mp  # noqa: E402

mp.dps = 50
iv.dps = 50
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W = Fr(1, 2)
w = Fr(1, 2)                                     # eta = 0


def exact_dE(p):
    return w * w * sum((l - 1) * n for l, n in enumerate(p, start=1)) - sum(Fr(n * (n * n - 1), 12) for n in p)


def energy_case(N, m, dE):
    Am = m + Fr(m * m) / ((N - m - 1) + Fr(1, N + 1))
    return Am, (m == 1 or 4 * Am * dE <= m * m)


def energy_le(N, m, dE, c):
    """Exact test z^E <= c of the energy layer bound; returns (holds, interval lower bound of c - z^E)."""
    Am, case1 = energy_case(N, m, dE)
    if case1:
        d = c - w - Fr(m - 1, 2)
        rad = dE / Am
    else:
        d = c - w
        rad = (1 / Am + Fr(m - 1, m)) * (dE + Fr(m * (m - 1), 4))
    holds = d >= 0 and d * d >= rad
    marg = iv.mpf(d.numerator) / d.denominator - iv.sqrt(iv.mpf(rad.numerator) / rad.denominator)
    return holds, float(marg.a)


def x1ub(N):
    return Fr(N - 3, 2 * N)                     # (N - 3) w / N - eps at eta = 0


def shat_upper(l, N):
    m, e = mp.mpf(shat_enclosure(l, N).b).man_exp   # exact rational value of the upper end of the enclosure
    return Fr(m) * Fr(2) ** e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "results", "screening"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    summary, bad = [], []
    shat_cache = {}

    def shat_hi(l, N):
        if (l, N) not in shat_cache:
            shat_cache[(l, N)] = shat_upper(l, N)
        return shat_cache[(l, N)]

    for N in range(5, 17):
        cut_d = float(pf.CERT_UB[N])
        cut = Fr(cut_d)                          # exact value of the double passed to the solver
        rows = []
        n_energy = n_nocand = n_layers = 0
        for p in pf.enumerate_profiles(N, 2):
            pt = "-".join(map(str, p))
            dE = exact_dE(p)
            ok_d, _, _ = pf.energy_profile_check(N, p)
            if ok_d != (dE >= 0):
                bad.append(("energy decision", N, pt))
            if dE < 0:
                n_energy += 1
                rows.append(["cert", pt, "", "energy", float(dE), cut_d])
                continue
            kept_d = split_candidate_layers(N, p, cut_d)
            kept_d = [l for l in kept_d if [t for t in range(3, l + 1) if p[t - 2] >= 2]
                      or 2.0 * float(w) * pf.spinal_bound(l, N) + 1e-9 > cut_d]
            if closure_reason(N, p, cut_d, True) is not None:
                n_nocand += 1
            for l in range(2, len(p) + 1):
                if l in kept_d:
                    continue
                n_layers += 1
                m = p[l - 1]
                base = x1ub(N) + W + 2 * w * (l - 1)
                if base <= cut:
                    rows.append(["cert", pt, l, "base", float(cut - base), cut_d]); continue
                holds, marg = energy_le(N, m, dE, cut)
                if holds:
                    rows.append(["cert", pt, l, "energy layer", marg, cut_d]); continue
                J = [t for t in range(3, l + 1) if p[t - 2] >= 2]
                if not J and 2 * w * shat_hi(l, N) <= cut:
                    rows.append(["cert", pt, l, "sp pure chain", float(cut - 2 * w * shat_hi(l, N)), cut_d]); continue
                bad.append(("dropped candidate layer not confirmed", N, pt, l))
        # non-spinal branches
        g_d = pf.spinal_opt(N) * (1.0 - 1e-2)   # threshold passed to the solver by ns_check.py (2w = 1)
        g = Fr(g_d)
        n_br = n_empty = n_conf = 0
        conflict = set()
        for p in pf.enumerate_profiles(N, 2):
            if exact_dE(p) < 0:
                continue
            pt = "-".join(map(str, p))
            for ls in range(2, len(p) + 1):
                n_br += 1
                if not [t for t in range(3, ls + 1) if p[t - 2] >= 2]:
                    n_empty += 1
                    continue
                found = None
                for l in range(1, ls):
                    lb = g - W - 2 * w * (ls - l)
                    ub = x1ub(N) + 2 * w * (l - 1)
                    if lb > ub:
                        found = (l, "rank box", float(lb - ub)); break
                    if l >= 2 and not [t for t in range(3, l + 1) if p[t - 2] >= 2]:
                        ub2 = 2 * w * shat_hi(l, N) - W
                        if lb > ub2:
                            found = (l, "sp pure chain", float(lb - ub2)); break
                if found:
                    n_conf += 1
                    conflict.add((pt, ls))
                    rows.append(["ns", pt, ls, f"bounds conflict at layer {found[0]} ({found[1]})", found[2], g_d])
        # archived non-spinal records closed without nodes
        arch = os.path.join(HERE, "results", "nonspinal", f"ns_N{N}_rel0.01.csv")
        raw = os.path.join(HERE, "results", "raw_logs", f"ns_N{N}", f"ns_N{N}_rel0.01.csv")
        zero_inf = zero_exact = 0
        if os.path.exists(raw):                  # rerun with closure kinds: check the bound closures directly
            for r in csv.DictReader(open(raw)):
                if r.get("closure") == "bounds before solving" and (r["profile"], int(r["lstar"])) not in conflict:
                    bad.append(("bounds-before-solving record without exact conflict", N, r["profile"], r["lstar"]))
        if os.path.exists(arch):
            seen = {}
            for r in csv.DictReader(open(arch)):
                seen[(r["profile"], int(r["lstar"]))] = r
            for k, r in seen.items():
                if r["status"] == "INFEASIBLE" and float(r.get("nodes") or 0) == 0:
                    zero_inf += 1
                    zero_exact += k in conflict
                    if (r.get("closure") or "") == "bounds before solving" and k not in conflict:
                        bad.append(("bounds-before-solving record without exact conflict", N, k))
        with open(os.path.join(a.out, f"screening_N{N}.csv"), "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["kind", "profile", "layer_or_lstar", "criterion", "margin_lower_bound", "threshold"])
            wr.writerows(rows)
        summary.append(f"N={N:2d} gamma_N={cut_d:.7f}: energy-pruned profiles {n_energy}, profiles without candidate "
                       f"layer {n_nocand}, dropped candidate layers {n_layers} (all confirmed exactly); "
                       f"gamma_ns={g_d:.10f}: branches {n_br}, J empty {n_empty}, exact bound conflicts {n_conf}; "
                       f"archived INFEASIBLE records without nodes {zero_inf}, of which {zero_exact} exact bound "
                       f"conflicts and {zero_inf - zero_exact} closed by the solver in presolve")
    summary.append("disagreements: " + (str(bad) if bad else "none"))
    txt = "\n".join(summary)
    open(os.path.join(a.out, "screening_summary.txt"), "w").write(txt + "\n")
    print(txt)
    if bad:
        sys.exit(1)
    print("ALL SCREENING DECISIONS CONFIRMED EXACTLY")


if __name__ == "__main__":
    main()
