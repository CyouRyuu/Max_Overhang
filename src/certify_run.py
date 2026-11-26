#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
certify_run.py
Monolithic certification runs and the shared model builders used by the layer-profile decomposition.

Models used in the paper
------------------------
M0   : baseline formulation M_0 (paper, Section 2), built by certify_m0.build_m0.
M1u  : the order-pruned formulation of ablation_m1_order_prune.py with the support-degree bound set to
       the geometric value sum_j T_ij <= 3 (<= 2 when eta > 0), which holds for every feasible stack.
       The layer-profile decomposition (profile_run.py) builds every subproblem from M1u.
Mrel : relative-coordinate formulation (optional, not used for the certificates in the paper).

The variants M2u, M3u, M1w, M2w, M3w, M1, M2, M3 of an earlier development stage need model files that are
not part of this repository; selecting them raises an error.

The M1u model file is used unchanged: the script intercepts its execute_and_report() call to obtain the
built model, then edits the right-hand side of the DAGSupMax_* constraints.

The --strengthen option adds the structural reductions of Sections 3 and 4 (TB, S1, S2P, S4, C2, ...);
the dictionaries KNOWN_LB and CERT_UB below supply the bounds needed by the base-singleton check (S2P).
"""

import os
import sys
import json
import time
import signal
import socket
import argparse
import platform
from datetime import datetime

import gurobipy as gp
from gurobipy import GRB

import certify_m0 as base

MODELS = ["M0", "M1u", "M2u", "M3u", "M1w", "M2w", "M3w", "M1", "M2", "M3", "Mrel"]
# M1/M2/M3 = earlier models with the two-support bound (<= 2) kept: CONDITIONAL
# (two-support subclass only). Offered only for controlled comparisons with M1u/M2u/M3u.
CONDITIONAL_MODELS = {"M1", "M2", "M3"}
GEOMETRIC_MAX_SUPPORTS = 3


def relax_two_support(model, N, log):
    model.update()
    count = 0
    for c in model.getConstrs():
        if c.ConstrName.startswith("DAGSupMax_"):
            c.RHS = float(GEOMETRIC_MAX_SUPPORTS)
            count += 1
    model.update()
    if count != N - 1:
        raise RuntimeError(f"Expected {N - 1} DAGSupMax constraints, found {count}; formulation changed?")
    log(f"INFO: relaxed {count} two-support rows (sum_j T_ij <= 2) to the geometric bound <= {GEOMETRIC_MAX_SUPPORTS}")


def add_window_cuts(model, N, log):
    """Valid unconditional cuts from same-layer spacing (centres >= 2W = 1 apart) and contact
    reach (|X_i - X_j| <= 1). Layers are monotone in the index, so blocks with indices between two
    same-layer blocks lie in that layer too, and X_j - X_k >= k - j for same-layer j < k.
      (a) down-window : T_ij + T_ik <= 1          for j < k < i,  k - j >= 3
          (two supports of i are within 1 of X_i, so X_j - X_k <= 2)
      (b) up-window   : T_ij + T_kj <= 1          for j < i < k,  k - i >= 3
          (two blocks resting on j are within 1 of X_j)
      (c) up-degree   : sum_{i > j} T_ij <= 3     (at most three unit blocks can rest on block j)"""
    model.update()
    T = {}
    for v in model.getVars():
        if v.VarName.startswith("T_"):
            _, a, b = v.VarName.split("_")
            T[int(a), int(b)] = v
    n_down = n_up = n_deg = 0
    for i in range(2, N + 1):
        for j in range(1, i):
            for k in range(j + 3, i):
                if (i, j) in T and (i, k) in T:
                    model.addConstr(T[i, j] + T[i, k] <= 1, name=f"WinDown_{i}_{j}_{k}")
                    n_down += 1
    for j in range(1, N + 1):
        uppers = [i for i in range(j + 1, N + 1) if (i, j) in T]
        for a_ in uppers:
            for b_ in uppers:
                if b_ - a_ >= 3:
                    model.addConstr(T[a_, j] + T[b_, j] <= 1, name=f"WinUp_{j}_{a_}_{b_}")
                    n_up += 1
        if len(uppers) > 3:
            model.addConstr(gp.quicksum(T[i, j] for i in uppers) <= 3, name=f"UpDeg_{j}")
            n_deg += 1
    model.update()
    log(f"INFO: added support-window cuts: {n_down} down-window, {n_up} up-window, {n_deg} up-degree")


def add_forced_3support(model, N, log):
    """Restrict to stacks in which at least one block touches three supports (Y3_i = 1)."""
    model.update()
    ys = []
    for i in range(4, N + 1):
        y = model.addVar(vtype=GRB.BINARY, name=f"Y3_{i}")
        model.addConstr(gp.quicksum(model.getVarByName(f"T_{i}_{j}") for j in range(1, i)) >= 3 * y,
                        name=f"Force3_{i}")
        ys.append(y)
    model.addConstr(gp.quicksum(ys) >= 1, name="AtLeastOne3Support")
    model.update()
    log(f"INFO: FORCED 3-SUPPORT model: at least one of blocks 4..{N} must touch three supports")


# ---------------------------------------------------------------------------------------------
# Structural strengthening (Sections 3 and 4 of the paper)
# ---------------------------------------------------------------------------------------------
from profiles import KNOWN_LB, CERT_UB  # noqa: E402,F401  (exact lower bounds and certified upper bounds)


def layer_capacity(n2, lmax=40):
    """Largest possible n_l under the layer-width doubling theorem n_{l+1} <= 2 n_l, n_l <= 2l - 2,
    with n_1 = 1 and n_2 fixed to n2 in {1, 2}.  Returns cumulative capacities C[l] = sum_{t<=l} n_t."""
    cap = {1: 1, 2: n2}
    for l in range(3, lmax + 1):
        # n_2 = 1: layers >= 2 form a single-root stack on block 2, so n_l <= 2(l - 1) - 2 = 2l - 4
        # (only changes anything for N >= 15)
        cap[l] = min(2 * cap[l - 1], 2 * l - 2 if n2 == 2 else 2 * l - 4)
    C, s = {}, 0
    for l in range(1, lmax + 1):
        s += cap[l]; C[l] = s
    return C


def min_layer(i, C):
    return min(l for l in C if C[l] >= i)


def apply_strengthening(model, name, N, eta, opts, log):
    """opts: set of strings from {S1, S2, S2P, S4, C1, C2, TB}.  All rows are valid for at least one
    optimal stack (S1, TB: for every optimal stack; S2: for every feasible stack; S2P: for every
    optimal stack, given the certified bound CERT_UB[N-1]; C1, C2: valid inequalities)."""
    if not opts:
        return
    model.update()
    W = 0.5; EPS = W * float(eta); w = W - EPS
    V = model.getVarByName
    rel = (name == "Mrel")
    X = {i: V(f"X_{i}") for i in range(1, N + 1)}
    L = {i: V(f"L_{i}") for i in range(1, N + 1)}
    b = {i: V(f"LayerJump_{i}") for i in range(1, N)}
    T = {(i, j): V(f"T_{i}_{j}") for i in range(2, N + 1) for j in range(1, i)}
    F = {(i, j): V(f"F_{i}_{j}") for i in range(2, N + 1) for j in range(1, i)}
    if rel:
        D = {(i, j): V(f"D_{i}_{j}") for (i, j) in T}
        MT = {(i, j): V(f"MuTop_{i}_{j}") for (i, j) in T}
        MB = {(i, j): V(f"MuBot_{i}_{j}") for (i, j) in T}
    else:
        P = {(i, j): V(f"P_{i}_{j}") for (i, j) in T}
    if any(v is None for v in list(X.values()) + list(L.values()) + list(T.values()) + list(F.values())):
        raise RuntimeError("variable names changed; cannot apply strengthening")
    if N < 5:
        raise RuntimeError("strengthening flags need N >= 5")
    added = []

    if "S2P" in opts:
        opts = set(opts) | {"S2"}
        if eta > 0:
            raise RuntimeError("S2+ is only certified for eta = 0")
        if (N - 1) not in CERT_UB or N not in KNOWN_LB:
            raise RuntimeError(f"S2+ needs a certified bound for N-1 = {N - 1}")
        margin = KNOWN_LB[N] - (CERT_UB[N - 1] + w / N)
        if margin <= 1e-3:
            raise RuntimeError(f"S2+ margin too small for N = {N}: {margin:.3e}")
        b[2].UB = 0.0
        # blocks 2 and 3 share layer 2, so L_j <= j - 1 for j >= 3
        for j in range(3, N + 1):
            L[j].UB = min(L[j].UB, j - 1)
        F[2, 1].UB = min(F[2, 1].UB, N - 2.0); F[3, 1].UB = min(F[3, 1].UB, N - 2.0)
        model.update()
        added.append(f"S2+ (n_2 = 2; base-singleton bound {CERT_UB[N - 1]:.7f} + 1/{2 * N} "
                     f"< {KNOWN_LB[N]:.7f}, margin {margin:.4f})")

    if "S1" in opts:
        # P_10 = -eps at every optimum (else shift the whole stack right); then sum_i X_i = -N eps.
        if rel:
            c = model.getConstrByName("BaseTableEdge"); c.Sense = GRB.EQUAL
        else:
            p10 = V("P_1_0"); p10.LB = -EPS; p10.UB = -EPS
            c = model.getConstrByName("Global_CoM"); c.Sense = GRB.EQUAL
        added.append("S1 (P_10 = -eps, sum X = -N eps)")

    if "S2" in opts:
        C0, C1_ = layer_capacity(2), layer_capacity(1)
        nfix = 0; lmin = {}
        for i in range(2, N + 1):
            l0, l1 = min_layer(i, C0), min_layer(i, C1_)
            L[i].LB = max(L[i].LB, l0)
            lmin[i] = max(L[i].LB, l0)
            if l1 > l0:   # n_2 = 1 (b_2 = 1) pushes block i one layer higher
                model.addConstr(L[i] >= l0 + (l1 - l0) * b[2], name=f"S2_Lmin_{i}")
        for (i, j) in T:
            if min(j, L[j].UB) + 1 < lmin[i]:   # L_j <= j (<= j - 1 under S2+), contact needs L_i = L_j + 1
                T[i, j].UB = 0.0; F[i, j].UB = 0.0; nfix += 1
        model.addConstr(T[2, 1] == 1, name="S2_T21")
        model.addConstr(T[3, 1] + b[2] == 1, name="S2_T31")
        model.addConstr(T[3, 2] == b[2], name="S2_T32")
        model.addConstr(F[2, 1] >= 1 + (N - 2) * b[2], name="S2_F21lo")
        model.addConstr(F[2, 1] <= (N - 2) + b[2], name="S2_F21hi")
        model.addConstr(F[3, 1] >= 1 - b[2], name="S2_F31lo")
        model.addConstr(F[3, 1] <= (N - 2) * (1 - b[2]), name="S2_F31hi")
        # n_2 = 2  =>  X_2 > X_1 > X_3 (endpoint exclusion + same-layer spacing)
        model.addConstr(X[2] - X[1] >= -2 * w * b[2], name="S2_sign2")
        model.addConstr(X[3] - X[1] <= 4 * w * b[2], name="S2_sign3")
        if rel:
            model.addConstr(D[2, 1] >= -2 * w * b[2], name="S2_D21")
            model.addConstr(D[3, 1] <= 0, name="S2_D31")
            model.addConstr(D[2, 1] - D[3, 1] >= 1 - (1 + 2 * w) * b[2], name="S2_D23")
        model.update()
        added.append(f"S2 (layer-width doubling: L_i lower bounds, {nfix} contact pairs fixed to 0, "
                     f"b_2 links, base force bounds, base sign rows)")

    if "TB" in opts and not rel:     # Mrel gets the tight box inside build_mrel (needs its big-M values)
        X[1].UB = min(X[1].UB, ((N - 1.0) / N) * w - EPS)   # macro-block envelope (as in build_mrel)
        if "S4" in opts and "S2P" in opts:                   # base-singleton corollary first, so the boxes below use it
            X[1].UB = min(X[1].UB, ((N - 3.0) / N) * w - EPS)
        model.update()
        xlb1 = -W; xub1 = X[1].UB
        for i in range(2, N + 1):
            lo = max(-0.5 * N, xlb1 - 2 * w * (i - 1)); hi = min(0.5 * (N - 1), xub1 + 2 * w * (i - 1))
            X[i].LB = max(X[i].LB, lo); X[i].UB = min(X[i].UB, hi)
        model.update()
        for i in range(2, N + 1):
            model.addConstr(X[i] - X[1] <= 2 * w * (L[i] - 1), name=f"TightUp_{i}")
            model.addConstr(X[1] - X[i] <= 2 * w * (L[i] - 1), name=f"TightDn_{i}")
        for (i, j) in T:
            plo = min(0.0, max(X[i].LB, X[j].LB) - w); phi = max(0.0, min(X[i].UB, X[j].UB) + w)
            plo = max(P[i, j].LB, plo); phi = min(P[i, j].UB, phi)
            P[i, j].LB = plo; P[i, j].UB = phi
            tq = V(f"Torque_{i}_{j}"); u = F[i, j].UB
            tq.LB = max(tq.LB, u * plo); tq.UB = min(tq.UB, u * phi)
        added.append("TB (tight coordinate / contact-point / torque boxes)")

    if "S4" in opts:
        # Base right-corner theorem (paper: Base Right-Corner Tightness theorem): at every optimum with
        # N >= 3, P_21 = X_1 + w.  With n_2 = 2 this gives X_1 <= (N - 3) w / N - eps (paper: S2+ base-singleton corollary).
        if N < 3:
            raise RuntimeError("S4 needs N >= 3")
        if "S2P" in opts:
            X[1].UB = min(X[1].UB, ((N - 3.0) / N) * w - EPS)
            model.update()
        if rel:
            model.addConstr(MT[2, 1] == w * F[2, 1], name="S4_BaseRightCorner")
        else:
            model.addConstr(P[2, 1] == X[1] + w, name="S4_BaseRightCorner")
            P[2, 1].LB = max(P[2, 1].LB, X[1].LB + w); P[2, 1].UB = min(P[2, 1].UB, X[1].UB + w)
        model.update()
        added.append("S4 (P_21 = X_1 + w" + (", X_1 <= (N-3)w/N - eps" if "S2P" in opts else "") + ")")

    if "C1" in opts:
        if not rel:
            raise RuntimeError("C1 is defined for Mrel only")
        n = 0
        for (i, j) in T:
            if T[i, j].UB < 0.5:
                continue
            u = F[i, j].UB; Z = MT[i, j] - MB[i, j]
            model.addConstr(Z >= u * D[i, j] + 2 * w * F[i, j] - 2 * w * u * T[i, j], name=f"C1lo_{i}_{j}")
            model.addConstr(Z <= u * D[i, j] - 2 * w * F[i, j] + 2 * w * u * T[i, j], name=f"C1hi_{i}_{j}")
            n += 2
        added.append(f"C1 ({n} lifted McCormick facets)")

    if "C2" in opts:
        if b[2].UB > 0.5:
            raise RuntimeError("C2 needs n_2 = 2 fixed (use S2P)")
        FL, FU = 1.0, float(N - 2)
        if rel:   # W23 = F21 * Delta, Delta = D21 - D31 = X2 - X3 in [1, 4w]
            Q = (MT[2, 1] - MB[2, 1]) + (MT[3, 1] - MB[3, 1]) - (N - 1) * D[3, 1]
            d = D[2, 1] - D[3, 1]; dL, dU = 1.0, 4 * w
        else:     # Q23 = tau10 - X1 - (N-1) P31 = F21 * delta, delta = P21 - P31 in [1 - 2w, 2w]
            Q = V("Torque_1_0") - X[1] - (N - 1) * P[3, 1]
            d = P[2, 1] - P[3, 1]; dL, dU = 1.0 - 2 * w, 2 * w
        f = F[2, 1]
        model.addConstr(Q >= FL * d + dL * f - FL * dL, name="C2_a")
        model.addConstr(Q >= FU * d + dU * f - FU * dU, name="C2_b")
        model.addConstr(Q <= FU * d + dL * f - FU * dL, name="C2_c")
        model.addConstr(Q <= FL * d + dU * f - FL * dU, name="C2_d")
        added.append("C2 (base-pair difference RLT, 4 rows)")

    model.update()
    for a in added:
        log(f"INFO: strengthening {a}")


def build_mrel(model, N, eta=0.0, rightmost=None, layer_lb=None, tight=False):
    """M_rel: centre-relative moment reformulation .
    Same feasible set as M1u (layer-ordered upper DAG, base linearisation, force capacities,
    geometric bound sum_j T_ij <= 3), but contact points P_ij and torques F_ij*P_ij are replaced by
      d_ij  = X_i - X_j (if T_ij = 1, else 0)          in [-2w, 2w]
      mub_ij = F_ij (P_ij - X_i), mut_ij = F_ij (P_ij - X_j)  with |mub|, |mut| <= w F_ij
      mut_ij - mub_ij = F_ij * d_ij                      (the only bilinear term, box width 4w)
    and the moment balance about each block centre: sum_j mub_ij = sum_k mut_ki.
    w = W - eps is the effective half width.  Base: F_10 = N, mub_10 = N (P_10 - X_1), P_10 <= -eps.
    Optional decomposition (rightmost = k): block k is a rightmost block, i.e. X_i <= X_k for all i,
    k is the first block of its layer (LayerJump_{k-1} = 1), objective X_k + W, and optionally
    L_k >= layer_lb (valid if layer_lb <= ceil(LB + 1/(2N)) for a known feasible value LB)."""
    W = 0.5; EPS = W * float(eta); w = W - EPS; MASS = 1.0
    # Coordinate box, valid for optimal stacks:
    #   X_1 >= -W: at an optimum P_10 = -eps (else shift right), and P_10 >= X_1 - w.
    #   X_1 <= ((N-1)/N) w - eps: base envelope (paper: Support Degrees and Envelopes lemma).
    #   X_i <= (N-1)/2: shortest-path argument from the rightmost block + sum X <= 0.
    #   X_i >= X_1 - (N-1)/2 >= -N/2: same argument applied to each balanced component resting on block 1.
    xlb = {1: -W}; xub = {1: ((N - 1.0) / N) * w - EPS}
    for i in range(2, N + 1):
        xlb[i] = -0.5 * float(N); xub[i] = 0.5 * float(N - 1)
        if tight:
            # Tight box (paper: Support Degrees and Envelopes lemma): L_i <= i and |X_p - X_q| <= 2w along a support
            # chain down to block 1, so |X_i - X_1| <= 2w (L_i - 1) <= 2w (i - 1).
            xlb[i] = max(xlb[i], xlb[1] - 2 * w * (i - 1))
            xub[i] = min(xub[i], xub[1] + 2 * w * (i - 1))
    MAX_X = 0.5 * float(N)
    X = {i: model.addVar(lb=xlb[i], ub=xub[i], name=f"X_{i}") for i in range(1, N + 1)}
    L = {1: model.addVar(lb=1, ub=1, vtype=GRB.INTEGER, name="L_1")}
    for i in range(2, N + 1):
        L[i] = model.addVar(lb=2, ub=i, vtype=GRB.INTEGER, name=f"L_{i}")
    if tight:
        for i in range(2, N + 1):
            model.addConstr(X[i] - X[1] <= 2 * w * (L[i] - 1), name=f"TightUp_{i}")
            model.addConstr(X[1] - X[i] <= 2 * w * (L[i] - 1), name=f"TightDn_{i}")
    T = {}; F = {}; D = {}; MB = {}; MT = {}
    # |mub10| <= (N-1) w from MomBal_1 and |mut_k1| <= w F_k1, sum_k F_k1 = N - 1.
    mub10 = model.addVar(lb=-(N - 1) * w, ub=(N - 1) * w, name="MuBot_1_0")
    model.addConstr(mub10 + N * X[1] <= -N * EPS, name="BaseTableEdge")   # P_10 <= -eps
    # Global identity sum_i X_i = mub10 + N X_1 (sum of all centre moment balances + force balance);
    # valid in the exact model but NOT implied by the LP relaxation, so it is added explicitly.
    model.addConstr(gp.quicksum(X[i] for i in range(1, N + 1)) == mub10 + N * X[1], name="GlobalTorqueCoupling")
    for i in range(2, N + 1):
        u = float(N - i + 1) * MASS
        for j in range(1, i):
            T[i, j] = model.addVar(vtype=GRB.BINARY, name=f"T_{i}_{j}")
            F[i, j] = model.addVar(lb=0.0, ub=u, name=f"F_{i}_{j}")
            D[i, j] = model.addVar(lb=-2 * w, ub=2 * w, name=f"D_{i}_{j}")
            MB[i, j] = model.addVar(lb=-u * w, ub=u * w, name=f"MuBot_{i}_{j}")
            MT[i, j] = model.addVar(lb=-u * w, ub=u * w, name=f"MuTop_{i}_{j}")
    for i in range(2, N + 1):
        u = float(N - i + 1) * MASS
        model.addConstr(gp.quicksum(T[i, j] for j in range(1, i)) >= 1, name=f"DAGSupMin_{i}")
        max_sup = 2 if eta > 0.0 else 3          # support-degree lemma (Section 3.1): <= 2 when eta > 0
        model.addConstr(gp.quicksum(T[i, j] for j in range(1, i)) <= max_sup, name=f"DAGSupMax_{i}")
        for j in range(1, i):
            model.addConstr(F[i, j] <= u * T[i, j])
            model.addConstr(MB[i, j] <= w * F[i, j]); model.addConstr(MB[i, j] >= -w * F[i, j])
            model.addConstr(MT[i, j] <= w * F[i, j]); model.addConstr(MT[i, j] >= -w * F[i, j])
            dlo = xlb[i] - xub[j]; dhi = xub[i] - xlb[j]
            model.addConstr(D[i, j] <= 2 * w * T[i, j]); model.addConstr(D[i, j] >= -2 * w * T[i, j])
            model.addConstr(D[i, j] - (X[i] - X[j]) <= -dlo * (1 - T[i, j]))
            model.addConstr(D[i, j] - (X[i] - X[j]) >= -dhi * (1 - T[i, j]))
            model.addConstr(MT[i, j] - MB[i, j] == F[i, j] * D[i, j], name=f"Bilinear_{i}_{j}")
            model.addConstr(L[i] - L[j] <= 1 + (i - 1) * (1 - T[i, j]))
            model.addConstr(L[i] - L[j] >= 1 - j * (1 - T[i, j]))
    b = {}
    for i in range(1, N):
        model.addConstr(L[i] <= L[i + 1], name=f"SymLayer_{i}")
        b[i] = model.addVar(vtype=GRB.BINARY, name=f"LayerJump_{i}")
        model.addConstr(L[i + 1] - L[i] >= b[i]); model.addConstr(L[i + 1] - L[i] <= i * b[i])
        if tight:   # X_i - X_{i+1} >= 2W - big_m always holds
            big_m = max(0.0, 2 * W + xub[i + 1] - xlb[i])
        else:       # original constant (kept so that untightened Mrel runs reproduce earlier models exactly)
            big_m = (xub[2] - xlb[2]) + 2.0 * W
        model.addConstr(X[i] - X[i + 1] >= 2 * W - big_m * b[i], name=f"NonOverlap_{i}")
    # force balance (block 1 carries everything: F_10 = N)
    model.addConstr(gp.quicksum(F[k, 1] for k in range(2, N + 1)) + MASS == float(N), name="ForceBal_1")
    model.addConstr(gp.quicksum(MT[k, 1] for k in range(2, N + 1)) == mub10, name="MomBal_1")
    for i in range(2, N + 1):
        model.addConstr(gp.quicksum(F[i, j] for j in range(1, i)) ==
                        MASS + gp.quicksum(F[k, i] for k in range(i + 1, N + 1)), name=f"ForceBal_{i}")
        model.addConstr(gp.quicksum(MB[i, j] for j in range(1, i)) ==
                        gp.quicksum(MT[k, i] for k in range(i + 1, N + 1)), name=f"MomBal_{i}")
    if rightmost is None:
        re = [model.addVar(lb=-MAX_X, ub=MAX_X, name=f"RE_{i}") for i in range(1, N + 1)]
        for i in range(1, N + 1):
            model.addConstr(re[i - 1] == X[i] + W)
        mo = model.addVar(lb=0.5 - 2.0 * EPS, ub=MAX_X, name="max_overhang")
        model.addGenConstrMax(mo, re, name="LinkMaxEdge")
        model.setObjective(mo, GRB.MAXIMIZE)
    else:
        k = rightmost
        for i in range(1, N + 1):
            if i != k:
                model.addConstr(X[i] <= X[k], name=f"Rightmost_{i}")
        if k >= 2:
            model.addConstr(b[k - 1] == 1, name="RightmostFirstInLayer")
        if layer_lb is not None and k >= 2:
            model.addConstr(L[k] >= layer_lb, name="RightmostLayerLB")
        model.setObjective(X[k] + W, GRB.MAXIMIZE)
    model.update()


def build_model(name, N, eta, log, rightmost=None, layer_lb=None, tight=False):
    if name == "Mrel":
        env = gp.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = gp.Model(f"Mrel_certify_N{N}", env=env)
        build_mrel(model, N, eta, rightmost=rightmost, layer_lb=layer_lb, tight=tight)
        if tight:
            log("INFO: Mrel with tight coordinate box |X_i - X_1| <= 2w (L_i - 1)")
        return model
    if name == "M0":
        env = gp.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = gp.Model(f"M0_certify_N{N}", env=env)
        base.build_m0(model, N, eta)
        return model

    captured = {}

    def _capture(model, *args, **kwargs):
        captured["model"] = model
        return None

    if name in ("M1u", "M1", "M1w"):
        import ablation_m1_order_prune as mod
        solve = mod.solve_m1_order_prune
    else:
        raise SystemExit(f"model {name} needs model files that are not part of this repository; use M0, M1u or Mrel")

    original = mod.execute_and_report
    mod.execute_and_report = _capture
    try:
        solve(N, 0, 42, threads=1, eta=eta)
    finally:
        mod.execute_and_report = original

    model = captured["model"]
    if name in CONDITIONAL_MODELS:
        log("WARNING: two-support bound kept (sum_j T_ij <= 2): results are CONDITIONAL "
            "(valid only within the two-support subclass)")
    else:
        relax_two_support(model, N, log)
    if name.endswith("w"):
        add_window_cuts(model, N, log)
    return model


def main():
    ap = argparse.ArgumentParser(description="Long-horizon unconditional certification (M0 / M1u / M2u / M3u; M1 / M2 / M3 = conditional originals)")
    ap.add_argument("N", type=int)
    ap.add_argument("--model", choices=MODELS, default="M0")
    ap.add_argument("--time-limit", type=float, default=259200, help="seconds (default 72 h)")
    ap.add_argument("--threads", "-t", type=int, default=4)
    ap.add_argument("--gap", type=float, default=1e-6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--eta", type=float, default=0.0)
    ap.add_argument("--soft-mem-limit", type=float, default=6.0, help="GB")
    ap.add_argument("--nodefile-start", type=float, default=2.0, help="GB")
    ap.add_argument("--log-interval", type=float, default=60.0)
    ap.add_argument("--start-json", type=str, default=None)
    ap.add_argument("--outdir", type=str, default="results_certify")
    ap.add_argument("--force-3support", action="store_true",
                    help="only stacks in which some block touches three supports (lemma probe)")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="Gurobi Cutoff: status CUTOFF proves no solution better than this value")
    ap.add_argument("--strengthen", type=str, default="",
                    help="comma list from S1,S2,S2P,S4,TB,C1,C2 (M1u / Mrel only), e.g. S1,S2P,S4,C2")
    ap.add_argument("--param", action="append", default=[],
                    help="extra Gurobi parameter NAME=VALUE (repeatable), e.g. --param MIPFocus=3")
    ap.add_argument("--ptag", type=str, default="",
                    help="label for --param / --branch-priority settings, appended to the run tag")
    ap.add_argument("--branch-priority", action="store_true",
                    help="BranchPriority attribute: LayerJump_i (i<=6) = 20-i, L_i = 10")
    ap.add_argument("--no-callback", action="store_true",
                    help="optimize without the Python callback (no trace/incumbent files)")
    args = ap.parse_args()
    opts = {s.strip().upper() for s in args.strengthen.split(",") if s.strip()}
    bad = opts - {"S1", "S2", "S2P", "S4", "TB", "C1", "C2"}
    if bad:
        sys.exit(f"ERROR: unknown --strengthen items {sorted(bad)}")
    if opts and args.model not in ("M1u", "Mrel"):
        sys.exit("ERROR: --strengthen is implemented for --model M1u and --model Mrel only")
    if "C1" in opts and args.model != "Mrel":
        sys.exit("ERROR: C1 (lifted McCormick facets) is defined for --model Mrel only")
    if "C2" in opts and "S2P" not in opts:
        sys.exit("ERROR: C2 needs S2P (n_2 = 2 fixed)")

    if args.model in ("M2u", "M3u", "M2w", "M3w", "M2", "M3") and args.N >= 14:
        sys.exit("ERROR: M2u/M3u use the spinal-based bound D*(N-1) for N >= 14, so they are not "
                 "unconditional there. Use --model M0 for N >= 14.")

    os.makedirs(args.outdir, exist_ok=True)
    eta_sfx = f"_eta{args.eta:.2f}" if args.eta > 0 else ""
    f3 = "F3" if args.force_3support else ""
    order = ["TB", "S1", "S2", "S2P", "S4", "C1", "C2"]
    stag = "".join("-" + o for o in order if o in opts)
    ptag = ("-" + args.ptag) if args.ptag else ""
    tag = f"{args.model}{stag}{ptag}{f3}cert_N{args.N}{eta_sfx}_s{args.seed}"
    tag_path = os.path.join(args.outdir, tag)
    nodedir = os.path.join(args.outdir, f"nodefiles_{tag}")
    os.makedirs(nodedir, exist_ok=True)

    log = base.TeeLog(f"{tag_path}.out")
    signal.signal(signal.SIGTERM, base._request_stop)
    signal.signal(signal.SIGINT, base._request_stop)

    model = build_model(args.model, args.N, args.eta, log, tight=("TB" in opts))
    apply_strengthening(model, args.model, args.N, args.eta, opts, log)
    if args.force_3support:
        if args.model in ("M0", "M1", "M2", "M3"):
            sys.exit("ERROR: --force-3support needs a model with the <= 3 bound (M1u/M2u/M3u/M1w/M2w/M3w)")
        add_forced_3support(model, args.N, log)

    model.setParam("OutputFlag", 1)
    model.setParam("LogToConsole", 0)
    model.setParam("LogFile", f"{tag_path}_gurobi.log")
    model.setParam("NonConvex", 2)
    model.setParam("Threads", args.threads)
    model.setParam("MIPGap", args.gap)
    model.setParam("Seed", args.seed)
    model.setParam("TimeLimit", float(args.time_limit))
    model.setParam("NodefileStart", args.nodefile_start)
    model.setParam("NodefileDir", nodedir)
    model.setParam("SoftMemLimit", args.soft_mem_limit)
    for kv in args.param:
        k, v = kv.split("=", 1)
        model.setParam(k.strip(), float(v) if any(ch in v for ch in ".eE") else int(v))
        log(f"INFO: Gurobi parameter {k.strip()} = {v}")
    if args.branch_priority:
        for i in range(1, args.N):
            var = model.getVarByName(f"LayerJump_{i}")
            if var is not None and i <= 6:
                var.BranchPriority = 20 - i
        for i in range(2, args.N + 1):
            var = model.getVarByName(f"L_{i}")
            if var is not None:
                var.BranchPriority = 10
        model.update()
        log("INFO: BranchPriority set (LayerJump_i = 20 - i for i <= 6, L_i = 10)")
    if args.cutoff is not None:
        model.setParam("Cutoff", args.cutoff)
        log(f"INFO: Cutoff = {args.cutoff} (status CUTOFF => no solution better than the cutoff exists)")

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
    log(f"META: Method={args.model}_Certify | N={args.N} | Eta={args.eta:.4f} | Seed={args.seed} | "
        f"TimeLimit={args.time_limit:.0f}s | Threads={args.threads} | MIPGap={args.gap}")
    log(f"  Host={socket.gethostname()} | Gurobi={ver} | Python={platform.python_version()} | "
        f"Start={datetime.now():%F %T} | PID={os.getpid()}")
    log(f"  SoftMemLimit={args.soft_mem_limit} GB | NodefileStart={args.nodefile_start} GB | "
        f"Vars={model.NumVars} | LinConstrs={model.NumConstrs} | QConstrs={model.NumQConstrs}")
    log("=" * 80)

    t0 = time.time()
    cb = base.CertifyCallback(args.N, model.getVars(), tag_path, t0, args.log_interval, log)
    try:
        if args.no_callback:
            log("INFO: optimizing WITHOUT the Python callback")
            model.optimize()
        else:
            model.optimize(cb)
    except gp.GurobiError as e:
        log(f"[{datetime.now():%F %T}] GUROBI ERROR: {e}")
        cb.close()
        log.close()
        base._cleanup_nodedir(nodedir)
        sys.exit(2)
    finally:
        if not cb.trace_f.closed:
            cb.close()
    runtime = time.time() - t0

    status = "CUTOFF" if model.Status == GRB.CUTOFF else base.STATUS_NAMES.get(model.Status, f"STATUS_{model.Status}")
    has_sol = model.SolCount > 0
    obj = model.ObjVal if has_sol else float("nan")
    bnd = model.ObjBound
    gap_pct = model.MIPGap * 100.0 if has_sol else float("nan")
    nodes = int(model.NodeCount)

    stack = None
    if has_sol:
        val = {v.VarName: v.X for v in model.getVars()}
        stack = base.stack_from_values(args.N, val)
        model.write(f"{tag_path}.sol")

    base.atomic_write_json(f"{tag_path}_result.json", {
        "method": f"{args.model}_Certify (" + ("CONDITIONAL: two-support subclass" if args.model in CONDITIONAL_MODELS else "unconditional") + ")", "N": args.N, "eta": args.eta,
        "seed": args.seed, "threads": args.threads, "mip_gap_target": args.gap,
        "strengthen": sorted(opts), "params": args.param, "branch_priority": args.branch_priority,
        "callback": not args.no_callback, "cutoff": args.cutoff,
        "time_limit_sec": args.time_limit, "gurobi_version": ver, "host": socket.gethostname(),
        "status": status, "objective": obj, "bound": bnd, "gap_pct": gap_pct,
        "runtime_sec": runtime, "nodes": nodes, "sol_count": model.SolCount,
        "first_sol_sec": cb.first_sol_time,
        "milestones": {f"{k:g}": m for k, m in sorted(cb.milestones.items(), reverse=True)},
        "tolerances": {"FeasibilityTol": model.Params.FeasibilityTol,
                       "IntFeasTol": model.Params.IntFeasTol, "OptimalityTol": model.Params.OptimalityTol},
        "stack": stack,
    })

    log("")
    log("=" * 80)
    log("SUMMARY_REPORT:")
    log(f"  Method           : {args.model}_Certify")
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
    base._cleanup_nodedir(nodedir)


if __name__ == "__main__":
    main()
