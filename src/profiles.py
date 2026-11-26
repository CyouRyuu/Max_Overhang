"""Solver-independent helpers shared by the certification scripts and the auditors.

Nothing here imports Gurobi, so the auditors (pd_certificate.py, ns_audit.py) and the exact checks run without a
solver licence.  The functions were moved unchanged from profile_run.py (enumerate_profiles, layers_of,
energy_profile_check, spinal_bound), ns_check.py (spinal_opt) and certify_run.py (KNOWN_LB, CERT_UB), which
re-export them under their old names.
"""

# Exactly verified lower bounds z_Q(N), truncated at the 7th decimal (main certificate table of the paper; exact
# fractions in data/rational_stacks.json), and certified upper bounds z_bar(N) of the unified tight chain (S2+ chain
# table of the paper:
# N = 4 from monolithic M0 at relative gap 1e-9, N >= 5 from the spinal-cut layer-profile certificates).
# They are used only by the base-singleton precondition of S2P,  CERT_UB[N - 1] + w / N < KNOWN_LB[N],  which never
# enters the optimization model.  (The original runs checked the same precondition with earlier, valid but looser
# constants; the models and hence all results are identical.)
KNOWN_LB = {4: 1.1678932, 5: 1.3045548, 6: 1.4367002, 7: 1.5300452, 8: 1.6315068, 9: 1.7152728,
            10: 1.7871258, 11: 1.8587795, 12: 1.9250913, 13: 1.9845089, 14: 2.0382169, 15: 2.0928962,
            16: 2.1438427}
CERT_UB = {4: 1.1678933, 5: 1.3045649, 6: 1.4367103, 7: 1.5300553, 8: 1.6315169, 9: 1.7152829,
           10: 1.7871359, 11: 1.8587896, 12: 1.9251014, 13: 1.9845190, 14: 2.0382270, 15: 2.0929064,
           16: 2.1438528}


def enumerate_profiles(N, n2_fixed=2):
    out = []

    def rec(prof, s):
        if s == N:
            out.append(tuple(prof)); return
        l = len(prof) + 1
        hi = min(2 * prof[-1], 2 * l - 2, N - s)
        for n in range(1, hi + 1):
            if l == 2 and n2_fixed is not None and n != n2_fixed:
                continue
            rec(prof + [n], s + n)
    rec([1], 1)
    return out


def layers_of(profile):
    lay, blocks, k = {}, {}, 1
    for l, n in enumerate(profile, start=1):
        blocks[l] = list(range(k, k + n))
        for i in blocks[l]:
            lay[i] = l
        k += n
    return lay, blocks


def energy_profile_check(N, profile, eta=0.0):
    """Second-moment energy test (paper: second-moment energy theorem and its two corollaries).
    Every feasible stack, translated so that P_10 = -eps, satisfies
        sum_b (X_b + eps)^2 + N (X_1 + eps)^2 <= w^2 sum_b (L_b - 1) =: E(n),
    and the layer of m blocks has within-layer second moment >= m(m^2-1)/12.  So dE < 0 => no feasible
    stack with this profile.  layer_ubs[l] bounds X_first(l) + W for stacks with sum_b X_b = -N eps (S1)."""
    import math
    W = 0.5; w = W * (1.0 - float(eta))
    E = w * w * sum((l - 1) * n for l, n in enumerate(profile, start=1))
    dE = E - sum(n * (n * n - 1) / 12.0 for n in profile)
    if dE < -1e-12:
        return False, dE, {}
    dE = max(0.0, dE); ubs = {}
    for l in range(2, len(profile) + 1):
        m = profile[l - 1]
        Am = m + m * m / ((N - m - 1.0) + 1.0 / (N + 1.0))
        if m == 1 or 2.0 * math.sqrt(Am * dE) <= m:
            ubs[l] = w + 0.5 * (m - 1) + math.sqrt(dE / Am)
        else:
            ubs[l] = w + math.sqrt((1.0 / Am + (m - 1.0) / m) * (dE + 0.25 * m * (m - 1)))
    return True, dE, ubs


def spinal_bound(r, F):
    """Paterson-Zwick loaded-spine optimum S*_r(F) (PZ 2009, Lemmas 3.4-3.5), in units of 2w.
    Unconstrained version (no t_i - t_{i-1} >= 1): S^_r(F) = max_{u_1..u_{r-1} > 0, u_r = F}
    sum_s [1 - (u_{s-1} + 1/2) / u_s], u_0 = 0.  Strictly convex in log u, unique stationary point
    u_{s+1} = u_s^2 / (u_{s-1} + 1/2); u_r is increasing in u_1, so bisection on u_1.
    Paper, spinal-cut proposition (family SP): if s_2..s_r are each supported ONLY by the rightmost
    block of the layer below, then X_{s_r} + W <= 2w * S^_r(N)."""
    F = float(F)
    if r <= 1:
        return 1.0 - 0.5 / F

    def last(u1):
        a, b = 0.0, u1
        for _ in range(r - 1):
            a, b = b, b * b / (a + 0.5)
        return b
    lo, hi = 0.0, 1.0
    while last(hi) < F:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if last(mid) < F:
            lo = mid
        else:
            hi = mid
    us = [0.0, 0.5 * (lo + hi)]
    for _ in range(r - 1):
        us.append(us[-1] * us[-1] / (us[-2] + 0.5))
    us[-1] = F
    return sum(1.0 - (us[s - 1] + 0.5) / us[s] for s in range(1, r + 1))


def spinal_opt(N):
    """S*(N) / (2w): the loaded-spine optimum, maximized over the spine height."""
    return max(spinal_bound(r, N) for r in range(2, N + 1))
