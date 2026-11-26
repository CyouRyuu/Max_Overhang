"""Independent exact checks that need no optimization solver.

1. Re-checks every stack built by verify_rational_stacks.reconstruct_* with a separate, minimal equilibrium
   checker (eta = 0, W = 1/2, table edge at 0, single root, layered; each contact is one resultant force
   F >= 0 applied at P), in exact rational arithmetic.
2. Distances of the exact overhangs z_Q(N) to the loaded-spine optimum S*(N) (40-digit mpmath).
3. The acyclic base-singleton (S2+) chain of the paper's S2+ table: theta_N = zbar(N-1) + 1/(2N) < z_Q(N).
4. The non-spinal chain: gamma_ns(N) = 0.99 S*(N) rounded upward at the 7th decimal, general condition theta_N < gamma_ns(N) and the sharper
   condition gamma_ns(N-1) + 1/(2N) < gamma_ns(N).
5. The energy pruning used before any model is built: the profile test is exact (energy margins are multiples
   of 1/4, every pruned profile has margin <= -1/4), and every candidate layer removed by the base envelope or
   by the energy layer bound at the cutoffs of the paper is removed with a margin far above rounding errors.
6. The loaded-spine constants: the double-precision values of S^_r(N) used by the cuts differ from a 60-digit
   evaluation by far less than the safety margin 1e-9, and for N <= 16 the optimal loads of the maximizing height
   have increments >= 1, so S^(N) = S*(N).
7. Directed rounding: the constants used in the cuts exceed rigorous interval upper bounds of 2w S^_r(N).

  python3 src/check_exact_stacks.py
"""
import os
import sys
import math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_rational_stacks as v
from mpmath import mp, mpf
mp.dps = 40


def shat(r, F):
    """Unconstrained loaded-spine optimum (Paterson-Zwick recursion u_{s+1} = u_s^2/(u_{s-1}+1/2)), 40 digits."""
    F = mpf(F)
    if r == 1:
        return 1 - 1 / (2 * F), None

    def roll(u1):
        us = [mpf(0), u1]
        for s in range(1, r):
            us.append(us[s] ** 2 / (us[s - 1] + mpf('0.5')))
        return us
    lo, hi = mpf(0), mpf(1)
    while roll(hi)[-1] < F:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if roll(mid)[-1] < F:
            lo = mid
        else:
            hi = mid
    us = roll((lo + hi) / 2); us[-1] = F
    return sum(1 - (us[s - 1] + mpf('0.5')) / us[s] for s in range(1, r + 1)), us

def check(N, layers, X, contacts):
    W = Fr(1, 2)
    assert sorted(X) == list(range(1, N + 1))
    # layer membership and spacing (sort each layer by X, not by id)
    by = {}
    for i, l in layers.items(): by.setdefault(l, []).append(i)
    assert by[1] == [1], "single root"
    for l, bs in by.items():
        xs = sorted((X[i] for i in bs), reverse=True)
        for a, b in zip(xs, xs[1:]): assert a - b >= 1, ("overlap", l)
    # contacts
    Fin = {i: Fr(0) for i in X}; Tin = {i: Fr(0) for i in X}
    Fout = {i: Fr(0) for i in X}; Tout = {i: Fr(0) for i in X}
    supp = set()
    for (i, j), (F, tau) in contacts.items():
        assert F >= 0
        if F == 0: continue
        P = tau / F
        if j == 0:
            assert i == 1 and P <= 0 and X[1] - W <= P <= X[1] + W
        else:
            assert layers[i] == layers[j] + 1
            assert max(X[i], X[j]) - W <= P <= min(X[i], X[j]) + W
            Fout[j] += F; Tout[j] += tau
        Fin[i] += F; Tin[i] += tau; supp.add(i)
    assert supp == set(X), "every block supported"
    for i in X:
        assert Fin[i] - Fout[i] == 1, ("force", i)
        assert Tin[i] - Tout[i] == X[i], ("moment", i)   # moment about x=0: weight 1 at X_i
    return max(X[i] + W for i in X)

def min_loaded_overlap(X, contacts):
    """Smallest overlap length over contacts with positive force (table: overlap of the base with x <= 0).
    Positive means no force passes through a point contact, as required by Paterson and Zwick."""
    ov = []
    for (i, j), (F, tau) in contacts.items():
        if F > 0:
            ov.append(-(X[i] - Fr(1, 2)) if j == 0 else 1 - abs(X[i] - X[j]))
    return min(ov)


def spinal_both(layers, X):
    """Spinality of the stack under the closed resting relation of the paper (|X_i - X_j| <= 1 between adjacent
    layers) and under Paterson and Zwick's strict relation (< 1).  Returns (closed, strict)."""
    zmax = max(X.values())
    lstar = min(layers[i] for i in X if X[i] == zmax)
    p = next(i for i in X if layers[i] == lstar and X[i] == zmax)
    out = []
    for strict in (False, True):
        S, todo = {p}, [p]
        while todo:
            i = todo.pop()
            for j in X:
                if layers[j] == layers[i] - 1 and j not in S and \
                        (abs(X[i] - X[j]) < 1 if strict else abs(X[i] - X[j]) <= 1):
                    S.add(j); todo.append(j)
        per = [sum(1 for j in S if layers[j] == t) for t in range(1, lstar + 1)]
        out.append(all(c == 1 for c in per))
    return tuple(out)


res, overlap, spinal = {}, {}, {}
specs = sorted(v.SPECS_4_9_SPINAL + v.SPECS_10_12 + v.SPECS_13_16)
for N, k, tl, ys, z1 in specs:
    oh, prof, layers, X, contacts, roles = v.reconstruct_spinal(N, k, tl, ys, z1)
    res[N] = (check(N, layers, X, contacts), prof)
    overlap[N] = min_loaded_overlap(X, contacts)
    spinal[N] = spinal_both(layers, X)
oh, prof, layers, X, contacts, roles = v.reconstruct_n5_bridgeless()
res[5] = (check(5, layers, X, contacts), prof)
overlap[5] = min_loaded_overlap(X, contacts)
spinal[5] = spinal_both(layers, X)
assert all(o > 0 for o in overlap.values()), overlap
assert all(a and b for a, b in spinal.values()), spinal
# Certified upper bounds zbar(N) (paper, main certificate table; identical to certify_run.CERT_UB).
CERT_UB = {4: "1.1678933", 5: "1.3045649", 6: "1.4367103", 7: "1.5300553", 8: "1.6315169", 9: "1.7152829",
           10: "1.7871359", 11: "1.8587896", 12: "1.9251014", 13: "1.9845190", 14: "2.0382270", 15: "2.0929064",
           16: "2.1438528"}

print("== 1-2. exact stacks (independent equilibrium check passed for every N)")
Sstar = {N: max(shat(r, N)[0] for r in range(2, N + 1)) for N in range(4, 17)}
for N in sorted(res):
    oh, prof = res[N]
    print(f"N={N:2d} {str(prof):34s} z_Q={float(oh):.10f} trunc7={math.floor(oh * 10**7) / 10**7:.7f} "
          f"S*-z_Q={float(Sstar[N] - mpf(oh.numerator) / oh.denominator):+.2e}  zbar-z_Q={float(Fr(CERT_UB[N]) - oh):.2e}")
print("N=4 closed form (15-4sqrt2)/8 =", (15 - 4 * 2 ** 0.5) / 8)
print("every stack is spinal under the closed resting relation and in the sense of Paterson and Zwick")
print(f"every loaded contact has positive overlap; minimum {float(min(overlap.values())):.6f} (N = {min(overlap, key=overlap.get)})")

print("== 3. base-singleton chain: theta_N = zbar(N-1) + 1/(2N) < z_Q(N)")
worst = None
for N in range(5, 17):
    theta = Fr(CERT_UB[N - 1]) + Fr(1, 2 * N)
    m = res[N][0] - theta
    assert m > 0, N
    worst = min(worst, (m, N)) if worst else (m, N)
    print(f"N={N:2d} theta={float(theta):.7f} z_Q={float(res[N][0]):.10f} margin={float(m):+.6f}")
print(f"minimum margin {float(worst[0]):.6f} at N = {worst[1]}")

print("== 4. non-spinal chain at rho = 1e-2, gamma_ns(N) = ceil_7(0.99 S*(N))")
from mpmath import ceil as mpceil
g = {N: mpceil(mpf("0.99") * Sstar[N] * 10**7) / 10**7 for N in range(4, 17)}
worst = None
for N in range(5, 17):
    gen = g[N] - (mpf(CERT_UB[N - 1]) + mpf(1) / (2 * N))
    sharp = g[N] - (g[N - 1] + mpf(1) / (2 * N)) if N >= 6 else None
    best = max(gen, sharp) if sharp is not None else gen
    assert best > 0, N
    worst = min(worst, (best, N)) if worst else (best, N)
    stxt = "" if sharp is None else f" sharper={float(sharp):+.6f}"
    print(f"N={N:2d} gamma_ns={float(g[N]):.7f} general={float(gen):+.6f}{stxt}")
print(f"minimum usable margin {float(worst[0]):.6f} at N = {worst[1]}")
print("S*(20) =", max((float(shat(r, 20)[0]), r) for r in range(2, 21)))
import profiles as pf
CUT = {N: float(CERT_UB[N]) for N in range(5, 17)}
neg, lay = None, None
for N in range(5, 17):
    for p in pf.enumerate_profiles(N, 2):
        dE = Fr(1, 4) * sum((l - 1) * n for l, n in enumerate(p, start=1)) - sum(Fr(n * (n * n - 1), 12) for n in p)
        ok, _, ubs = pf.energy_profile_check(N, p)
        assert ok == (dE >= 0), (N, p)
        if dE < 0:
            neg = dE if neg is None else max(neg, dE)
            continue
        for l in range(2, len(p) + 1):
            base = (N - 3) / (2 * N) + 0.5 + (l - 1) - CUT[N]
            m = -base if base <= 0 else (CUT[N] - ubs[l] if ubs[l] <= CUT[N] else None)
            if m is not None:
                lay = m if lay is None else min(lay, m)
assert neg <= Fr(-1, 4) and lay > 1e-6
print(f"== 5. energy pruning: largest pruned profile margin {neg}; smallest margin of a removed candidate layer {lay:.6f}")
from mpmath import mp
mp.dps = 60


def shat_mp(r, F):
    """S^_r(F) and its maximizing loads at 60 digits (same recurrence and bisection as profiles.spinal_bound)."""
    F = mpf(F)
    if r <= 1:
        return 1 - mpf(1) / 2 / F, [mpf(0), F]

    def last(u1):
        a, b = mpf(0), u1
        for _ in range(r - 1):
            a, b = b, b * b / (a + mpf(1) / 2)
        return b
    lo, hi = mpf(0), mpf(1)
    while last(hi) < F:
        hi *= 2
    for _ in range(400):
        mid = (lo + hi) / 2
        if last(mid) < F:
            lo = mid
        else:
            hi = mid
    us = [mpf(0), (lo + hi) / 2]
    for _ in range(r - 1):
        us.append(us[-1] ** 2 / (us[-2] + mpf(1) / 2))
    us[-1] = F
    return sum(1 - (us[t - 1] + mpf(1) / 2) / us[t] for t in range(1, r + 1)), us


dev = max(abs(mpf(pf.spinal_bound(r, N)) - shat_mp(r, N)[0]) for N in range(4, 21) for r in range(1, N + 1))
inc = min(min(b - a for a, b in zip(u[:-1], u[1:]))
          for u in (shat_mp(max(range(2, N + 1), key=lambda r: shat_mp(r, N)[0]), N)[1] for N in range(4, 17)))
assert dev < 1e-12 and inc >= 1
print(f"== 6. loaded-spine constants: max |double - 60 digits| = {float(dev):.1e} (margin used 1e-9); "
      f"smallest load increment at the maximizing height, N <= 16: {float(inc):.3f}")
# 7. Directed check with rigorous interval enclosures (mpmath.iv, outward rounding; the enclosure of Appendix C):
#    the constant written into the cuts, 2w S^_r(N) + 1e-9 in double precision (2w = 1 at eta = 0), exceeds the
#    upper end of the enclosure of 2w S^_r(N) for every r <= N <= 20.
from spine_faces import shat_enclosure  # noqa: E402
exc = min(mpf(pf.spinal_bound(r, N) + 1e-9) - mpf(shat_enclosure(r, N).b) for N in range(4, 21) for r in range(1, N + 1))
assert exc > 0
print(f"== 7. directed rounding of the loaded-spine constants: the constants used exceed rigorous upper bounds "
      f"of 2w S^_r(N) for all r <= N <= 20, by at least {float(exc):.2e}")
print("ALL CHECKS PASSED")
