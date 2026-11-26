"""Independent check of the exact spinal optima (Appendix C): S(N) = S*(N) for N = 4, 6, 8..19 and the
exception values at N = 5, 7.  Needs sympy and mpmath; no optimization solver.

  python3 src/check_algebraic_spines.py [--nmax 16]

For each N, the stack is the shielded construction of Section 5.2 (Table A1 seeds; same formulas as
verify_rational_stacks.py) with the exact loads of the loaded-spine optimum: t_s = u_s(x) for s < k and
t_k = N, where u_0 = 0, u_1 = x, u_{s+1} = u_s^2 / (u_{s-1} + 1/2) and x is the unique positive root of
u_k(x) = N (loaded-spinal chain bound of Section 4.2, part (i)).  Every quantity is an element of the rational function field Q(x).

  1. Balance: the force and moment residual of every block is the zero element of Q(x) (exact arithmetic in
     the rational function field of SymPy), so balance holds for every x, in particular at the root; the
     overhang equals the sum of the steps, which is S_k(N).
  2. Inequalities: every strict or weak inequality g(x) >= 0 (positive forces, positive added loads,
     contact points inside the overlaps, positive overlap of every loaded contact, unit same-layer spacing,
     single root, unique principal block, spine blocks resting only on the next spine block) is decided on
     an isolating interval [L, R] of the root: either a rigorous interval evaluation (mpmath.iv, outward
     rounding) of the numerator and the denominator of g over [L, R] excludes 0, or, failing that, [L, R] is
     refined until Sturm counting shows that neither has a real root in it and g(L) gives the sign; g
     identically 0 is accepted for weak inequalities (contact points at block edges).
  3. Height: with rigorous interval arithmetic (mpmath.iv), the loaded-spine value at height k exceeds the
     values at all other heights, and all added loads at height k are positive, so S*(N) = S*_k(N).
For N = 5 the construction is eq. (66) with t_2 = sqrt(15/2); for N = 7 it is the height-3 construction
with loads (0, x, 2x^2, 7), 8x^4 - 14x - 7 = 0.  The upper bounds that make these the spinal optima at
N = 5, 7 (other heights and the faces with a zero added load) are checked by spine_faces.py.
"""
import argparse
import sys
from fractions import Fraction as Fr

import sympy as sp
from mpmath import iv, mp
from sympy import QQ
from sympy.polys.fields import field

from verify_rational_stacks import SPECS_4_9_SPINAL, SPECS_10_12, SPECS_13_16

X = sp.Symbol("x")
F, XF = field("x", QQ)                     # the rational function field Q(x); XF is its generator
HALF = F(Fr(1, 2))
EXTRA = {  # N = 17..19 (seeds of Appendix C.2)
    17: (6, ["-2/5", "-4/5", "-6/5", "-13/10", "-8/5"], "3/5"),
    18: (6, ["-3/10", "-4/5", "-11/10", "-7/5", "-8/5"], "3/5"),
    19: (7, ["-9/40", "-29/40", "-43/40", "-27/20", "-3/2", "-13/8"], "13/20"),
}


def R(s):
    return F(Fr(str(s)))


def loads_symbolic(N, k):
    u = [F.zero, XF]
    for _ in range(1, k):
        u.append(u[-1] ** 2 / (u[-2] + HALF))
    return u[:k] + [F(N)], u[k] - N


def shielded(N, k, t, ys, z1):
    """The construction of verify_rational_stacks.reconstruct_spinal, over Q(x)."""
    w = [None] + [t[r] - t[r - 1] - 1 for r in range(1, k + 1)]
    d = [None] + [(w[r] + HALF) / t[r] for r in range(1, k + 1)]
    x = [None] * (k + 1)
    x[k] = d[k] - 1
    for r in range(k - 1, 0, -1):
        x[r] = x[r + 1] + d[r]
    y = [None] + [R(s) for s in ys]
    u, z = [None] * (k + 1), [None] * (k + 1)
    u[k] = z[k] = F.zero
    for r in range(k - 1, 1, -1):
        u[r] = u[r + 1] + w[r + 1] - 1
        z[r] = (u[r + 1] * z[r + 1] + w[r + 1] * x[r + 1] - (y[r] + HALF)) / u[r]
    v1 = N - 2 * k
    u[1] = 1 - w[1]
    m1 = u[2] * z[2] + w[2] * x[2] - (y[1] + HALF)
    if v1 == 0:
        z[1], pv = m1 / u[1], None
    else:
        z[1] = R(z1)
        pv = (m1 - u[1] * z[1]) / v1
    bridge = u[1] * z[1] + w[1] * x[1]          # centre of the bridge block
    C, L = {}, {}                               # centres, layers
    for r in range(1, k + 1):
        C[f"B{r}"], L[f"B{r}"] = x[r] + HALF, k + 1 - r
    for r in range(1, k):
        C[f"H{r}"], L[f"H{r}"] = y[r] + HALF, k + 1 - r
    C["G"], L["G"] = bridge, k + 1
    for s in range(1, v1 + 1):
        C[f"T{s}"], L[f"T{s}"] = pv, k + s
    K = []                                      # (upper, lower, force, point)
    K.append((f"B{k}", "table", F(N), F.zero))
    for r in range(1, k):
        K.append((f"B{r}", f"B{r + 1}", t[r], x[r + 1] + 1))
        K.append((f"H{r}", f"B{r + 1}", w[r + 1], x[r + 1]))
    K.append(("G", "B1", w[1], x[1]))
    K.append(("G", "H1", u[1], z[1]))
    for r in range(2, k):
        K.append((f"H{r - 1}", f"H{r}", u[r], z[r]))
    for s in range(1, v1 + 1):
        K.append((f"T{s}", "H1" if s == 1 else f"T{s - 1}", F(v1 - s + 1), pv))
    return C, L, K, w, d


def bridgeless_five(t2):
    """Eq. (66) with X(C_2) = -18/25 (centre), the N = 5 exception."""
    d = {1: HALF, 2: 1 - F(Fr(3, 2)) / t2, 3: (F(Fr(9, 2)) - t2) / 5}
    x3 = d[3] - 1
    x2, x1 = x3 + d[2], x3 + d[2] + d[1]
    w2, w3, u2 = t2 - 2, 4 - t2, 3 - t2
    c2 = F(Fr(-18, 25))
    m = w3 * x3 - c2
    c1 = w2 * x2 + m
    C = {"B3": x3 + HALF, "B2": x2 + HALF, "B1": x1 + HALF, "C2": c2, "C1": c1}
    L = {"B3": 1, "B2": 2, "B1": 3, "C2": 2, "C1": 3}
    K = [("B3", "table", F(5), F.zero), ("B2", "B3", t2, x3 + 1), ("C2", "B3", w3, x3),
         ("B1", "B2", F.one, x2 + 1), ("C1", "B2", w2, x2), ("C1", "C2", u2, m / u2)]
    w = [None, F.zero, t2 - 2, 4 - t2]
    return C, L, K, w, d


def horner(p, xint):
    """Rigorous interval enclosure of the univariate polynomial p (a SymPy PolyElement over QQ) on xint."""
    r = iv.mpf(0)
    for c in p.to_dense():
        c = Fr(int(c.numerator), int(c.denominator))
        r = r * xint + iv.mpf(c.numerator) / c.denominator
    return r


class Root:
    """An algebraic number given by a squarefree factor P of its polynomial and an isolating interval."""

    def __init__(self, P, lo, hi):
        self.P = sp.Poly(P, X)
        q = lambda v: sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else sp.Rational(v)
        self.lo, self.hi = q(lo), q(hi)
        assert self.P.count_roots(self.lo, self.hi) == 1 and self.P.eval(self.lo) != 0 and self.P.eval(self.hi) != 0

    def sign(self, g):
        """Sign of the element g of Q(x) at the root (0 only if g is identically zero)."""
        g = F(g)
        if g == F.zero:
            return 0
        xint = iv.mpf([iv.mpf(int(self.lo.p)) / int(self.lo.q), iv.mpf(int(self.hi.p)) / int(self.hi.q)])
        en, ed = horner(g.numer, xint), horner(g.denom, xint)
        if (en.a > 0 or en.b < 0) and (ed.a > 0 or ed.b < 0):          # fast path: both signs fixed on [L, R]
            return (1 if en.a > 0 else -1) * (1 if ed.a > 0 else -1)
        num, den = g.numer.as_expr(), g.denom.as_expr()
        for q in (num, den):
            q = sp.Poly(q, X)
            if q.degree() <= 0:
                continue
            while q.count_roots(self.lo, self.hi) > 0:   # refine until the root interval avoids q's roots
                if sp.gcd(q, self.P).degree() > 0:
                    raise SystemExit("inequality degenerates at the root (common factor)")
                self.refine()
        v = sp.Poly(num, X).eval(self.lo) * sp.Poly(den, X).eval(self.lo)
        assert v != 0
        return 1 if v > 0 else -1

    def refine(self):
        mid = (self.lo + self.hi) / 2
        if self.P.eval(mid) == 0:
            raise SystemExit("rational root")
        if self.P.count_roots(self.lo, mid) == 1:
            self.hi = mid
        else:
            self.lo = mid


def check(N, C, L, K, w, d, root, kk, strict_w):
    out = {"min_force": None}
    # 1. balance (identities in Q(x))
    f = {b: -F.one for b in C}
    m = {b: -C[b] for b in C}
    for (a, b, Fc, P) in K:
        f[a] += Fc
        m[a] += Fc * P
        if b != "table":
            f[b] -= Fc
            m[b] -= Fc * P
    for b in C:
        assert f[b] == F.zero and m[b] == F.zero, f"balance of {b}"
    assert sum(C.values(), F.zero) == F.zero, "centroid identity"
    assert C["B1"] + HALF - sum((d[r] for r in range(1, kk + 1)), F.zero) == F.zero, "overhang identity"
    # 2. inequalities
    pos = lambda g, what: (root.sign(g) > 0) or sys.exit(f"N={N}: {what} fails")
    weak = lambda g, what: (root.sign(g) >= 0) or sys.exit(f"N={N}: {what} fails")
    for r in range(1, kk + 1):
        (pos if strict_w(r) else weak)(w[r], f"added load {r}")
    assert [b for b in L if L[b] == 1] == [f"B{kk}"], "single root"
    for (a, b, Fc, P) in K:
        pos(Fc, f"force {a}/{b}")
        if b == "table":
            weak(P - (C[a] - HALF), "table contact left"); weak(C[a] + HALF - P, "table contact right")
            pos(-(C[a] - HALF), "table overlap")
            continue
        assert L[a] == L[b] + 1, f"{a} on {b}: layers"
        for q in (a, b):
            weak(P - (C[q] - HALF), f"point of {a}/{b} in {q} (left)")
            weak(C[q] + HALF - P, f"point of {a}/{b} in {q} (right)")
        s = root.sign(C[a] - C[b])
        pos(1 - s * (C[a] - C[b]) if s else F.one, f"positive overlap {a}/{b}")
    for l in set(L.values()):
        bs = [b for b in L if L[b] == l]
        for i in range(len(bs)):
            for j in range(i + 1, len(bs)):
                s = root.sign(C[bs[i]] - C[bs[j]])
                assert s != 0
                weak(s * (C[bs[i]] - C[bs[j]]) - 1, f"spacing {bs[i]}/{bs[j]}")
    for b in C:
        if b != "B1":
            pos(C["B1"] - C[b], f"principal B1 vs {b}")
    for r in range(1, kk):
        for b in C:
            if L[b] == L[f"B{r}"] - 1 and b != f"B{r + 1}":
                s = root.sign(C[f"B{r}"] - C[b])
                pos(s * (C[f"B{r}"] - C[b]) - 1, f"B{r} does not rest on {b}")
    return True


def iv_shat(K, N, lo, hi):
    """Enclosure of the unconstrained loaded-spine value at height K for x in [lo, hi] (mpmath.iv)."""
    x = iv.mpf([lo, hi])
    u = [iv.mpf(0), x]
    for _ in range(1, K):
        u.append(u[-1] ** 2 / (u[-2] + iv.mpf(0.5)))
    terminal = u[K]
    u[K] = iv.mpf(N)
    val = sum(1 - (u[s - 1] + iv.mpf(0.5)) / u[s] for s in range(1, K + 1))
    return terminal, val, u


def bracket(K, N):
    """Rational bracket [lo, hi] of the root of u_K(x) = N, with rigorous end signs."""
    if K == 1:
        return None
    lo, hi = Fr(1, 2), Fr(N)
    for _ in range(110):
        mid = (lo + hi) / 2
        t, _, _ = iv_shat(K, N, mp.mpf(mid.numerator) / mid.denominator, mp.mpf(mid.numerator) / mid.denominator)
        if t.b < N:
            lo = mid
        elif t.a > N:
            hi = mid
        else:
            break
    tl, _, _ = iv_shat(K, N, mp.mpf(lo.numerator) / lo.denominator, mp.mpf(lo.numerator) / lo.denominator)
    th, _, _ = iv_shat(K, N, mp.mpf(hi.numerator) / hi.denominator, mp.mpf(hi.numerator) / hi.denominator)
    assert tl.b < N < th.a
    return lo, hi


def heights(N):
    vals = {1: iv.mpf(1) - iv.mpf(0.5) / N}
    br = {}
    for K in range(2, N + 1):
        lo, hi = bracket(K, N)
        br[K] = (lo, hi)
        _, v, u = iv_shat(K, N, mp.mpf(lo.numerator) / lo.denominator, mp.mpf(hi.numerator) / hi.denominator)
        vals[K] = (v, u)
    best = max(range(2, N + 1), key=lambda K: vals[K][0].a)
    others = [vals[K][0].b if K > 1 else vals[1].b for K in range(1, N + 1) if K != best]
    assert vals[best][0].a > max(others), f"N={N}: height not separated"
    incs_pos = all((vals[best][1][s] - vals[best][1][s - 1] - 1).a > 0 for s in range(1, best + 1))
    return best, br[best], float(vals[best][0].a - max(others)), incs_pos, vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=19)
    a = ap.parse_args()
    mp.dps = iv.dps = 60
    specs = {n: (k, ys, z) for (n, k, _, ys, z) in SPECS_4_9_SPINAL + SPECS_10_12 + SPECS_13_16}
    specs.update(EXTRA)
    for N in range(a.nmin, a.nmax + 1):
        if N == 5:
            P = 2 * X ** 2 - 15
            root = Root(P, Fr(27386, 10000), Fr(27387, 10000))
            C, L, K, w, d = bridgeless_five(XF)
            check(N, C, L, K, w, d, root, 3, lambda r: r != 1)
            val = (12 - sp.sqrt(30)) / 5
            assert sp.simplify((C["B1"] + HALF).as_expr().subs(X, sp.sqrt(sp.Rational(15, 2))) - val) == 0
            print(f"N= 5: spinal stack of overhang (12 - sqrt(30))/5 = {float(val):.12f}: exact balance, "
                  f"all inequalities verified (upper bound: spine_faces.py)")
            continue
        k, ys, z = specs[N]
        t, Pk = loads_symbolic(N, k)
        num = sp.Poly(Pk.numer.as_expr(), X)
        best, (lo, hi), gap, incs_pos, _ = heights(N)
        if N != 7:
            assert best == k and incs_pos, f"N={N}: construction height {k}, loaded optimum at {best}"
        lo2, hi2 = bracket(k, N)
        fac = [f for f, _ in sp.factor_list(num)[1] if f.count_roots(sp.Rational(lo2.numerator, lo2.denominator),
                                                                     sp.Rational(hi2.numerator, hi2.denominator)) == 1][0]
        root = Root(fac.as_expr(), lo2, hi2)
        C, L, K, w, d = shielded(N, k, t, ys, z)
        check(N, C, L, K, w, d, root, k, lambda r: True)
        xr = iv.mpf([iv.mpf(int(root.lo.p)) / int(root.lo.q), iv.mpf(int(root.hi.p)) / int(root.hi.q)])
        S = horner((C["B1"] + HALF).numer, xr) / horner((C["B1"] + HALF).denom, xr)
        S = mp.mpf(S.a)
        if N == 7:
            print(f"N= 7: spinal stack at height 3 of overhang {float(S):.12f} (8x^4 - 14x - 7 = 0): exact "
                  f"balance, all inequalities verified (upper bound: spine_faces.py)")
        else:
            print(f"N={N:2d}: S(N) = S*(N) = S*_{k}(N) ~ {float(S):.12f}; t_1 of degree {fac.degree()}; exact "
                  f"balance, all inequalities verified; height gap {gap:.5f}")
    print(f"ALL ALGEBRAIC SPINAL STACKS VERIFIED (N = {a.nmin}..{a.nmax})")


if __name__ == "__main__":
    main()
