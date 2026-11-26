"""Turn a floating-point witness stack into an exactly balanced rational stack.

Usage:  python3 src/exact_witness.py WITNESS.json [INDEX] [--out EXACT.json]   (construct; needs scipy)
        python3 src/exact_witness.py --check EXACT.json                     (verify only; standard library)

The witness (a *_witnesses.json file of ns_check.py, or a stack with "contacts") gives the layers, the contact
pairs and the contact forces.  The forces are rounded to the grid 1e-9 and made to satisfy vertical force
balance exactly (from the top layer down, the last support of each block takes the remainder).  With the
forces fixed, every remaining condition is linear in the centers X and the contact points P: moment balance
of each block, contact points inside the closed contact intervals, unit same-layer spacing in the witness
order, and the table contact of block 1 at or left of the table edge 0.  This linear program is solved with
the principal block's overhang as objective; the active constraints of the floating-point optimum are then
solved again in exact rational arithmetic, and the result is verified with fractions only: zero residual in
force and moment balance, positive forces, contact points inside the closed contact intervals, same-layer
spacing at least 1, and single root.  Units: block length 1, table edge at 0.
"""
import json
import sys
from fractions import Fraction as Fr

GRID = 10 ** 9


def load(path, idx):
    d = json.load(open(path))
    return d["witnesses"][idx] if "witnesses" in d else d


def exact_forces(stack):
    L = {r["block"]: r["layer"] for r in stack}
    sup = {r["block"]: [(c["on"], c["F"]) for c in r["contacts"]] for r in stack}
    F = {}
    for b in sorted(L, key=lambda b: -L[b]):
        total = 1 + sum(F[(a, bb)] for (a, bb) in F if bb == b)
        rest = total
        for k, (j, f) in enumerate(sup[b]):
            if k < len(sup[b]) - 1:
                F[(b, j)] = Fr(round(f * GRID), GRID)
                rest -= F[(b, j)]
            else:
                F[(b, j)] = rest
    assert all(f > 0 for f in F.values()), "non-positive force after rounding"
    return L, F


def build(L, F, stack):
    blocks = sorted(L)
    contacts = sorted(F)
    nx, npv = len(blocks), len(contacts)
    xi = {b: k for k, b in enumerate(blocks)}
    pi = {c: nx + k for k, c in enumerate(contacts)}
    nv = nx + npv
    A_eq, b_eq, A_ub, b_ub = [], [], [], []

    def row(coefs):
        r = [Fr(0)] * nv
        for k, v in coefs:
            r[k] += v
        return r
    for b in blocks:   # sum F P (below) - sum F P (above) - X_b = 0
        co = [(pi[c], F[c]) for c in contacts if c[0] == b] + [(pi[c], -F[c]) for c in contacts if c[1] == b]
        A_eq.append(row(co + [(xi[b], Fr(-1))]))
        b_eq.append(Fr(0))
    half = Fr(1, 2)
    for (i, j) in contacts:   # P within [X_i - 1/2, X_i + 1/2] and the same for j; table edge at 0
        p = pi[(i, j)]
        for q in ([i] if j == 0 else [i, j]):
            A_ub.append(row([(xi[q], Fr(1)), (p, Fr(-1))])); b_ub.append(half)    # X_q - P <= 1/2
            A_ub.append(row([(p, Fr(1)), (xi[q], Fr(-1))])); b_ub.append(half)    # P - X_q <= 1/2
        if j == 0:
            A_ub.append(row([(p, Fr(1))])); b_ub.append(Fr(0))
    X0 = {r["block"]: r["center_x"] for r in stack}
    for l in set(L.values()):   # same-layer spacing, in the order of the witness
        bs = sorted((b for b in blocks if L[b] == l), key=lambda b: X0[b])
        for a, c in zip(bs[:-1], bs[1:]):
            A_ub.append(row([(xi[a], Fr(1)), (xi[c], Fr(-1))])); b_ub.append(Fr(-1))
    top = max(blocks, key=lambda b: X0[b])   # principal block of the witness
    obj = [Fr(0)] * nv
    obj[xi[top]] = Fr(-1)
    return blocks, contacts, xi, pi, A_eq, b_eq, A_ub, b_ub, obj, top


def solve_exact(A_eq, b_eq, A_ub, b_ub, obj):
    import numpy as np
    from scipy.optimize import linprog
    f = lambda M: np.array([[float(v) for v in r] for r in M])
    res = linprog(np.array([float(v) for v in obj]), A_ub=f(A_ub), b_ub=[float(v) for v in b_ub],
                  A_eq=f(A_eq), b_eq=[float(v) for v in b_eq], bounds=[(None, None)] * len(obj), method="highs")
    if res.status != 0:
        raise SystemExit(f"linear program failed: {res.message}")
    x = res.x
    active = [k for k, (r, rhs) in enumerate(zip(A_ub, b_ub))
              if abs(sum(float(a) * v for a, v in zip(r, x)) - float(rhs)) <= 1e-9]
    rows = [list(r) for r in A_eq] + [list(A_ub[k]) for k in active]
    rhs = list(b_eq) + [b_ub[k] for k in active]
    return exact_solve(rows, rhs, len(obj))


def exact_solve(rows, rhs, n):
    """Exact solution of a consistent linear system with a unique solution (Gaussian elimination)."""
    M = [r[:] + [v] for r, v in zip(rows, rhs)]
    piv, r0 = [], 0
    for c in range(n):
        k = next((k for k in range(r0, len(M)) if M[k][c] != 0), None)
        if k is None:
            continue
        M[r0], M[k] = M[k], M[r0]
        M[r0] = [v / M[r0][c] for v in M[r0]]
        for k2 in range(len(M)):
            if k2 != r0 and M[k2][c] != 0:
                fac = M[k2][c]
                M[k2] = [a - fac * b for a, b in zip(M[k2], M[r0])]
        piv.append(c)
        r0 += 1
    if len(piv) < n:
        raise SystemExit("active constraints do not determine a unique point")
    assert all(all(v == 0 for v in r) for r in M[r0:]), "inconsistent active set"
    sol = [Fr(0)] * n
    for k, c in enumerate(piv):
        sol[c] = M[k][n]
    return sol


def verify(L, F, X, P):
    half = Fr(1, 2)
    for b in L:
        below = [c for c in F if c[0] == b]
        above = [c for c in F if c[1] == b]
        assert sum(F[c] for c in below) - sum(F[c] for c in above) == 1, f"force balance, block {b}"
        assert sum(F[c] * P[c] for c in below) - sum(F[c] * P[c] for c in above) == X[b], f"moment, block {b}"
    for (i, j), f in F.items():
        assert f > 0
        if j == 0:
            assert X[i] - half <= P[(i, j)] <= min(X[i] + half, 0), "table contact"
        else:
            assert max(X[i], X[j]) - half <= P[(i, j)] <= min(X[i], X[j]) + half, f"contact {i} on {j}"
    for l in set(L.values()):
        bs = sorted((b for b in L if L[b] == l), key=lambda b: X[b])
        assert all(X[c] - X[a] >= 1 for a, c in zip(bs[:-1], bs[1:])), f"spacing, layer {l}"
    assert [b for b in L if L[b] == 1] == [1], "single root"


def spinal(L, X, strict):
    """Support set of the principal block under the closed (|dX| <= 1) or strict (|dX| < 1) resting relation."""
    top = max(X, key=lambda b: X[b])
    rests = lambda i, j: L[i] == L[j] + 1 and (abs(X[i] - X[j]) < 1 if strict else abs(X[i] - X[j]) <= 1)
    S, todo = {top}, [top]
    while todo:
        i = todo.pop()
        for j in X:
            if j not in S and rests(i, j):
                S.add(j)
                todo.append(j)
    per = {}
    for b in S:
        per[L[b]] = per.get(L[b], 0) + 1
    return all(v == 1 for v in per.values())


def check(path):
    d = json.load(open(path))
    X = {b["block"]: Fr(b["center_x"]) for b in d["blocks"]}
    L = {b["block"]: b["layer"] for b in d["blocks"]}
    F = {(c["block"], c["on"]): Fr(c["F"]) for c in d["contacts"]}
    P = {(c["block"], c["on"]): Fr(c["P"]) for c in d["contacts"]}
    verify(L, F, X, P)
    z = max(X.values()) + Fr(1, 2)
    assert str(z) == d["overhang"], "stored overhang differs"
    ov = min(1 - abs(X[i] - X[j]) for (i, j) in F if j)
    print(f"{path}: exact balance verified; overhang {float(z):.12f}; smallest loaded overlap {float(ov):.3e}; "
          f"spinal (closed / Paterson-Zwick): {spinal(L, X, False)} / {spinal(L, X, True)}")


def main():
    if "--check" in sys.argv:
        check(sys.argv[sys.argv.index("--check") + 1])
        return
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    if out in args:
        args.remove(out)
    w = load(args[0], int(args[1]) if len(args) > 1 else 0)
    stack = w["stack"]
    L, F = exact_forces(stack)
    blocks, contacts, xi, pi, A_eq, b_eq, A_ub, b_ub, obj, top = build(L, F, stack)
    sol = solve_exact(A_eq, b_eq, A_ub, b_ub, obj)
    X = {b: sol[xi[b]] for b in blocks}
    P = {c: sol[pi[c]] for c in contacts}
    try:
        for r, rhs in zip(A_ub, b_ub):
            assert sum(a * v for a, v in zip(r, sol)) <= rhs, "inequality violated in exact arithmetic"
        verify(L, F, X, P)
    except AssertionError as e:
        raise SystemExit(f"no exactly balanced stack obtained with the forces of this witness ({e})")
    z = max(X.values()) + Fr(1, 2)
    print(f"exactly balanced rational stack: overhang {z} = {float(z):.12f}")
    for (i, j) in contacts:
        if j and abs(X[i] - X[j]) == 1:
            print(f"  point contact: block {i} on block {j}, force {float(F[(i, j)]):.6f}")
    if out:
        rec = {"overhang": str(z), "overhang_float": float(z),
               "blocks": [{"block": b, "layer": L[b], "center_x": str(X[b])} for b in blocks],
               "contacts": [{"block": c[0], "on": c[1], "F": str(F[c]), "P": str(P[c])} for c in contacts]}
        json.dump(rec, open(out, "w"), indent=1)
        print(f"written {out}")


if __name__ == "__main__":
    main()
