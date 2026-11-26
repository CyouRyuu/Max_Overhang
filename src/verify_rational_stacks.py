#!/usr/bin/env python3
"""Exact rational verifier and closed-form constructor for single-root overhang stacks.
Zero floating-point arithmetic is used in construction or verification (eta = 0, W = 1/2).
"""
import argparse
import json
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 30


def verify_rational_stack(N, layers, X, contacts):
    """Independent exact-rational verifier for an N-block single-root stack (eta=0, W=1/2).

    Parameters
    ----------
    N : int
        Number of blocks (indexed 1..N).
    layers : dict[int, int]
        Mapping block_id -> layer_id (1-indexed, layer 1 is the base).
    X : dict[int, Fraction]
        Mapping block_id -> exact rational horizontal center X_i.
    contacts : dict[tuple[int, int], tuple[Fraction, Fraction]]
        Mapping (i, j) -> (F_{ij}, tau_{ij}) with tau_{ij} = F_{ij} * P_{ij}.
        Base-to-table contact is keyed by (1, 0).

    Returns
    -------
    tuple[Fraction, tuple[int, ...]]
        (exact_overhang, layer_profile)
    """
    W = Fraction(1, 2)
    assert len(X) == N and len(layers) == N, "Block count mismatch."
    H = max(layers.values())
    by_layer = {l: [] for l in range(1, H + 1)}
    for i in range(1, N + 1):
        assert isinstance(X[i], Fraction), f"X[{i}] is not a Fraction."
        assert 1 <= layers[i] <= H, f"Invalid layer for block {i}."
        by_layer[layers[i]].append(i)

    # 1. Single-root base and canonical ordering + same-layer non-overlap (>= 1)
    assert by_layer[1] == [1], f"Layer 1 must be [1], got {by_layer[1]}."
    profile = tuple(len(by_layer[l]) for l in range(1, H + 1))
    expected_id = 1
    for l in range(1, H + 1):
        bs = by_layer[l]
        assert bs == list(range(expected_id, expected_id + len(bs))), (
            f"Non-canonical block IDs in layer {l}: {bs}."
        )
        expected_id += len(bs)
        for a, b in zip(bs[:-1], bs[1:]):
            assert X[a] - X[b] >= 1, f"Same-layer overlap in L{l}: X[{a}] - X[{b}] = {X[a] - X[b]}."

    # 2. Contact intervals, table clearance, adjacency, and support existence
    assert (1, 0) in contacts, "Missing table contact (1, 0)."
    has_lower = {i: False for i in range(1, N + 1)}
    for (i, j), (F, tau) in contacts.items():
        assert isinstance(F, Fraction) and isinstance(tau, Fraction), "Forces/torques must be Fraction."
        assert F >= 0, f"Negative contact force F[{i},{j}] = {F}."
        if F == 0:
            assert tau == 0, f"Nonzero torque for zero force on ({i},{j})."
            continue
        P = tau / F
        if j == 0:
            assert i == 1, f"Only block 1 may touch the table, got ({i}, 0)."
            assert P <= 0, f"Table contact P_{{1,0}} = {P} > 0."
            assert X[1] - W <= P <= X[1] + W, f"P_{{1,0}} = {P} outside base block."
            has_lower[1] = True
        else:
            assert 1 <= j < i <= N, f"Invalid contact indices ({i}, {j})."
            assert layers[i] == layers[j] + 1, f"Non-adjacent contact ({i}, {j})."
            lo = max(X[i] - W, X[j] - W)
            hi = min(X[i] + W, X[j] + W)
            assert lo <= P <= hi, f"P_{i},{j} = {P} outside [{lo}, {hi}]."
            has_lower[i] = True
    assert all(has_lower.values()), f"Unsupported block detected: {has_lower}."

    # 3. Exact per-block vertical force balance (= 1) and moment balance (= X_i)
    for i in range(1, N + 1):
        net_F = (
            sum(F for (u, d), (F, _) in contacts.items() if u == i)
            - sum(F for (u, d), (F, _) in contacts.items() if d == i)
        )
        net_tau = (
            sum(tau for (u, d), (_, tau) in contacts.items() if u == i)
            - sum(tau for (u, d), (_, tau) in contacts.items() if d == i)
        )
        assert net_F == 1, f"Force imbalance on block {i}: {net_F} != 1."
        assert net_tau == X[i], f"Moment imbalance on block {i}: {net_tau} != {X[i]}."

    overhang = max(X[i] + W for i in range(1, N + 1))
    return overhang, profile


def reconstruct_spinal(N, k, t_list, y_shields, z1_val):
    """Pure-Fraction closed-form constructor for single-tower shielded spinal stacks."""
    W = Fraction(1, 2)
    t = [Fraction(0)] + [Fraction(s) for s in t_list] + [Fraction(N)]
    w = [None] + [t[r] - t[r - 1] - 1 for r in range(1, k + 1)]
    d = [None] + [(w[r] + W) / t[r] for r in range(1, k + 1)]

    x = [None] * (k + 1)
    x[k] = d[k] - 1
    for r in range(k - 1, 0, -1):
        x[r] = x[r + 1] + d[r]

    y = [None] * (k + 1)
    for r in range(1, k):
        y[r] = Fraction(y_shields[r - 1])

    u = [None] * (k + 1)
    z = [None] * (k + 1)
    u[k] = Fraction(0)
    z[k] = Fraction(0)
    for r in range(k - 1, 1, -1):
        u[r] = u[r + 1] + w[r + 1] - 1
        z[r] = (u[r + 1] * z[r + 1] + w[r + 1] * x[r + 1] - (y[r] + W)) / u[r]

    v1 = N - 2 * k
    u[1] = Fraction(1) - w[1]
    M1_rhs = u[2] * z[2] + w[2] * x[2] - (y[1] + W)
    if v1 == 0:
        z[1] = M1_rhs / u[1]
        p_v1 = None
    else:
        z[1] = Fraction(z1_val)
        p_v1 = (M1_rhs - u[1] * z[1]) / v1
    y[0] = u[1] * z[1] + w[1] * x[1] - W

    raw = []
    for r in range(k, 0, -1):
        raw.append((k + 1 - r, x[r] + W, ("spine", r)))
    for r in range(k - 1, 0, -1):
        raw.append((k + 1 - r, y[r] + W, ("shield", r)))
    raw.append((k + 1, y[0] + W, ("bridge", 0)))
    for s in range(1, v1 + 1):
        raw.append((k + s, p_v1, ("tower1", s)))
    raw.sort(key=lambda item: (item[0], -item[1]))

    tag_to_id, layers, X, roles = {}, {}, {}, {}
    for idx, (lay, cx, tag) in enumerate(raw, start=1):
        tag_to_id[tag] = idx
        layers[idx] = lay
        X[idx] = cx
        roles[idx] = f"{tag[0]}_{tag[1]}"

    contacts = {}

    def add_c(ut, dt, F, P):
        ui = tag_to_id[ut]
        di = 0 if dt is None else tag_to_id[dt]
        contacts[(ui, di)] = (F, F * P)

    add_c(("spine", k), None, Fraction(N), Fraction(0))
    for r in range(1, k):
        add_c(("spine", r), ("spine", r + 1), t[r], x[r + 1] + 1)
    add_c(("bridge", 0), ("spine", 1), w[1], x[1])
    add_c(("bridge", 0), ("shield", 1), u[1], z[1])
    for r in range(1, k):
        add_c(("shield", r), ("spine", r + 1), w[r + 1], x[r + 1])
    for r in range(2, k):
        add_c(("shield", r - 1), ("shield", r), u[r], z[r])
    if v1 > 0:
        add_c(("tower1", 1), ("shield", 1), Fraction(v1), p_v1)
        for s in range(2, v1 + 1):
            add_c(("tower1", s), ("tower1", s - 1), Fraction(v1 - s + 1), p_v1)

    oh, prof = verify_rational_stack(N, layers, X, contacts)
    return oh, prof, layers, X, contacts, roles


SPECS_13_16 = [
    (13, 5, ["1263/1000", "319/100", "2309/400", "903/100"], ["-3/5", "-31/50", "-11/10", "-29/20"], "2/5"),
    (14, 5, ["1277/1000", "1631/500", "1497/250", "19061/2000"], ["-1/2", "-3/5", "-21/20", "-7/5"], "1/2"),
    (15, 6, ["1239/1000", "6141/2000", "2711/500", "3293/400", "2289/200"], ["-1/2", "-101/200", "-19/20", "-13/10", "-8/5"], "1/2"),
    (16, 6, ["156/125", "779/250", "2777/500", "17061/2000", "24041/2000"], ["-9/20", "-1/2", "-9/10", "-5/4", "-31/20"], "11/20"),
]


SPECS_4_9_SPINAL = [
    (4, 2, ["707/500"], ["-1"], None),
    (6, 3, ["639/500", "1633/500"], ["-1", "-6/5"], None),
    (7, 3, ["2679/2000", "7177/2000"], ["-3/4", "-6/5"], "1/5"),
    (8, 4, ["154/125", "6069/2000", "2127/400"], ["-17/20", "-1", "-7/5"], None),
    (9, 4, ["631/500", "6371/2000", "11519/2000"], ["-13/20", "-1", "-7/5"], "11/40"),
]


def reconstruct_n5_bridgeless(t2_str="5477/2000", X3_str="-18/25"):
    """Pure-Fraction closed-form constructor for the N=5 Paterson-Zwick exception (profile (1,2,2))."""
    N = 5
    W = Fraction(1, 2)
    t1 = Fraction(1)
    t2 = Fraction(t2_str)
    t3 = Fraction(5)
    w2 = t2 - t1 - 1
    w3 = t3 - t2 - 1
    d1 = W
    d2 = (w2 + W) / t2
    d3 = (w3 + W) / t3

    x3 = d3 - 1
    x2 = x3 + d2
    x1 = x2 + d1

    u2 = Fraction(1) - w2
    X1 = x3 + W
    X2 = x2 + W
    X3 = Fraction(X3_str)
    X4 = x1 + W
    tau_53 = w3 * x3 - X3
    X5 = w2 * x2 + tau_53

    layers = {1: 1, 2: 2, 3: 2, 4: 3, 5: 3}
    X = {1: X1, 2: X2, 3: X3, 4: X4, 5: X5}
    roles = {1: "spine_3", 2: "spine_2", 3: "cw_2", 4: "spine_1", 5: "cw_1"}
    contacts = {
        (1, 0): (Fraction(5), Fraction(0)),
        (2, 1): (t2, t2 * (x3 + 1)),
        (3, 1): (w3, w3 * x3),
        (4, 2): (Fraction(1), x2 + 1),
        (5, 2): (w2, w2 * x2),
        (5, 3): (u2, tau_53),
    }
    oh, prof = verify_rational_stack(N, layers, X, contacts)
    return oh, prof, layers, X, contacts, roles


SPECS_10_12 = [
    (10, 4, ["129/100", "416/125", "6187/1000"], ["-11/20", "-1", "-7/5"], "2/5"),
    (11, 5, ["154/125", "607/200", "5317/1000", "8"], ["-11/20", "-9/10", "-6/5", "-303/200"], "2/5"),
    (12, 5, ["156/125", "623/200", "111/20", "213/25"], ["-1/2", "-4/5", "-6/5", "-3/2"], "17/40"),
]


def verify_all_n4_to_n16():
    """Run exact Fraction zero-residual verification across all N = 4..16 (including N = 5, 7 exceptions)."""
    all_spinal = sorted(SPECS_4_9_SPINAL + SPECS_10_12 + SPECS_13_16, key=lambda item: item[0])
    results = {}
    for N, k, tl, ys, z1 in all_spinal:
        oh, prof, layers, X, contacts, roles = reconstruct_spinal(N, k, tl, ys, z1)
        results[N] = (oh, prof, layers, X, contacts, roles)
    results[5] = reconstruct_n5_bridgeless()
    for N in range(4, 17):
        oh, prof = results[N][0], results[N][1]
        oh_dec = Decimal(oh.numerator) / Decimal(oh.denominator)
        print(f"[PASS] N={N:2d} | profile={prof} | overhang={oh} ({oh_dec:.12f})")
    return results



def main():
    parser = argparse.ArgumentParser(description="Exact rational verification of the feasible stacks (lower-bound witnesses), N = 4..16.")
    parser.add_argument("--print-all", action="store_true", help="Print full exact fractions for every block and contact.")
    parser.add_argument("--dump-json", type=str, default=None, help="Export exact verified stacks to JSON.")
    args = parser.parse_args()

    specs = {item[0]: item for item in SPECS_4_9_SPINAL + SPECS_10_12 + SPECS_13_16}
    exported = {}
    for N in range(4, 17):
        if N == 5:
            oh, prof, layers, X, contacts, roles = reconstruct_n5_bridgeless()
            k = None
        else:
            _, k, tl, ys, z1 = specs[N]
            oh, prof, layers, X, contacts, roles = reconstruct_spinal(N, k, tl, ys, z1)
        oh_dec = Decimal(oh.numerator) / Decimal(oh.denominator)
        ktxt = "bridgeless" if k is None else f"k={k}"
        print(f"[PASS] N={N:2d} | {ktxt:10s} | profile={prof} | overhang={oh} ({oh_dec:.12f})")
        if args.print_all:
            for i in range(1, N + 1):
                print(f"  Block {i:2d} (L{layers[i]}, {roles[i]:8s}): X = {X[i]} ({float(X[i]):+.10f})")
            for (i, j) in sorted(contacts):
                F, tau = contacts[(i, j)]
                print(f"  Contact ({i:2d},{j:2d}): F = {F}, tau = {tau}, P = {float(tau / F):+.10f}")
        exported[str(N)] = {
            "N": N,
            "k": k,
            "profile": list(prof),
            "overhang_frac": str(oh),
            "overhang_dec12": f"{oh_dec:.12f}",
            "blocks": [
                {"id": i, "layer": layers[i], "role": roles[i], "X": str(X[i])}
                for i in range(1, N + 1)
            ],
            "contacts": [
                {"i": i, "j": j, "F": str(contacts[(i, j)][0]), "tau": str(contacts[(i, j)][1])}
                for (i, j) in sorted(contacts)
            ],
        }

    if args.dump_json:
        with open(args.dump_json, "w", encoding="utf-8") as f:
            json.dump(exported, f, indent=2)
        print(f"Exported exact verified stacks to {args.dump_json}")


if __name__ == "__main__":
    main()
