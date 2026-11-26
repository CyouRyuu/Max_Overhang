"""Check a feasible stack reported by a solver run (a witness), in floating point, without a solver.

Reads either a *_witnesses.json file of ns_check.py or a *_result.json file of profile_run.py / certify_run.py
(field "stack") and reports, for zero stability margin (W = w = 1/2), the largest violation of
  * vertical force balance and moment balance of every block (unit weight at the block centre),
  * the contact intervals |P - X_i| <= 1/2 and |P - X_j| <= 1/2 (table: P <= 0),
  * nonnegative forces, adjacent-layer contacts and unit spacing within each layer,
together with the overhang max_i X_i + 1/2.  Violations of the order of the solver tolerances (1e-6) are expected.

  python3 src/check_witness.py results/reverse_controls/ns_N9_g1.6000000_witnesses.json
"""
import sys
import json


def check(stack):
    X = {r["block"]: r["center_x"] for r in stack}
    L = {r["block"]: r["layer"] for r in stack}
    up = {i: 0.0 for i in X}; mom = {i: 0.0 for i in X}
    viol = {"force_balance": 0.0, "moment_balance": 0.0, "contact_interval": 0.0, "negative_force": 0.0,
            "layer_link": 0.0, "same_layer_spacing": 0.0}
    for r in stack:
        i = r["block"]
        for c in r.get("contacts", []):
            j, F, P = c["on"], c["F"] or 0.0, c["P"]
            if P is None:
                raise SystemExit("witness without contact points; rerun with the current code")
            viol["negative_force"] = max(viol["negative_force"], -F)
            up[i] += F; mom[i] += F * P
            if j == 0:
                viol["contact_interval"] = max(viol["contact_interval"], P, abs(P - X[i]) - 0.5)
                viol["layer_link"] = max(viol["layer_link"], abs(L[i] - 1))
            else:
                up[j] -= F; mom[j] -= F * P
                viol["contact_interval"] = max(viol["contact_interval"], abs(P - X[i]) - 0.5, abs(P - X[j]) - 0.5)
                viol["layer_link"] = max(viol["layer_link"], abs(L[i] - L[j] - 1))
    for i in X:
        viol["force_balance"] = max(viol["force_balance"], abs(up[i] - 1.0))
        viol["moment_balance"] = max(viol["moment_balance"], abs(mom[i] - X[i]))
    for i in X:
        for k in X:
            if i < k and L[i] == L[k]:
                viol["same_layer_spacing"] = max(viol["same_layer_spacing"], 1.0 - abs(X[i] - X[k]))
    return max(X.values()) + 0.5, viol, spinality(X, L)


def spinality(X, L, tol=1e-6):
    """Spinal status under the paper's closed resting relation (|X_i - X_j| <= 1, up to tol) and under
    Paterson and Zwick's strict relation (|X_i - X_j| < 1 - tol)."""
    z = max(X.values())
    lstar = min(L[i] for i in X if X[i] >= z - tol)
    p = max((i for i in X if L[i] == lstar), key=lambda i: X[i])
    out = {}
    for name, rel in (("closed", lambda d: d <= 1 + tol), ("strict", lambda d: d < 1 - tol)):
        S, todo = {p}, [p]
        while todo:
            i = todo.pop()
            for j in X:
                if L[j] == L[i] - 1 and j not in S and rel(abs(X[i] - X[j])):
                    S.add(j); todo.append(j)
        out[name] = all(sum(1 for j in S if L[j] == t) == 1 for t in range(1, lstar + 1))
    return out


def main():
    d = json.load(open(sys.argv[1]))
    items = [(w["profile"], w["lstar"], w["stack"]) for w in d["witnesses"]] if "witnesses" in d else \
            [(d.get("profile", "-"), "-", d["stack"])]
    for prof, lstar, stack in items[: int(sys.argv[2]) if len(sys.argv) > 2 else 3]:
        if not stack:
            print(f"{prof} {lstar}: no feasible stack recorded")
            sys.exit(1)
        z, v, sp = check(stack)
        print(f"profile {prof} l*={lstar}: overhang {z:.7f}; spinal (closed / Paterson-Zwick): "
              f"{sp['closed']} / {sp['strict']}; max violations " +
              ", ".join(f"{k} {max(val, 0.0):.1e}" for k, val in v.items()))


if __name__ == "__main__":
    main()
