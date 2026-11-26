"""Verify a finite loaded-spine construction from the submitted overhang paper.

Usage: python3 overhang_spine_check.py   (standard library only)

This independent check concerns eta=0 and N=4,6,8,...,16 (including odd N>=9).
It uses Proposition 22(i), equations (48)-(50) and (60)-(65), and Table A1.
For each N it isolates the maximizing height and its algebraic load vector.
It retains the printed shield and bridge seeds and verifies that these exact
loads can be realized by ordinary unit blocks. Every comparison used as a
certificate has outward-rounded rational endpoints. Decimal arithmetic only
suggests a starting bracket and formats the results.

Equilibrium follows algebraically from the construction. In spine row r the
force identity is t[r] = 1 + omega[r] + t[r-1]. The moment equation reduces
to t[r]*d[r] = omega[r] + 1/2. For each intermediate shield the upward
force and moment from the shield above it are defined by the two balance
equations (62). At the top shield, summing the load increments gives
u[2] + omega[2] - 1 = u[1] + v, where v=N-2k and u[1]=1-omega[1].
The bridge and tower centers solve the two remaining moment equations
(63)-(65). These identities also hold at the isolated algebraic root.
The interval checks below establish feasibility of the resulting geometry,
force signs, the unique principal block and its geometric support chain.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations


HALF = Fraction(1, 2)
GRID = 10**60


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


class Interval:
    """Exact rational endpoints rounded outward to a common decimal grid."""

    def __init__(self, lo, hi=None):
        if isinstance(lo, Interval) and hi is None:
            self.lo, self.hi = lo.lo, lo.hi
            return
        lo = Fraction(lo)
        hi = lo if hi is None else Fraction(hi)
        require(lo <= hi, "reversed interval")
        self.lo = Fraction((lo * GRID) // 1, GRID)
        self.hi = Fraction(-((-hi * GRID) // 1), GRID)

    def __add__(self, other):
        other = Interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-Interval(other))

    def __rsub__(self, other):
        return Interval(other) + (-self)

    def __mul__(self, other):
        other = Interval(other)
        values = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Interval(other)
        require(other.lo > 0 or other.hi < 0, "divisor may contain zero")
        return self * Interval(1 / other.hi, 1 / other.lo)


def recurrence(first, height, half):
    values = [first * 0, first]
    for _ in range(1, height):
        previous, current = values[-2:]
        values.append(current * current / (previous + half))
    return values


def loaded_spine(height, weight):
    if height == 1:
        return [Interval(0), Interval(weight)], Interval(1) - Interval(HALF) / weight

    # Produce a candidate interval, then check its endpoints independently.
    with localcontext() as ctx:
        ctx.prec = 52
        lo, hi = Decimal("0.5"), Decimal(weight)
        for _ in range(150):
            mid = (lo + hi) / 2
            if recurrence(mid, height, Decimal("0.5"))[-1] < weight:
                lo = mid
            else:
                hi = mid
        denominator = 10**26
        guess = int((lo + hi) * denominator / 2)
    left = Fraction(guess - 1, denominator)
    right = Fraction(guess + 2, denominator)
    terminal_lo = recurrence(Interval(left), height, HALF)[-1]
    terminal_hi = recurrence(Interval(right), height, HALF)[-1]
    require(terminal_lo.hi < weight < terminal_hi.lo, "root is not bracketed")

    # Proposition 22 gives a unique root between these rational endpoints.
    loads = recurrence(Interval(left, right), height, HALF)
    loads[-1] = Interval(weight)
    value = sum(1 - (a + HALF) / b for a, b in zip(loads, loads[1:]))
    return loads, value


# Table A1. Height equals one plus the number of shield edges.
# Columns are shield left edges, bridge contact (when supplied), and the
# manuscript's rational lower bound for comparison.
SEEDS = {
    4:  (("-1",), None, "1651401/1414000"),
    6:  (("-1", "-6/5"), None, "21115183/14697000"),
    8:  (("-17/20", "-1", "-7/5"), None, "247112620283/151462819200"),
    9:  (("-13/20", "-1", "-7/5"), "11/40", "158860140542171/92615086838000"),
    10: (("-11/20", "-1", "-7/5"), "2/5", "1483400102771/830047920000"),
    11: (("-11/20", "-9/10", "-6/5", "-303/200"), "2/5", "335947378939/180735464000"),
    12: (("-1/2", "-4/5", "-6/5", "-3/2"), "17/40", "12287515301/6382821900"),
    13: (("-3/5", "-31/50", "-11/10", "-29/20"), "2/5", "120400867309071/60670358899150"),
    14: (("-1/2", "-3/5", "-21/20", "-7/5"), "1/2", "23072966927875223/11320172285996000"),
    15: (("-1/2", "-101/200", "-19/20", "-13/10", "-8/5"), "1/2",
         "58035616013404873/27729809664879000"),
    16: (("-9/20", "-1/2", "-9/10", "-5/4", "-31/20"), "11/20",
         "41648849852234795719/19427194065076512000"),
}


def check_stack(n, loads, shield_edges, bridge_seed):
    k = len(loads) - 1
    tower_size = n - 2 * k
    omega = {r: loads[r] - loads[r-1] - 1 for r in range(1, k+1)}
    require(all(w.lo > 0 for w in omega.values()), "nonpositive spine load increment")
    step = {r: (omega[r] + HALF) / loads[r] for r in omega}
    spine_edge = {k: step[k] - 1}
    for r in reversed(range(1, k)):
        spine_edge[r] = spine_edge[r+1] + step[r]

    shield_edge = {r: Interval(y) for r, y in enumerate(shield_edges, start=1)}
    force = {k: Interval(0), 1: 1 - omega[1]}
    torque = {k: Interval(0)}
    for r in reversed(range(2, k)):
        force[r] = force[r+1] + omega[r+1] - 1
        torque[r] = torque[r+1] + omega[r+1] * spine_edge[r+1] - (shield_edge[r] + HALF)
    demand = torque[2] + omega[2] * spine_edge[2] - (shield_edge[1] + HALF)
    bridge_point = Interval(bridge_seed) if tower_size else demand / force[1]
    bridge_center = force[1] * bridge_point + omega[1] * spine_edge[1]

    # Work in left-edge coordinates. Each unit block occupies [edge, edge+1].
    edges, layers = {}, {}
    def place(name, layer, left_edge):
        edges[name] = left_edge
        layers[name] = layer

    def spine(r):
        return ("spine", r)

    def shield(r):
        return ("shield", r)

    for r, edge in spine_edge.items():
        place(spine(r), k+1-r, edge)
    for r, edge in shield_edge.items():
        place(shield(r), k+1-r, edge)
    bridge = ("bridge", 0)
    place(bridge, k+1, bridge_center - HALF)
    if tower_size:
        tower_center = (demand - force[1] * bridge_point) / tower_size
        for j in range(1, tower_size+1):
            place(("tower", j), k+j, tower_center - HALF)
    require(len(edges) == n, "wrong number of blocks")

    def check_contact(upper, lower, load, point=None, lower_side=None):
        require(load.lo > 0 and layers[upper] == layers[lower] + 1, "invalid contact")
        a, b = edges[upper], edges[lower]
        # All force-bearing contacts must have positive overlap (also in PZ's model).
        require((1+a-b).lo > 0 and (1+b-a).lo > 0, "contact overlap may be zero")
        if lower_side is None:
            require((point-b).lo >= 0 and (b+1-point).lo >= 0, "point outside lower block")
        else:
            # By construction this coordinate is identically the indicated endpoint.
            require(lower_side in (0, 1) and point is None, "bad endpoint specification")
            point = b + lower_side
        require((point-a).lo >= 0 and (a+1-point).lo >= 0, "point outside upper block")

    base = edges[spine(k)]
    require(base.hi < 0 and (base+1).lo >= 0, "table contact at zero unavailable")
    for r in range(1, k):
        check_contact(spine(r), spine(r+1), loads[r], lower_side=1)
        check_contact(shield(r), spine(r+1), omega[r+1], lower_side=0)
    check_contact(bridge, spine(1), omega[1], lower_side=0)
    check_contact(bridge, shield(1), force[1], point=bridge_point)
    for r in range(2, k):
        check_contact(shield(r-1), shield(r), force[r], point=torque[r]/force[r])
    for j in range(1, tower_size+1):
        lower = shield(1) if j == 1 else ("tower", j-1)
        check_contact(("tower", j), lower, Interval(tower_size-j+1), point=tower_center)

    for a, b in combinations(edges, 2):
        if layers[a] == layers[b]:
            dist = edges[a] - edges[b]
            require(max(dist.lo, -dist.hi) >= 1, "same-layer blocks may overlap")
    top = spine(1)
    require(all((edges[top]-e).lo > 0 for name, e in edges.items() if name != top),
            "the principal block may not be unique")
    for r in range(1, k):
        s = spine(r)
        for name, edge in edges.items():
            if layers[name] == layers[s]-1 and name != spine(r+1):
                dist = edges[s] - edge
                require(max(dist.lo, -dist.hi) > 1, "spine may rest on another block")


def main():
    print("N   maximizing height   loaded value (approx.)   loaded minus Table A1 (approx.)")
    for n, (shields, bridge, zq) in SEEDS.items():
        k = len(shields) + 1
        candidates = {h: loaded_spine(h, n) for h in range(1, n+1)}
        loads, best = candidates[k]
        require(all(best.lo > s.hi for h, (_, s) in candidates.items() if h != k),
                f"N={n}: maximizing height not separated")
        check_stack(n, loads, shields, bridge)
        gap = best - Fraction(zq)
        require(0 <= gap.lo <= gap.hi <= Fraction(31, 10**9), "stated rational gap exceeded")
        print(f"{n:2d}  {k:8d}            {float(best.hi):.12f}             {float(gap.hi):.6g}")
    print("All eleven finite spinal realizations passed.")


if __name__ == "__main__":
    main()
