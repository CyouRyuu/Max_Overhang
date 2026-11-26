"""Upper bounds that make the exception values optimal among spinal stacks at N = 5 and N = 7 (Appendix C).
Needs mpmath only (rigorous interval arithmetic); no optimization solver.

  python3 src/spine_faces.py

A spinal stack with spine height k and spine loads 0 = t_0 < t_1 < ... < t_k = N (t_s is the load that spine
block s passes down, counted from the top) has overhang at most  f(t) = k - sum_s (t_{s-1} + 1/2) / t_s,  and its
added loads w_s = t_s - t_{s-1} - 1 are nonnegative.  A positive w_s needs a block outside the spine resting on
spine block s, in the layer above it, so at most N - k added loads are positive.  We check:
  * every height K other than the realized one has  max f <= S^_K(N) < cut  (interval enclosure of S^_K);
  * at the height k* of the loaded-spine optimum (3 for N = 5, 4 for N = 7), where k* > N - k*, every face
    {w_j = 0} has  max f < cut  by rational box branch-and-bound, except the face w_1 = 0 at N = 5, on which
    f = 12/5 - 3/(2 t_2) - t_2/5 is concave with maximum (12 - sqrt 30)/5 at t_2 = sqrt(15/2);
  * stacks with several blocks on the table: the spine root passes at most N - 1, so their overhang is at most
    S*(N - 1) < cut.
With cut = 1.30455 (N = 5) and 1.53004 (N = 7), below the overhangs (12 - sqrt 30)/5 and S*_3(7) of the stacks
verified by check_algebraic_spines.py, these are the largest overhangs of spinal stacks.
"""
from fractions import Fraction as Fr

from mpmath import iv, mp

mp.dps = iv.dps = 50


def shat_enclosure(K, N, loads=False):
    """Rigorous enclosure of max f over the loads with t_K = N and no constraint on the increments
    (with loads=True, also the enclosures of the maximizing loads t_0..t_K)."""
    if K == 1:
        v = iv.mpf(1) - iv.mpf(0.5) / N
        return (v, [iv.mpf(0), iv.mpf(N)]) if loads else v
    lo, hi = mp.mpf(0.5), mp.mpf(N)
    for _ in range(200):
        mid = (lo + hi) / 2
        u = [iv.mpf(0), iv.mpf(mid)]
        for _ in range(1, K):
            u.append(u[-1] ** 2 / (u[-2] + iv.mpf(0.5)))
        if u[K].b < N:
            lo = mid
        elif u[K].a > N:
            hi = mid
        else:
            break
    x = iv.mpf([lo, hi])
    u = [iv.mpf(0), x]
    for _ in range(1, K):
        u.append(u[-1] ** 2 / (u[-2] + iv.mpf(0.5)))
    u[K] = iv.mpf(N)
    v = sum(1 - (u[s - 1] + iv.mpf(0.5)) / u[s] for s in range(1, K + 1))
    return (v, u) if loads else v


def face_bound(N, k, j, cut):
    """Box branch-and-bound over t_1..t_{k-1} with w_j = 0; returns (boxes closed, largest bound)."""
    free = list(range(1, k))
    box0 = [(Fr(s), Fr(N - (k - s))) for s in free]

    def bound(box):
        lo = [Fr(0)] + [b[0] for b in box] + [Fr(N)]
        hi = [Fr(0)] + [b[1] for b in box] + [Fr(N)]
        for s in range(1, k + 1):           # feasibility of the increments
            if hi[s] < lo[s - 1] + 1:
                return None
            if s == j and (lo[s] > hi[s - 1] + 1 or hi[s] < lo[s - 1] + 1):
                return None
        return k - sum((lo[s - 1] + Fr(1, 2)) / hi[s] for s in range(1, k + 1))

    work, closed, worst = [box0], 0, None
    while work:
        box = work.pop()
        b = bound(box)
        if b is None:
            continue
        if b < cut:
            closed += 1
            worst = b if worst is None else max(worst, b)
            continue
        i = max(range(len(box)), key=lambda i: box[i][1] - box[i][0])
        a, c = box[i]
        assert c - a > Fr(1, 10 ** 9), f"N={N}, face w_{j} = 0: cut not reached"
        m = (a + c) / 2
        work += [box[:i] + [(a, m)] + box[i + 1:], box[:i] + [(m, c)] + box[i + 1:]]
    return closed, worst


def main():
    for N, cut, kstar, realized in ((5, Fr("1.30455"), 3, 3), (7, Fr("1.53004"), 4, 3)):
        c = iv.mpf(cut.numerator) / cut.denominator
        # the loaded-spine optimum S*(N) sits at height kstar, with all added loads positive, above S(N)
        vstar, ustar = shat_enclosure(kstar, N, loads=True)
        assert all((ustar[s] - ustar[s - 1] - 1).a > 0 for s in range(1, kstar + 1)), "added loads"
        assert vstar.a > max(shat_enclosure(K, N).b for K in range(1, N + 1) if K != kstar), "height"
        for K in range(1, N + 1):
            if K in (kstar, realized):
                continue
            v = shat_enclosure(K, N)
            assert v.b < c, f"N={N}, height {K}"
        prev = max(shat_enclosure(K, N - 1).b for K in range(1, N))
        assert prev < c, "several roots"
        line = [f"N={N}: S*(N) = S*_{kstar}(N) > {mp.nstr(mp.floor(mp.mpf(vstar.a) * 10 ** 10) / 10 ** 10, 11)}; other heights and stacks with "
                f"several roots below {float(cut)}"]
        for j in range(1, kstar + 1):
            if N == 5 and j == 1:
                # w_1 = 0: t_1 = 1, t_3 = 5, f = 12/5 - 3/(2 t_2) - t_2/5 on [2, 4], concave, max at sqrt(15/2)
                continue
            closed, worst = face_bound(N, kstar, j, cut)
            up = -((-worst * 10 ** 7) // 1) / 10 ** 7     # rounded upward at the 7th decimal
            line.append(f"face w_{j} = 0: {closed} boxes, largest bound {float(up):.7f} (rounded upward)")
        print("; ".join(line))
    print("N=5: face w_1 = 0 is concave with maximum (12 - sqrt 30)/5 = 1.304554884990 at t_2 = sqrt(15/2)")
    print("SPINAL OPTIMA AT N = 5, 7 VERIFIED")


if __name__ == "__main__":
    main()
