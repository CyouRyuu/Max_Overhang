"""Independent replay of the face certificates of the exception sizes N = 5 and N = 7 (Appendix C).

src/spine_faces.py bounds the loaded-spine objective f(t) = k - sum_{s=1}^{k} (t_{s-1} + 1/2) / t_s, with t_0 = 0
and t_k = N, on every face w_j = t_j - t_{j-1} - 1 = 0 of the load polytope by bisecting the initial rational box
t_s in [s, N - (k - s)] (s = 1..k-1) and writes every leaf of its bisection tree to face_leaves_N<N>.csv.  This
script re-reads the leaves and checks, in exact rational arithmetic and without sharing code with spine_faces.py:

  1. the leaves of each face form a complete bisection tree of the initial box: every split halves one coordinate,
     each internal node has both children, no leaf is an ancestor of another, and the box of every leaf equals the
     box obtained by applying its path to the initial box;
  2. a leaf marked empty contains no point with t_s - t_{s-1} >= 1 for all s and t_j - t_{j-1} = 1, shown by one of
     the linear conditions  hi_s < lo_{s-1} + 1  or (for s = j)  lo_j > hi_{j-1} + 1;
  3. on a leaf marked bound, the monotone bound  k - sum_s (lo_{s-1} + 1/2) / hi_s  (f decreases in t_{s-1} and
     increases in t_s) equals the stored value and lies strictly below the cut;
  4. the threshold is the one of the paper (1.30455 at N = 5, 1.53004 at N = 7), and every face j = 1..k is covered, except the face w_1 = 0 at N = 5, which spine_faces.py treats analytically
     (f = 12/5 - 3/(2 t_2) - t_2/5 on [2, 4], concave, with maximum (12 - sqrt 30)/5 at t_2 = sqrt(15/2)).

  python3 src/check_face_leaves.py [results/exact_spines/face_leaves]
"""
import csv
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPECT = {5: (3, {2, 3}, Fr("1.30455")), 7: (4, {1, 2, 3, 4}, Fr("1.53004"))}   # height k, faces j, threshold


def parse_box(txt):
    return [tuple(Fr(x) for x in iv.split(":")) for iv in txt.split(";")] if txt else []


def parse_path(txt):
    return tuple((int(t[:-1]), "LR".index(t[-1])) for t in txt.split(";")) if txt else ()


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results", "exact_spines", "face_leaves")
    for N, (k, faces, want_cut) in EXPECT.items():
        rows = list(csv.DictReader(open(os.path.join(d, f"face_leaves_N{N}.csv"))))
        cut = {Fr(r["cut"]) for r in rows}
        assert cut == {want_cut}, ("threshold of the face certificate", N, cut, want_cut)
        cut = want_cut
        assert {int(r["height"]) for r in rows} == {k}
        seen = {int(r["face_j"]) for r in rows}
        assert seen == faces, (N, seen)
        counts = []
        for j in sorted(faces):
            leaves = [(parse_path(r["path"]), parse_box(r["box_t1_to_tkminus1"]), r["kind"], r["bound"])
                      for r in rows if int(r["face_j"]) == j]
            box0 = [(Fr(s), Fr(N - (k - s))) for s in range(1, k)]
            paths = {p for p, _, _, _ in leaves}
            assert len(paths) == len(leaves), "duplicate leaf"
            internal = {p[:t] for p in paths for t in range(len(p))}
            assert not (internal & paths), "a leaf is an ancestor of another leaf"
            for q in internal:                   # both children, split on the same coordinate
                kids = [p for p in paths | internal if len(p) == len(q) + 1 and p[:-1] == q]
                assert len(kids) == 2 and kids[0][-1][0] == kids[1][-1][0] and {kids[0][-1][1], kids[1][-1][1]} == {0, 1}, q
            closed = 0
            for path, box, kind, b in leaves:
                bx = list(box0)
                for i, side in path:
                    lo, hi = bx[i]
                    m = (lo + hi) / 2
                    bx[i] = (lo, m) if side == 0 else (m, hi)
                assert bx == box, ("box does not match its path", N, j, path)
                lo = [Fr(0)] + [x[0] for x in bx] + [Fr(N)]
                hi = [Fr(0)] + [x[1] for x in bx] + [Fr(N)]
                empty = any(hi[s] < lo[s - 1] + 1 for s in range(1, k + 1)) or lo[j] > hi[j - 1] + 1 \
                    or hi[j] < lo[j - 1] + 1
                if kind == "empty":
                    assert empty, ("leaf marked empty is not empty", N, j, path)
                    continue
                assert kind == "bound" and not empty
                val = k - sum((lo[s - 1] + Fr(1, 2)) / hi[s] for s in range(1, k + 1))
                assert val == Fr(b) and val < cut, ("bound", N, j, path)
                closed += 1
            counts.append(f"face w_{j} = 0: {len(leaves)} leaves, {closed} closed by the bound")
        print(f"N={N}, cut {cut} = {float(cut)}: " + "; ".join(counts))
    print("FACE CERTIFICATES OF N = 5, 7 REPLAYED")


if __name__ == "__main__":
    main()
