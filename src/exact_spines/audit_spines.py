"""Produce exact rational enclosures, geometric margins, and exception-face bounds.

Run next to overhang_spine_check.py using python3 audit_spines.py --out results.
This verifier uses only the Python standard library and the supplied original.
Block names B, H, G, T stand for spine, shield, bridge and tower, respectively.
"""
import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt
from pathlib import Path

import overhang_spine_check as original
from overhang_spine_check import HALF, Interval as I, require


EXTRA_SEEDS = {
    17: (("-2/5", "-4/5", "-6/5", "-13/10", "-8/5"), "3/5"),
    18: (("-3/10", "-4/5", "-11/10", "-7/5", "-8/5"), "3/5"),
    19: (("-9/40", "-29/40", "-43/40", "-27/20", "-3/2", "-13/8"), "13/20"),
}
FAMILY7 = (("-3/4", "-6/5"), "1/5")


@dataclass
class Contact:
    upper: str
    lower: str
    force: object
    point: object
    lower_side: str = ""


def shielded_stack(n, t, ys, z_seed):
    """Eqs. (60)-(65); works for real intervals and exact rational functions."""
    k = len(t)-1
    v = n-2*k
    zero = t[1]*0
    omega = {s: t[s]-t[s-1]-1 for s in range(1,k+1)}
    d = {s: (omega[s]+HALF)/t[s] for s in omega}
    x = {k: d[k]-1}
    for s in reversed(range(1,k)):
        x[s] = x[s+1]+d[s]
    y = dict(enumerate(ys,1))
    u, m = {k: zero, 1: 1-omega[1]}, {k: zero}
    for s in reversed(range(2,k)):
        u[s] = u[s+1]+omega[s+1]-1
        m[s] = m[s+1]+omega[s+1]*x[s+1]-(y[s]+HALF)
    demand = m[2]+omega[2]*x[2]-(y[1]+HALF)
    zp = z_seed if v else demand/u[1]
    bridge_center = u[1]*zp+omega[1]*x[1]
    blocks = {f"B{s}": (k+1-s, e) for s,e in x.items()}
    blocks.update({f"H{s}": (k+1-s, e) for s,e in y.items()})
    blocks["G"] = (k+1, bridge_center-HALF)
    contacts = [Contact(f"B{k}", "table", zero+n, zero, "right")]
    for s in range(1,k):
        contacts.append(Contact(f"B{s}", f"B{s+1}", t[s], x[s+1]+1, "right"))
        contacts.append(Contact(f"H{s}", f"B{s+1}", omega[s+1], x[s+1], "left"))
    contacts.extend([Contact("G", "B1", omega[1], x[1], "left"),
                     Contact("G", "H1", u[1], zp)])
    for s in range(2,k):
        contacts.append(Contact(f"H{s-1}", f"H{s}", u[s], m[s]/u[s]))
    if v:
        cx = (demand-u[1]*zp)/v
        for j in range(1,v+1):
            blocks[f"T{j}"] = (k+j, cx-HALF)
            contacts.append(Contact(f"T{j}", "H1" if j == 1 else f"T{j-1}", zero+v-j+1, cx))
    return blocks, contacts, d


def bridgeless_five(t2):
    """Eq. (66), with the Table A1 coordinate X(C2)=-18/25."""
    zero = t2*0
    d = {1: zero+HALF, 2: 1-(zero+1+HALF)/t2, 3: (4+HALF-t2)/5}
    x3 = d[3]-1
    x2 = x3+d[2]
    x1 = x2+d[1]
    w2, w3, u2 = t2-2, 4-t2, 3-t2
    c2 = zero+Q(-18,25)
    m = w3*x3-c2
    c1 = w2*x2+m
    blocks = {"B3": (1,x3), "B2": (2,x2), "B1": (3,x1),
              "C2": (2,c2-HALF), "C1": (3,c1-HALF)}
    contacts = [Contact("B3","table",zero+5,zero,"right"),
                Contact("B2","B3",t2,x3+1,"right"), Contact("C2","B3",w3,x3,"left"),
                Contact("B1","B2",zero+1,x2+1,"right"), Contact("C1","B2",w2,x2,"left"),
                Contact("C1","C2",u2,m/u2)]
    return blocks, contacts, d


def exact_sqrt_five():
    scale = 10**26
    a = isqrt(15*scale**2//2)
    lo, hi = Q(a,scale), Q(a+1,scale)
    require(2*lo**2 < 15 < 2*hi**2, "bad square-root enclosure")
    return I(lo,hi)


def imin(*values):
    return I(min(v.lo for v in values), min(v.hi for v in values))


def absolute(v):
    lower = max(v.lo, -v.hi, 0)
    return I(lower, max(-v.lo, v.hi))


def audit_geometry(n, k, t, blocks, contacts):
    rows = []
    def record(group, item, value, strict=False):
        require(value.lo > 0 if strict else value.lo >= 0, f"{n}: {group}, {item}")
        rows.append({"N": n, "group": group, "item": item, "lo": value.lo, "hi": value.hi})

    require(len(blocks)==n and [s for s,(h,_) in blocks.items() if h==1]==[f"B{k}"], "single root")
    require(len({(c.upper,c.lower) for c in contacts})==len(contacts), "repeated contact")
    require({c.upper for c in contacts}==set(blocks), "unsupported block")
    for s in range(1,k+1):
        record("added_load", f"B{s}", t[s]-t[s-1]-1, strict=(n!=5 or s!=1))
    for c in contacts:
        name = c.upper+"/"+c.lower
        a, p = blocks[c.upper][1], c.point
        record("force",name,c.force,True)
        if c.lower=="table":
            require(blocks[c.upper][0]==1 and p.lo==p.hi==0 and (a+1).lo>=0, "table contact")
            record("overlap",name,-a,True)
            record("point_side",name+":left",-a)
            record("point_side",name+":right",I(0))
            continue
        h,b = blocks[c.lower]
        require(blocks[c.upper][0]==h+1,"nonadjacent contact")
        record("overlap",name,1-absolute(a-b),True)
        # A designated lower-edge application is an identity in its lower block.
        # Otherwise enclose both lower-block distances by interval arithmetic.
        if c.lower_side=="left":
            left_b,right_b = I(0),I(1)
        elif c.lower_side=="right":
            left_b,right_b = I(1),I(0)
        else:
            left_b,right_b = p-b,b+1-p
        record("point_side",name+":left",imin(p-a,left_b))
        record("point_side",name+":right",imin(a+1-p,right_b))
    for a,b in combinations(blocks,2):
        if blocks[a][0]==blocks[b][0]:
            record("same_layer",a+"/"+b,absolute(blocks[a][1]-blocks[b][1])-1)
    for name,(_,edge) in blocks.items():
        if name!="B1":
            record("principal",name,blocks["B1"][1]-edge,True)
    for s in range(1,k):
        key = f"B{s}"
        h,edge = blocks[key]
        for other,(ho,e) in blocks.items():
            if ho==h-1 and other!=f"B{s+1}":
                record("other_support",key+"/"+other,absolute(edge-e)-1,True)
    return rows


def face_upper(n,k,missing,cut):
    """Cover the entire simplex with omega[missing]=0 by rational rectangles."""
    live = [j for j in range(1,k+1) if j!=missing]
    rest = Q(n-k)
    work = [[(Q(0),rest) for _ in range(len(live)-1)]]
    closed, outside, maximum = 0, 0, None
    while work:
        box = work.pop()
        if any(b[0]>c[1] for b,c in zip(box,box[1:])):
            outside += 1
            continue
        prefixes = [(Q(0),Q(0)),*box,(rest,rest)]
        lo, hi = [Q(0)], [Q(0)]
        for s in range(1,k+1):
            a,b = prefixes[sum(r<=s for r in live)]
            lo.append(s+a)
            hi.append(s+b)
        bound = k-sum((lo[s-1]+HALF)/hi[s] for s in range(1,k+1))
        if bound < cut:
            closed += 1
            maximum = bound if maximum is None else max(maximum,bound)
            continue
        j = max(range(len(box)),key=lambda i: box[i][1]-box[i][0])
        a,b = box[j]
        require(b-a>Q(1,10**8), "the proposed cutoff fails to separate this face")
        midpoint = (a+b)/2
        work.append(box[:j]+[(a,midpoint)]+box[j+1:])
        work.append(box[:j]+[(midpoint,b)]+box[j+1:])
    return {"N":n, "loaded_height":k, "zero_added_load":missing,
            "cut":cut, "closed_boxes":closed, "outside_boxes":outside, "max_closed_upper":maximum}


def dump(path,rows):
    with path.open("w",newline="",encoding="utf-8") as out:
        writer = csv.DictWriter(out,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("results"))
    folder = parser.parse_args().out
    folder.mkdir(parents=True,exist_ok=True)
    heights, load_rows, optimal_loads, detail, summary, faces = [],[],[],[],[],[]
    block_rows, contact_rows = [],[]
    previous_relaxed_upper = None
    for n in range(4,20):
        scores = {h:original.loaded_spine(h,n) for h in range(1,n+1)}
        order = sorted(scores,key=lambda h:scores[h][1].lo,reverse=True)
        winner, runner = order[:2]
        winner_t, best = scores[winner]
        rival_ub = max(v.hi for h,(_,v) in scores.items() if h!=winner)
        require(best.lo>rival_ub and scores[runner][1].hi==rival_ub,"height comparison")
        require(all((b-a-1).lo>0 for a,b in zip(winner_t,winner_t[1:])),"winning height infeasible")
        for s,load in enumerate(winner_t):
            optimal_loads.append({"N":n,"winning_loaded_height":winner,"s":s,
                                  "t_lo":load.lo,"t_hi":load.hi})
        for h,(loads,value) in scores.items():
            sign_lo = sign_hi = ""
            if h>1:
                first = loads[1]
                sign_lo = n-original.recurrence(I(first.lo),h,HALF)[-1].hi
                sign_hi = original.recurrence(I(first.hi),h,HALF)[-1].lo-n
                require(sign_lo>0 and sign_hi>0,"terminal signs")
            heights.append({"N":n,"K":h,"first_lo":loads[1].lo,"first_hi":loads[1].hi,
                            "left_sign_margin":sign_lo,"right_sign_margin":sign_hi,
                            "S_hat_lo":value.lo,"S_hat_hi":value.hi})
        if n==5:
            k = 3
            t = [I(0),I(1),exact_sqrt_five(),I(5)]
            blocks,cts,steps = bridgeless_five(t[2])
        else:
            if n==7:
                ys,z = FAMILY7
            elif n in original.SEEDS:
                ys,z,_ = original.SEEDS[n]
            else:
                ys,z = EXTRA_SEEDS[n]
            k = len(ys)+1
            t = scores[k][0]
            seeds = [I(y) for y in ys]
            zp = I(z) if z is not None else None
            blocks,cts,steps = shielded_stack(n,t,seeds,zp)
            original.check_stack(n,t,ys,z)
        value = sum(steps.values())
        require(k==winner or n in (5,7), "unexpected winning height")
        for block,(layer,edge) in sorted(blocks.items()):
            block_rows.append({"N":n,"block":block,"layer":layer,
                               "left_edge_lo":edge.lo,"left_edge_hi":edge.hi})
        for c in cts:
            contact_rows.append({"N":n,"upper":c.upper,"lower":c.lower,
                                 "force_lo":c.force.lo,"force_hi":c.force.hi,
                                 "point_lo":c.point.lo,"point_hi":c.point.hi})
        for s,load in enumerate(t):
            load_rows.append({"N":n,"witness_height":k,"s":s,"t_lo":load.lo,"t_hi":load.hi})
        rows = audit_geometry(n,k,t,blocks,cts)
        detail.extend(rows)
        summary_row = {"N":n,"witness_height":k,"winning_loaded_height":winner,
                       "runner_up_relaxed_height":runner,"height_gap_lower":best.lo-rival_ub,
                       "winning_root_lo":winner_t[1].lo,"winning_root_hi":winner_t[1].hi,
                       "value_lo":value.lo,"value_hi":value.hi,
                       "exception_multibase_upper":previous_relaxed_upper if n in (5,7) else ""}
        for group in ["added_load","force","overlap","point_side","positive_point_side",
                      "same_layer","principal","other_support"]:
            selected = ([r for r in rows if r["group"]=="point_side" and r["lo"]>0]
                        if group=="positive_point_side" else [r for r in rows if r["group"]==group])
            if selected:
                small = min(selected,key=lambda r:r["lo"])
                summary_row[group+"_lo"] = small["lo"]
                summary_row[group+"_hi"] = min(r["hi"] for r in selected)
                summary_row[group+"_item"] = small["item"]
            else:
                summary_row[group+"_lo"] = summary_row[group+"_hi"] = summary_row[group+"_item"] = ""
        sides = [r for r in rows if r["group"]=="point_side"]
        summary_row["zero_point_sides"] = sum(r["lo"]==r["hi"]==0 for r in sides)
        require(all(r["lo"]>0 or r["lo"]==r["hi"]==0 for r in sides),"undecided endpoint")
        summary.append(summary_row)
        if n in (5,7):
            cut = Q("1.30455") if n==5 else Q("1.53004")
            # A second base carries at least its own weight. For the spine
            # root the transmitted load is thus at most n-1. Increasing only
            # its final load to n-1 cannot decrease the loaded upper bound.
            require(previous_relaxed_upper<cut,"multiple-base spine upper bound")
            require(value.lo>cut and all(v.hi<cut for h,(_,v) in scores.items()
                                          if h not in (winner,k)),"exception comparison")
            for j in range(1,winner+1):
                if n==5 and j==1:
                    continue  # stationary point in the closed-form calculation
                faces.append(face_upper(n,winner,j,cut))
        previous_relaxed_upper = best.hi
        print(f"N={n:2d}   loaded k={winner}   realized k={k}   overhang~{float(value.lo):.12f}")
    for name,rows in [("heights",heights),("loads",load_rows),("optimal_loads",optimal_loads),
                      ("blocks",block_rows),("contacts",contact_rows),("margins",detail),
                      ("summary",summary),("exception_faces",faces)]:
        dump(folder/(name+".csv"),rows)


if __name__=="__main__":
    main()
