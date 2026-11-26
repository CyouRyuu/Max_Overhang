"""Independently check per-block balances in Q(x) and compute root polynomials.

Optional cross-check: Python with SymPy. The interval verifier is stand-alone.
"""
import argparse
from fractions import Fraction
from pathlib import Path

from sympy import Poly, QQ
from sympy.polys.fields import field

import audit_spines as audit
import overhang_spine_check as original


def polynomial(n, role, height, variable, relation):
    # The field reduces every rational function before its numerator is taken.
    p = Poly(relation.numer.as_expr()).clear_denoms(convert=True)[1].primitive()[1]
    if p.LC() < 0:
        p = -p
    original.require(p.is_irreducible, f"reducible polynomial for {role}, N={n}")
    return {"N": n, "role": role, "K": height, "variable": variable,
            "degree": p.degree(), "irreducible_over_Q": True,
            "coefficients_descending": " ".join(map(str,p.all_coeffs()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("results"))
    folder = parser.parse_args().out
    folder.mkdir(parents=True,exist_ok=True)
    F,x = field("x",QQ)
    half = F(Fraction(1,2))
    rec = original.recurrence(x,7,half)
    polynomials,identities = [],[]
    for n in range(4,20):
        winner = max(range(1,n+1), key=lambda h: original.loaded_spine(h,n)[1].lo)
        polynomials.append(polynomial(n,"loaded optimum",winner,"t1",rec[winner]-n))
        if n==5:
            k = 3
            blocks,contacts,step = audit.bridgeless_five(x)
            polynomials.append(polynomial(n,"ordinary witness",k,"t2",2*x*x-15))
        else:
            if n==7:
                ys,z = audit.FAMILY7
            elif n in original.SEEDS:
                ys,z,_ = original.SEEDS[n]
            else:
                ys,z = audit.EXTRA_SEEDS[n]
            k = len(ys)+1
            t = rec[:k]+[F(n)]
            blocks,contacts,step = audit.shielded_stack(
                n,t,[F(Fraction(s)) for s in ys],None if z is None else F(Fraction(z)))
            if n==7:
                polynomials.append(polynomial(n,"ordinary witness",k,"t1",rec[k]-n))
        original.require(len(blocks)==n,"wrong block count")
        # Each block starts with its own downward unit weight and moment.
        # At each interface add the upward support to the upper block and its
        # opposite to the lower block. The table has no balance equation.
        force = {b:-F.one for b in blocks}
        moment = {b:-(edge+half) for b,(_,edge) in blocks.items()}
        for c in contacts:
            torque = c.force*c.point
            force[c.upper] += c.force
            moment[c.upper] += torque
            if c.lower!="table":
                force[c.lower] -= c.force
                moment[c.lower] -= torque
        for b in sorted(blocks):
            original.require(force[b]==moment[b]==F.zero,f"block equilibrium: {n} {b}")
            identities.append({"N":n,"K":k,"object":b,"force_residual":0,"moment_residual":0,
                               "centroid_sum":"","overhang_minus_sum_of_steps":""})
        centroid = sum(edge+half for _,edge in blocks.values())
        telescoping = blocks["B1"][1]+1-sum(step.values())
        original.require(centroid==telescoping==F.zero,f"global identities: {n}")
        identities.append({"N":n,"K":k,"object":"whole stack","force_residual":"",
                           "moment_residual":"","centroid_sum":0,"overhang_minus_sum_of_steps":0})
        print(f"N={n:2d}: exact balance on {n} blocks; primitive root polynomial verified.")
    audit.dump(folder/"polynomials.csv",polynomials)
    audit.dump(folder/"identities.csv",identities)


if __name__=="__main__":
    main()
