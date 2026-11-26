"""Summary of the audited chain N = 5..16 (output of scripts/audit_archive.sh); needs no solver.

  python3 src/chain_summary.py runs/audit

Reads cert_N<N>.json (pd_certificate.py) and ns_N<N>.json (ns_audit.py) for N = 5..16 and prints, per N, the
certificate, the S2+ condition, the non-spinal exclusion and its condition for stacks with n_2 = 1, together
with the sources of the bounds of size N - 1.  Exits 1 unless every entry holds.
"""
import json
import os
import sys


def main():
    d = sys.argv[1]
    ok_all = True
    print(f"{'N':>3} {'cert':>5} {'S2+ margin':>11} {'zbar(N-1) from':<16} {'ns':>5} {'n2=1 condition (margins)':<44}"
          f" {'analytic closures (energy / no cand.)':<20}")
    for N in range(5, 17):
        c = json.load(open(os.path.join(d, f"cert_N{N}.json")))
        n = json.load(open(os.path.join(d, f"ns_N{N}.json")))
        src = os.path.basename(c["S2plus_upper_bound_source"])
        ch = n["n2_equals_1_condition"]
        parts = []
        for k in ("general", "sharper"):
            if f"{k}_margin" in ch:
                parts.append(f"{k} {ch[f'{k}_margin']:+.5f}{'' if ch[f'{k}_ok'] else ' (fails)'}")
        ok = c["certified"] and c["S2plus_ok"] and n["certified"] and ch.get("ok", False)
        ok_all &= ok
        an = c["analytic_closures"]
        print(f"{N:>3} {str(c['certified']):>5} {c['S2plus_margin']:>+11.5f} {src:<16} {str(n['certified']):>5} "
              f"{'; '.join(parts):<44} {an['energy']:>5} / {an['no_candidate_layer']}")
    print("CHAIN N = 5..16 " + ("HOLDS" if ok_all else "DOES NOT HOLD"))
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
