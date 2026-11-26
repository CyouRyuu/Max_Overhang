"""Negative tests of the two auditors, pd_certificate.py and ns_audit.py.

Each test copies an archived record (results/raw_logs, results/nonspinal), changes one field, and runs the auditor,
which must reject the tampered record (nonzero exit code); the untouched records must pass.  Needs no solver.

  python3 src/test_auditors.py
"""
import csv, glob, json, os, shutil, subprocess, sys, tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = f"{REPO}/results"
SRC = f"{REPO}/src"
ANCHOR = f"{R}/anchor_n4/M0cert_N4_s42_result.json"
CERT6 = f"{R}/certificates/cert_N6.json"
CERT15 = f"{R}/certificates/cert_N15.json"
results = []


def run(args):
    p = subprocess.run([sys.executable, *args], capture_output=True, text=True, cwd=REPO)
    return p.returncode


def tamper_csv(src_csv, dst_dir, edit):
    os.makedirs(dst_dir, exist_ok=True)
    base = os.path.basename(src_csv)
    rows = list(csv.DictReader(open(src_csv)))
    fields = list(rows[0].keys())
    edit(rows)
    dst = os.path.join(dst_dir, base)
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    rj = src_csv.replace("_profiles.csv", "_result.json")
    if rj != src_csv and os.path.exists(rj):
        shutil.copy(rj, dst_dir)
    return dst


def expect(name, rc, want_fail=True):
    ok = (rc != 0) if want_fail else (rc == 0)
    results.append((name, rc, ok))
    print(f"{'PASS' if ok else 'FAIL'}  rc={rc}  {name}")


tmp = tempfile.mkdtemp()
c5 = glob.glob(f"{R}/raw_logs/cert_N5/*_profiles.csv")[0]
c7 = glob.glob(f"{R}/raw_logs/cert_N7/*_profiles.csv")[0]
ns5 = f"{R}/raw_logs/ns_N5/ns_N5_rel0.01.csv"
ns7 = f"{R}/raw_logs/ns_N7/ns_N7_rel0.01.csv"
ns16 = f"{R}/nonspinal/ns_N16_rel0.01.csv"
t016 = f"{R}/raw_logs/ns_N16_t0/ns_N16_rel0.01_only_t0.csv"
out = os.path.join(tmp, "o")


def pd(N, csvp, cut, prev):
    return run([f"{SRC}/pd_certificate.py", str(N), csvp, "--cutoff", str(cut), "--prev-cert", prev, "--out", out])


def ns(N, csvp, *extra):
    return run([f"{SRC}/ns_audit.py", str(N), csvp, *extra, "--out", out + ".json"])


# controls: untouched inputs pass
expect("control pd N=5", pd(5, c5, 1.3045649, ANCHOR), want_fail=False)
expect("control pd N=7", pd(7, c7, 1.5300553, CERT6), want_fail=False)
expect("control ns N=5", ns(5, ns5, "--prev-cert", ANCHOR), want_fail=False)
expect("control ns N=7", ns(7, ns7, "--prev-cert", CERT6), want_fail=False)
expect("control ns N=16 split", ns(16, ns16, "--split", t016), want_fail=False)

# BYPASS-A: M0 result of N = 4 passed as the previous certificate of N = 16, and an M0 result without N
expect("A: ns N=16 with the N=4 M0 result as prev-cert", ns(16, ns16, "--split", t016, "--prev-cert", ANCHOR))
m0 = json.load(open(ANCHOR)); m0.pop("N")
noN = os.path.join(tmp, "m0_noN.json"); json.dump(m0, open(noN, "w"))
expect("A: pd N=5 with an M0 result lacking N", pd(5, c5, 1.3045649, noN))
expect("A: ns N=5 with an M0 result lacking N", ns(5, ns5, "--prev-cert", noN))


def setf(i, **kv):
    def e(rows):
        for k, v in kv.items():
            rows[i][k] = v
    return e


# BYPASS-C/F: NaN and infinite bounds
expect("C: CUTOFF bound nan (second row)", pd(5, tamper_csv(c5, f"{tmp}/c1", setf(1, bound="nan")), 1.3045649, ANCHOR))
expect("C: CUTOFF bound -inf", pd(5, tamper_csv(c5, f"{tmp}/c2", setf(0, bound="-inf")), 1.3045649, ANCHOR))
expect("C: CUTOFF bound empty", pd(5, tamper_csv(c5, f"{tmp}/c3", setf(1, bound="")), 1.3045649, ANCHOR))

# BYPASS-D/E1: blank vars, wrong closure kinds
expect("D: vars blank on a solver row", pd(5, tamper_csv(c5, f"{tmp}/d1", setf(0, vars="")), 1.3045649, ANCHOR))
rows7 = list(csv.DictReader(open(c7)))
ia = next(i for i, r in enumerate(rows7) if r["closure"].startswith("analytic"))
isv = next(i for i, r in enumerate(rows7) if r["closure"] == "solver")
expect("E1: analytic row labelled solver", pd(7, tamper_csv(c7, f"{tmp}/e1", setf(ia, closure="solver")), 1.5300553, CERT6))
expect("E1: analytic row with the wrong analytic kind",
       pd(7, tamper_csv(c7, f"{tmp}/e2", setf(ia, closure="analytic: no candidate layer")), 1.5300553, CERT6))
expect("E1: solver row labelled analytic",
       pd(7, tamper_csv(c7, f"{tmp}/e3", setf(isv, closure="analytic: energy test")), 1.5300553, CERT6))
expect("E1: solver row with closure blank", pd(7, tamper_csv(c7, f"{tmp}/e4", setf(isv, closure="")), 1.5300553, CERT6))

# BYPASS-B/E2: Gurobi status, objective bound, cutoff columns
expect("B: CUTOFF with Gurobi status 9", pd(5, tamper_csv(c5, f"{tmp}/b1", setf(0, gurobi_status="9")), 1.3045649, ANCHOR))
expect("B: CUTOFF with obj_bound 9.99", pd(5, tamper_csv(c5, f"{tmp}/b2", setf(0, obj_bound="9.99")), 1.3045649, ANCHOR))
expect("B: CUTOFF with obj_bound nan", pd(5, tamper_csv(c5, f"{tmp}/b3", setf(0, obj_bound="nan")), 1.3045649, ANCHOR))
expect("B: recorded cutoff 9.99", pd(5, tamper_csv(c5, f"{tmp}/b4", setf(0, cutoff="9.99", obj_bound="1.3")), 1.3045649, ANCHOR))
expect("B: Gurobi status on an analytic row",
       pd(7, tamper_csv(c7, f"{tmp}/b5", setf(ia, gurobi_status="3")), 1.5300553, CERT6))


def ns_tamper(src, d, pick, **kv):
    rows = list(csv.DictReader(open(src)))
    i = next(i for i, r in enumerate(rows) if pick(r))
    return tamper_csv(src, d, setf(i, **kv))


cut = lambda r: r["status"] == "CUTOFF"
inf = lambda r: r["status"] == "INFEASIBLE" and r["closure"] == "solver"
emp = lambda r: r["status"] == "EMPTY"
expect("E2: ns CUTOFF with Gurobi status 9", ns(7, ns_tamper(ns7, f"{tmp}/n1", cut, gurobi_status="9"), "--prev-cert", CERT6))
expect("E2: ns CUTOFF with obj_bound 9.99", ns(7, ns_tamper(ns7, f"{tmp}/n2", cut, obj_bound="9.99"), "--prev-cert", CERT6))
expect("E2: ns CUTOFF with recorded cutoff 9.99",
       ns(7, ns_tamper(ns7, f"{tmp}/n3", cut, cutoff="9.99", obj_bound="1.0"), "--prev-cert", CERT6))
expect("E2: ns INFEASIBLE with Gurobi status 6", ns(5, ns_tamper(ns5, f"{tmp}/n4", inf, gurobi_status="6"), "--prev-cert", ANCHOR))
expect("E2: ns EMPTY with a Gurobi status", ns(5, ns_tamper(ns5, f"{tmp}/n5", emp, gurobi_status="6"), "--prev-cert", ANCHOR))
expect("E2: ns solver row relabelled bounds-before-solving as CUTOFF",
       ns(7, ns_tamper(ns7, f"{tmp}/n6", cut, closure="bounds before solving", gurobi_status=""), "--prev-cert", CERT6))
expect("E2: ns t0 split CUTOFF with Gurobi status 9",
       ns(16, ns16, "--split", ns_tamper(t016, f"{tmp}/n7", cut, gurobi_status="9")))

# regressions of earlier attacks: EMPTY / TIE loopholes
expect("old: ns CUTOFF row relabelled EMPTY", ns(7, ns_tamper(ns7, f"{tmp}/o1", cut, status="EMPTY", closure="analytic: J empty",
                                                               gurobi_status=""), "--prev-cert", CERT6))
expect("old: t0 split CUTOFF row relabelled TIE",
       ns(16, ns16, "--split", ns_tamper(t016, f"{tmp}/o2", cut, status="TIE", closure="analytic: t0 = l*", gurobi_status="")))

bad = [r for r in results if not r[2]]
print(f"\n{len(results) - len(bad)} / {len(results)} tests behave as expected")
shutil.rmtree(tmp)
sys.exit(1 if bad else 0)
