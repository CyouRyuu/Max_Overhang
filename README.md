# Max_Overhang

[![DOI](https://zenodo.org/badge/1471033814.svg)](https://doi.org/10.5281/zenodo.24518950)

Code, certificates and exact data for the paper

> A. Zha, Q. Zhuang, Q. Chang.
> **Certified Global Optimization for the Single-Root Maximum Overhang Problem.** 2026 (preprint).

The repository contains everything needed to reproduce the paper's results:

* **certified upper bounds** $`\overline{z}(N)`$ for $`N = 5,\dots,16`$ (spinal-cut layer-profile decomposition, solved
  with Gurobi, valid within its numerical tolerances);
* **non-spinal exclusion** above $`\gamma_{\mathrm{ns}}(N)`$, the value $`0.99\,S^*(N)`$ rounded upward at the 7th decimal, for $`N = 5,\dots,16`$;
* **exact rational feasible stacks** (the lower-bound witnesses) for $`N = 4,\dots,16`$, verified in exact arithmetic with zero residual
  (Appendix A), together with independent checks of the base-singleton (S2+) and non-spinal chains;
* **exact spinal optima** (Appendix C): stacks with exact algebraic loads that attain the loaded-spine bound, so
  $`S(N) = S^*(N)`$ for $`N = 4, 6, 8, \dots, 19`$, and the closed forms $`S(5) = (12 - \sqrt{30})/5`$ and
  $`S(7) = \hat S_3(7)`$ at the exception sizes, checked by two independent verifiers without any solver;
* the **baselines and cross-checks** of the methodological-progression table ($`N = 12, 13`$);
* the **archived solver records** of all runs reported in the paper.

## Layout

```
src/                 Python sources
  certify_m0.py        baseline formulation M0 (monolithic runs of the progression table)
  certify_run.py       shared model builders (M0, M1u) and the structural reductions; N = 4 anchor; it also contains
                       an experimental centre-relative moment model (Mrel) that the paper does not use
  ablation_m1_order_prune.py, maxOverhang_proposed.py
                       order-pruned formulation (the base of model M1u) and shared solve/report utilities
  profile_run.py       layer-profile decomposition over P_2(N), with energy filter (--energy) and spinal cuts (--spine)
  pd_certificate.py    assembles and audits a decomposition certificate: every expected profile present and closed,
                       no unexpected profile, repeated records reported, run cutoff within the certificate cutoff,
                       S2+ condition against the audited certificate of N - 1, analytic closures listed separately
                       and recomputed, closure kind, Gurobi status, ObjBound and cutoff of every record that stores
                       them, no NaN or infinite bound on a closed record
  ns_audit.py          audits non-spinal exclusion records: expected branches (and t0 sub-branches) present and
                       closed, no unexpected branch, solver threshold <= gamma_ns(N), condition for n_2 = 1 against
                       the audited results of N - 1, analytic statuses recomputed, Gurobi status, ObjBound and
                       cutoff of every record that stores them
  test_auditors.py     negative tests of the two auditors: 25 tampered copies of archived records, each rejected
  chain_summary.py     summary of the audited chain N = 5..16 (exit 1 unless every link holds)
  make_index.py        writes results/INDEX.csv (one row per archived run)
  exact_witness.py     turns a solver witness into an exactly balanced rational stack (--check: verify only)
  ns_check.py          non-spinal exclusion branches S_ns(n, l*, gamma_ns); --split-t0 for the t0 split
  pd_race.py           re-solves open profiles of a run in a separate process (more threads / other seed)
  pd_split.py          rightmost-layer split of one profile (independent cross-check at N = 13)
  verify_rational_stacks.py
                       exact rational construction and zero-residual verification of the stacks N = 4..16
  check_witness.py     re-checks a feasible stack reported by a solver run (forces, moments, contacts, spinality)
  check_exact_stacks.py
                       independent exact re-check of the stacks, distances to S*(N), S2+ and non-spinal chains
  summarize_certify.py collects result JSON files of monolithic runs into one table
  profiles.py          solver-free helpers (profile enumeration, energy test, loaded-spine bound, bound tables)
  exact_spines/        first verifier of the exact spinal optima (interval arithmetic, standard library; README inside)
  check_algebraic_spines.py
                       second verifier: the stacks of Appendix C in the field Q(x) (SymPy); signs by interval
                       enclosure on the isolating interval, else by exact Sturm counting
  spine_faces.py       second verifier: upper bounds at the exception sizes N = 5, 7 (fixed-load bound, box covering)
scripts/             one shell script per experiment (see the table below)
data/                rational_stacks.json: exact fractions of all blocks and contact forces, N = 4..16
results/             archived results of the runs reported in the paper (see results/README.md)
```

## Requirements

* Python 3.10 (the paper used 3.10.18), `numpy`, `scipy`, `mpmath` (see `requirements.txt` or
  `environment.yml`).
* **Gurobi** with `gurobipy` (the paper used Gurobi 13.0.1). A license is required:
  free academic licenses are available from Gurobi. Without a license, the `gurobipy` package from PyPI runs
  in a size-limited mode that is sufficient for the certificates and the non-spinal checks for $`N \le 9`$
  (tested with `gurobipy` 14.0.0); larger sizes need a full license.
* The exact-arithmetic checks (`scripts/verify_exact.sh`) and the auditors (`scripts/audit_archive.sh`,
  `src/pd_certificate.py`, `src/ns_audit.py`) need **no** optimization solver: only Python with `mpmath`.
  The verifiers of the exact spinal optima (`scripts/exact_spines.sh`) also need `sympy`.

```bash
conda env create -f environment.yml && conda activate maxoverhang
# or: pip install -r requirements.txt
```

## Quick start

```bash
scripts/verify_exact.sh          # exact stacks N = 4..16, S2+ chain and non-spinal chain (seconds, no solver)
scripts/exact_spines.sh          # exact spinal optima N = 4..19, Appendix C (seconds, no solver; needs sympy)
scripts/anchor_n4.sh             # N = 4 anchor with M0 at gap 1e-9
scripts/certify.sh 9 4           # tight certificate for N = 9 with 4 parallel workers (seconds)
scripts/nonspinal.sh 9 4         # non-spinal exclusion for N = 9
```

New outputs are written to `runs/` (change with `OUT=/path scripts/...`); the archived results in
`results/` are never overwritten. All scripts use `python3` (change with `PYTHON=...`).

## Paper ↔ scripts ↔ archived results

| Paper item | Script | Archived result |
|---|---|---|
| Exact rational stacks, Appendix A | `scripts/verify_exact.sh` | `data/rational_stacks.json` |
| S2+ chain table; non-spinal chain margins (Section 6.1) | `scripts/verify_exact.sh` | — (exact computation) |
| Exact spinal optima, Theorem 30(a) and Appendix C | `scripts/exact_spines.sh` (no solver) | `results/exact_spines/` |
| Anchor $`\overline{z}(4) = 1.1678933`$ | `scripts/anchor_n4.sh` | `results/anchor_n4/` |
| Main certificate table, tight certificates $`N = 5..16`$ | `scripts/certify.sh N [workers]` | `results/certificates/` |
| Main certificate table, non-spinal exclusion $`N = 5..16`$ | `scripts/nonspinal.sh N [workers]` | `results/nonspinal/` |
| $`N = 16`$ sub-branch table (Appendix B) | `scripts/nonspinal_n16_split.sh` | `results/nonspinal/` |
| Re-audit of all archived certificates and exclusions, including the S2+ and non-spinal conditions between consecutive sizes, the analytic closures and the chain summary | `scripts/audit_archive.sh` (no solver) | `results/` |
| Negative tests of the auditors: tampered copies of archived records (status, bound, closure kind, Gurobi status, ObjBound, cutoff, previous certificate) must be rejected | `python3 src/test_auditors.py` (no solver; also run by `scripts/audit_archive.sh`) | — |
| Index of every archived run (solver version, options, parameters, cutoff or threshold, files) | `python3 src/make_index.py` | `results/INDEX.csv` |
| Positive reverse controls with witnesses, and the exactly balanced rational non-spinal stack at $`N = 9`$ (Section 6) | `scripts/reverse_controls.sh`; `src/exact_witness.py --check` (no solver) | `results/reverse_controls/` |
| Stability under tighter tolerances (Section 6.1): hardest certificate profiles at $`N = 14, 15`$; $`N = 16`$ non-spinal sub-branches | `scripts/tolerance_check.sh`; `scripts/tolerance_check_ns16.sh` | `results/tolerance_check/` |
| Raw Gurobi logs and model files (MPS) of key instances: certificates and non-spinal exclusions at $`N = 5, 7`$, and the 11 $`N = 16`$ sub-branches that need the solver; every record also stores the closure kind (analytic or solver), the raw Gurobi status, ObjBound and the cutoff | `scripts/raw_logs.sh [workers]` | `results/raw_logs/` |
| Matched comparison with and without spinal cuts ($`N = 11, 12, 13`$, Section 6.2) | `scripts/ablation_sp.sh [workers]` | `results/ablation_sp/` |
| Methodological-progression table ($`N = 12, 13`$) | `scripts/baselines.sh {m0_12, m0_13, pd_12, pd_13, pd_13_race, pd_13_split, pd_13_cert}` | `results/baselines/` |

### Reference effort (single-thread subproblems, Intel Core Ultra 7 265K)

| | $`N \le 12`$ | $`N = 13`$ | $`N = 14`$ | $`N = 15`$ | $`N = 16`$ |
|---|---|---|---|---|---|
| tight certificate (thread-s) | $`\le 37`$ | 247 | 1,279 | 8,641 | 72,855 |
| non-spinal exclusion (thread-s) | $`\le 46`$ | 320 | 1,344 | 16,499 | 96,756 + 16,244 (t0 split) |

The monolithic baselines of the progression table run for 48–72 hours each.

## Settings used in the paper

* Tight cutoff $`\gamma_N = \underline{z}_{\mathbb{Q}}(N) + 10^{-5}`$, rounded upward at the 7th decimal
  (`CUT` in `scripts/common.sh`); at $`N = 15`$ the runs used the slightly larger cutoff 2.0929064 instead of 2.0929063.
* Certificate options (identical for every $`N`$):
  `profile_run.py N --model M1u --strengthen TB,S1,S2P,S4,C2 --rank-box --c2-layer --energy --spine
  --cutoff γ_N --param Heuristics=0.0 --param IntFeasTol=1e-9`; one thread per subproblem.
* Non-spinal exclusion: `ns_check.py N --rel 1e-2` (per-branch limit 7,200 s), the same options plus the
  branch constraints; at $`N = 16`$ the three branches left open are split by the lowest impure layer $`t_0`$.

Option names follow the paper: `TB` coordinate and contact boxes, `S1` table-edge equality, `S2P` base-singleton
reduction S2+ ($`n_2 = 2`$), `S4` base right-corner equality, `C2` / `--c2-layer` cut-set cuts C2 and
C2-Layer, `--rank-box` relative-rank boxes, `--energy` energy filter, `--spine` spinal cuts (family SP).

## Notes

* **Base-singleton precondition.** `certify_run.KNOWN_LB` and `certify_run.CERT_UB` hold the exact lower bounds
  $`\underline{z}_{\mathbb{Q}}(N)`$ (truncated) and the certified upper bounds $`\overline{z}(N)`$ of the paper. They
  are used only to check the S2+ precondition $`\overline{z}(N-1) + 1/(2N) < \underline{z}_{\mathbb{Q}}(N)`$, which
  never enters the optimization model. The archived runs checked the same precondition with earlier, valid but
  looser constants; the models, and hence all results, are identical.
* **Other Gurobi versions.** Node counts and run times depend on the Gurobi version and the machine; the
  certificates themselves (every profile and branch closed) do not. For example, with `gurobipy` 14.0.0 the
  certificates for $`N = 5..9`$ close with 12, 3, 109, 98 and 539 nodes (paper: 12, 3, 204, 85, 486), and the
  $`N = 7`$ reverse control recovers the same non-spinal stack of overhang 1.3842371. With default tolerances,
  Gurobi 14.0.0 reports 1.1678937 for the $`N = 4`$ anchor; adding
  `--param FeasibilityTol=1e-9 --param IntFeasTol=1e-9` gives 1.1678932, as with Gurobi 13.0.1.
* **Rounding.** Profiles are pruned by the energy test in exact arithmetic (every pruned profile has margin
  at most -1/4), and candidate layers are removed with margins of at least 0.0165 (`check_exact_stacks.py`,
  part 5). The loaded-spine constants are computed in double precision and increased by 1e-9; they differ from
  a 60-digit evaluation by at most 2e-15 for all $`\ell \le N \le 20`$ (1.2e-15 with Python 3.10; Python 3.12 and later
  sum floating-point numbers more accurately and give 4.7e-16; the archived runs used Python 3.10), and for $`N \le 16`$ the optimal loads of the
  maximizing height have increments of at least 1.19, so $`\hat{S}(N) = S^*(N)`$ (part 6). The non-spinal threshold is $`\gamma_{\mathrm{ns}}(N)`$, the value $`0.99\,S^*(N)`$ rounded upward at the
  7th decimal; the double-precision threshold used by the solver does not exceed it, and every closed record
  stores the threshold of its run (`ns_audit.py`). The archived `ns_summary.txt` prints the solver threshold
  rounded to the nearest 7th decimal, so it can show 1e-7 less than $`\gamma_{\mathrm{ns}}(N)`$.
* **Exact arithmetic.** The lower bounds of the paper come from `verify_rational_stacks.py`, which uses Python
  `fractions.Fraction` only (no floating point); `check_exact_stacks.py` re-checks every stack with a separate
  equilibrium checker.

## License

Code (`src/`, `scripts/`): MIT License, see `LICENSE`.
Data (`data/`, `results/`): CC BY 4.0, see `LICENSE-DATA`.

## Citation

Please cite the paper and this archive; see `CITATION.cff`.
Archived on Zenodo; all versions: [10.5281/zenodo.24518950](https://doi.org/10.5281/zenodo.24518950). The paper uses
version 1.0.4 (this version). Version 1.0.0, [10.5281/zenodo.24518951](https://doi.org/10.5281/zenodo.24518951), is the
first release; version 1.0.1, [10.5281/zenodo.24546538](https://doi.org/10.5281/zenodo.24546538), adds the solver-free
auditing tools and the records of the matched comparison, the tolerance check and the reverse controls; version 1.0.2,
[10.5281/zenodo.24554261](https://doi.org/10.5281/zenodo.24554261), adds the two verifiers of the exact spinal optima
(Appendix C) with their outputs, the exactly balanced rational non-spinal stack at $`N = 9`$, the tolerance check of the
$`N = 16`$ non-spinal sub-branches, and the audit of the conditions that link consecutive sizes; version 1.0.3,
[10.5281/zenodo.24559614](https://doi.org/10.5281/zenodo.24559614), adds the closure kind, raw status and ObjBound in
every new record, the mechanical audit of all closures without a solver run, and a directed interval check of the
constants of the spinal cuts; version 1.0.4 adds the raw Gurobi logs and model files of key instances
(`results/raw_logs/`), the audit of the raw status, ObjBound and cutoff of every record that stores them, stricter
checks of bounds, closure kinds and previous certificates in both auditors, and negative tests of the auditors
(`src/test_auditors.py`).
