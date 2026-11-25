# results/

Archived records of the runs reported in the paper (Gurobi 13.0.1, Python 3.10.18, Intel Core Ultra 7 265K,
50 GB, Linux). New runs of the scripts go to `runs/`, never here.

| directory | content | paper |
|---|---|---|
| `anchor_n4/` | monolithic M0 at N = 4, relative gap 1e-9 | anchor zbar(4) = 1.1678933 |
| `certificates/` | per N = 5..16: run summary (`*_result.json`), per-profile records (`*_profiles.csv`), assembled certificate (`cert_N*.json`, readable `cert_N*.md`) | main certificate table (tight certificates) |
| `nonspinal/` | `ns_summary.txt` and per-branch records `ns_N*_rel0.01.csv`; `ns_N16_rel0.01_only_t0.csv` holds the 14 sub-branches of the t0 split at N = 16 | main certificate table (non-spinal), Appendix B |
| `baselines/m0/` | monolithic M0 at N = 12 (72 h) and N = 13 (48 h): result, best incumbent, console log | progression table |
| `baselines/nosp_N12/` | layer-profile decomposition without spinal cuts, N = 12 | progression table |
| `baselines/nosp_N13/` | the same at N = 13: profile records, the hard profile (1,2,2,2,3,2,1) solved on 3 threads (`race/`), its rightmost-layer split (`split_1-2-2-2-3-2-1/`), and the assembled certificate | progression table and cross-check |

Notes

* File names carry the run options: `M1u` model, `TB-S1-S2P-S4-C2` structural reductions, `RB` rank boxes,
  `EN` energy filter, `SP` spinal cuts, `CUT<value>` the cutoff, `C2L` C2-Layer cuts, `s42` the seed.
  The labels `SP1` and `PD6` only name the parameter set (`Heuristics=0`, plus `IntFeasTol=1e-9` for `SP1`);
  `scripts/baselines.sh` uses the label `NOSP` for the runs without spinal cuts.
* The runs without spinal cuts used the cutoff LB (1 + 5e-4) with the feasible values available at the time,
  i.e. 1.9260577 at N = 12 and 1.9855056 at N = 13. The latter LB, 1.9845133, is a solver incumbent of the
  monolithic N = 13 run; it exceeds the exact optimum 1.9845089 by 4.4e-6 because of the solver's feasibility
  tolerances. It only sets that cutoff and plays no role in the tight certificates.
* The certificate files `cert_N*.json/md` were assembled by `src/pd_certificate.py` from the archived
  `*_profiles.csv` with the constants of the paper (exact z_Q(N), tight zbar(N - 1)), so they show the S2+
  margins of the paper's chain table; the profile records themselves are the original solver output.
* `nodes` and `runtime_sec` are per subproblem; thread-seconds in the paper are sums of `runtime_sec`
  over single-thread subproblems.
