# results/

Archived records of the runs reported in the paper (Gurobi 13.0.1, Python 3.10.18, Intel Core Ultra 7 265K,
50 GB, Linux). New runs of the scripts go to `runs/`, never here.

| directory | content | paper |
|---|---|---|
| `anchor_n4/` | monolithic M0 at N = 4, relative gap 1e-9 | anchor zbar(4) = 1.1678933 |
| `certificates/` | per N = 5..16: run summary (`*_result.json`), per-profile records (`*_profiles.csv`), assembled certificate (`cert_N*.json`, readable `cert_N*.md`) | main certificate table (tight certificates) |
| `nonspinal/` | `ns_summary.txt` and per-branch records `ns_N*_rel0.01.csv`; `ns_N16_rel0.01_only_t0.csv` holds the 14 sub-branches of the t0 split at N = 16 | main certificate table (non-spinal), Appendix B |
| `baselines/m0/` | monolithic M0 at N = 12 (72 h time limit) and N = 13 (stopped manually after 48 h, status INTERRUPTED): result, best incumbent, console log | progression table |
| `baselines/nosp_N12/` | layer-profile decomposition without spinal cuts, N = 12 | progression table |
| `baselines/nosp_N13/` | the same at N = 13: profile records, the hard profile (1,2,2,2,3,2,1) solved on 3 threads (`race/`), its rightmost-layer split (`split_1-2-2-2-3-2-1/`), and the assembled certificate; the main run was stopped manually after 9.2 h while this profile was still open on one thread, so its record is missing from the profile file and the certificate takes it from `race/` | progression table and cross-check |
| `ablation_sp/N11..N13/` | matched comparison without spinal cuts: the certificate options of `certificates/` minus family SP, same cutoff, 7 workers; time limit 3,600 s (N = 11, 12) and 14,400 s (N = 13) per profile; per-profile records, run summary, console log | matched comparison (Section 6.2); N = 11 is a complete SP-free certificate, N = 12 and 13 leave 2 and 4 profiles open |
| `reverse_controls/` | positive reverse controls (`scripts/reverse_controls.sh`, log `rc.log`): N = 13 spinal stack 1.9845088923 at cutoff 1.9844089 (`n13_spinal/`; `n13_certificate_settings.log` is the same control with the certificate settings, no incumbent in 1 h, dual bound 1.9845089); N = 7 at gamma 1.38 and N = 9 at gamma 1.60, open, with witnesses `*_witnesses.json` (overhangs 1.3842371 and 1.6318871, both non-spinal under the closed resting relation through a point contact, spinal in the sense of Paterson and Zwick); N = 8 at gamma 1.60, all 48 branches closed (`n8/`; it excludes N = 9 non-spinal stacks with n_2 = 1 above 1.60 + 1/18); N = 9 at rho = 3e-2, all 92 branches closed | Section 6 |
| `tolerance_check/N14, N15/` | the three most expensive profiles of the tight certificates at N = 14 and 15 re-solved with FeasibilityTol = OptimalityTol = IntFeasTol = 1e-9 and NumericFocus = 3 (label `TOL`): all six close by cutoff | Section 6.1 |

Notes

* File names carry the run options: `M1u` model, `TB-S1-S2P-S4-C2` structural reductions, `RB` rank boxes,
  `EN` energy filter, `SP` spinal cuts, `CUT<value>` the cutoff, `C2L` C2-Layer cuts, `s42` the seed.
  The labels `SP1` and `PD6` only name the parameter set (`Heuristics=0`, plus `IntFeasTol=1e-9` for `SP1`);
  `scripts/baselines.sh` uses the label `NOSP` for the runs without spinal cuts, `scripts/ablation_sp.sh` the
  label `NOSPABL` for the matched comparison (same parameters as `SP1`).
* The runs without spinal cuts used the cutoff LB (1 + 5e-4) with the feasible values available at the time,
  i.e. 1.9260577 at N = 12 and 1.9855056 at N = 13. The latter LB, 1.9845133, is a solver incumbent of the
  monolithic N = 13 run; it exceeds the exact optimum 1.9845089 by 4.4e-6 because of the solver's feasibility
  tolerances. It only sets that cutoff and plays no role in the tight certificates.
* The certificate reports (`certificates/cert_N*.json/md`, `baselines/nosp_N13/certificate_N13.*`) were
  regenerated from the archived records with the final auditor `src/pd_certificate.py` and the constants of the
  paper (exact z_Q(N), tight zbar(N - 1)), so they show the S2+ margins of the paper's chain table; the profile,
  race, split and branch records themselves are the original solver output (`scripts/audit_archive.sh`).
* Console logs: the host name and home paths are replaced by `workstation` and `~`, and the Gurobi licence
  lines are removed.
* `nodes` and `runtime_sec` are per subproblem; thread-seconds in the paper are sums of `runtime_sec`
  over single-thread subproblems.
