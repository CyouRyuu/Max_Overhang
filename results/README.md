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
| `reverse_controls/` | positive reverse controls (`scripts/reverse_controls.sh`, log `rc.log`): N = 13 spinal stack 1.9845088923 at cutoff 1.9844089 (`n13_spinal/`; `n13_certificate_settings.log` is the same control with the certificate settings, no incumbent in 1 h, dual bound 1.9845089); N = 7 at gamma 1.38 and N = 9 at gamma 1.60, open, with witnesses `*_witnesses.json` (overhangs 1.3842371 and 1.6318871, both non-spinal under the closed resting relation through a point contact, spinal in the sense of Paterson and Zwick); N = 8 at gamma 1.60, all 48 branches closed (`n8/`; it excludes N = 9 non-spinal stacks with n_2 = 1 above 1.60 + 1/18); N = 9 at rho = 3e-2, all 92 branches closed; `n9/ns_N9_exact_stack.json`: the N = 9 witness made exactly balanced in rational arithmetic (`src/exact_witness.py`), overhang 1.6318870324..., every loaded contact of positive length, non-spinal in both senses | Section 6 |
| `tolerance_check/N14, N15/` | the three most expensive profiles of the tight certificates at N = 14 and 15 re-solved with FeasibilityTol = OptimalityTol = IntFeasTol = 1e-9 and NumericFocus = 3 (label `TOL`): all six close by cutoff | Section 6.1 |
| `tolerance_check/ns16/` | the 14 t0 sub-branches of the N = 16 non-spinal exclusion (sub-branch table of Appendix B) re-solved with FeasibilityTol = OptimalityTol = 1e-9, IntFeasTol = 1e-9 and NumericFocus = 3 (`scripts/tolerance_check_ns16.sh 11`): 11 close by cutoff, 3 need no computation (t0 = l*); 56,439,244 nodes, 23,443 thread-s (16,244 at default settings); `ns_summary.txt`, per-sub-branch records, and the audit output `tol16.log` (Gurobi licence lines removed) ending in NON-SPINAL EXCLUSION CERTIFIED | Section 6.1 |
| `exact_spines/` | exact spinal optima (`scripts/exact_spines.sh`, no solver): the CSV files of the first verifier (`src/exact_spines/README.md` describes them) and the output of its initial script (`original_stdout.txt`); the leaves of the box covering at N = 5, 7 (`face_leaves/`, replayed by `src/check_face_leaves.py`); the logs of the second verifier (`spine_faces.log`, `check_algebraic_spines.log`) | exact spinal optima theorem, Appendix C |
| `raw_logs/` | rerun of key instances with the settings of the archived runs (`scripts/raw_logs.sh`): tight certificates and non-spinal exclusions at N = 5, 7, and the 14 t0 sub-branches of the N = 16 exclusion (11 solved). One Gurobi log and one model file (MPS, gzip) per solved subproblem or branch; every record also stores the closure kind (analytic or solver), the raw Gurobi status code, ObjBound and the cutoff. Node counts equal those of the archived runs (N = 5: 12, N = 7: 204, N = 16 sub-branches: 58,300,685). In the 11 N = 16 models the nonzero linear constraint coefficients lie in [0.973, 240], all quadratic coefficients equal 1, the nonzero variable bounds lie in [0.377, 35.2] and the nonzero right-hand sides in [0.5, 59.5] (the Gurobi log rounds these ranges to one significant figure: [1, 200], [0.4, 40], [0.5, 60]); no log reports numerical issues. `raw_logs_console.log` is the console output, including the audits | Section 5.1, Section 6.1, Appendix B |
| `screening/` | exact rational re-check of every analytic screening decision (`src/check_screening.py`): per N, one row per energy-pruned profile, per candidate layer dropped from a certificate subproblem (criterion base envelope, energy layer bound or pure-chain bound, exact margin, threshold) and per non-spinal branch with a bound conflict; `screening_summary.txt` gives the counts and splits the archived INFEASIBLE branch records without nodes into exact bound conflicts and presolve closures of the solver | Sections 5.1, 6.1 |
| `INDEX.csv` | one row per archived run: record, kind, N, Gurobi version, model options and parameters, cutoff or threshold, threads, time limit, status, bound, nodes, files (`src/make_index.py`) | all |

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
  paper (exact z_Q(N), tight zbar(N - 1)), so they show the S2+ margins of the paper's constants table; the profile,
  race, split and branch records themselves are the original solver output (`scripts/audit_archive.sh`).
* Console logs: the host name and home paths are replaced by `workstation` and `~`, and the Gurobi licence
  lines are removed.
* The `MIPGap` values in `INDEX.csv` (5e-4 for the profile runs, 1e-4 for the non-spinal branches) play no role in
  any certificate: the auditors accept a subproblem or branch only with status `CUTOFF` (no solution above the
  cutoff or threshold exists) or `INFEASIBLE`, or, for non-spinal branches, `EMPTY` or `TIE` (closed analytically without a solver run);
  a run that stops on the gap with a solution above the cutoff would end with status `OPTIMAL` and fail the audit.
* `nodes` and `runtime_sec` are per subproblem; thread-seconds in the paper are sums of `runtime_sec`
  over single-thread subproblems.

Fields stored by each archive

* `certificates/`, `nonspinal/`, `baselines/`, `ablation_sp/`, `tolerance_check/`, `reverse_controls/`: the status
  of each subproblem or branch as summarized by the run (`CUTOFF`, `INFEASIBLE`, `EMPTY`, `TIE`, open statuses),
  the recorded bound, nodes and run time, and in the run summaries the cutoff or threshold and all parameters.
  These archives predate the extra fields below; the closure kind of their analytic records is recomputed by the
  auditors and, exactly, by `src/check_screening.py`.
* `raw_logs/` (N = 5, 7 and the N = 16 t0 split, rerun with the archived settings): in addition the closure kind
  (`analytic: ...`, `bounds before solving`, `solver`), the raw Gurobi status code, `ObjBound` and the cutoff of
  each record, one Gurobi log and one model file (MPS) per subproblem solved, checked by the auditors.
* A record with status `INFEASIBLE` closed by the solver stores the bound `-inf`, the raw value of the run. Since the
  cutoff enters presolve and the non-spinal branches contain rows that hold only above the threshold, such a record
  proves only that no point of the subproblem has objective above the cutoff; the paper, and the certificate
  reports, read every closed record as the bound given by the cutoff. Only the energy test proves a profile empty.
* No archived branch was closed by `bounds before solving`: `src/check_screening.py` finds no exact bound conflict,
  and every archived `INFEASIBLE` branch record without nodes is a presolve closure of the solver.
