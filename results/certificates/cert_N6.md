# Layer-profile certificate, N = 6

- certified: **True**
- feasible value (lower bound): 1.4367002
- upper bound (cutoff): 1.43671030  (relative gap 7.03e-06)
- profiles: 4 / 4 closed; status counts {'CUTOFF': 3, 'INFEASIBLE': 1}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.4367103 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(5) + 1/(2N) < lower bound: margin 0.04880 (holds); zbar(5) = 1.3045649 from cert_N5.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 1 by the energy test (largest Delta E -0.5), 0 with no candidate layer; 3 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 4 records, 0 mismatches; conflicting repeated records: 0
- bounds that are NaN, or not finite on a CUTOFF or OPTIMAL record: 0; Gurobi status, objective bound and cutoff of the records not recorded by this run
- total nodes 3; total thread-seconds 0

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-3 | INFEASIBLE | 0 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.4367103-C2L-SP1-PDcert_N6_s42_profiles.csv |
| 1-2-2-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.4367103-C2L-SP1-PDcert_N6_s42_profiles.csv |
| 1-2-1-1-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.4367103-C2L-SP1-PDcert_N6_s42_profiles.csv |
| 1-2-1-2 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.4367103-C2L-SP1-PDcert_N6_s42_profiles.csv |
