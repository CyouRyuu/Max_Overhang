# Layer-profile certificate, N = 7

- certified: **True**
- feasible value (lower bound): 1.5300452
- upper bound (cutoff): 1.53005530  (relative gap 6.60e-06)
- profiles: 7 / 7 closed; status counts {'CUTOFF': 6, 'INFEASIBLE': 1}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.5300553 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(6) + 1/(2N) < lower bound: margin 0.02191 (holds); zbar(6) = 1.4367103 from cert_N6.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 1 by the energy test (largest Delta E -3.0), 0 with no candidate layer; 6 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 7 records, 0 mismatches; conflicting repeated records: 0
- bounds that are NaN, or not finite on a CUTOFF or OPTIMAL record: 0; Gurobi status, objective bound and cutoff of the records not recorded by this run
- total nodes 204; total thread-seconds 0

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2 | CUTOFF | 198 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-2-1-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-1-2-1 | CUTOFF | 3 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-4 | INFEASIBLE | 0 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-3-1 | CUTOFF | 0 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-1-1-2 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
| 1-2-1-1-1-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.5300553-C2L-SP1-PDcert_N7_s42_profiles.csv |
