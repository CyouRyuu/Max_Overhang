# Layer-profile certificate, N = 8

- certified: **True**
- feasible value (lower bound): 1.6315068
- upper bound (cutoff): 1.63151690  (relative gap 6.19e-06)
- profiles: 12 / 12 closed; status counts {'CUTOFF': 11, 'INFEASIBLE': 1}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.6315169 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(7) + 1/(2N) < lower bound: margin 0.03895 (holds); zbar(7) = 1.5300553 from cert_N7.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 1 by the energy test (largest Delta E -2.25), 0 with no candidate layer; 11 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 12 records, 0 mismatches; conflicting repeated records: 0
- bounds that are NaN, or not finite on a CUTOFF or OPTIMAL record: 0; Gurobi status, objective bound and cutoff of the records not recorded by this run
- total nodes 85; total thread-seconds 0

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2-1 | CUTOFF | 60 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-1-2-2 | CUTOFF | 8 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-1-2-1-1 | CUTOFF | 3 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-2-1-1-1 | CUTOFF | 4 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-2-3 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-2-1-2 | CUTOFF | 5 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-1-1-2-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
| 1-2-3-2 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.6315169-C2L-SP1-PDcert_N8_s42_profiles.csv |
