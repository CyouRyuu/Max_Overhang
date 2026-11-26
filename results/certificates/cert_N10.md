# Layer-profile certificate, N = 10

- certified: **True**
- feasible value (lower bound): 1.7871258
- upper bound (cutoff): 1.78713590  (relative gap 5.65e-06)
- profiles: 39 / 39 closed; status counts {'CUTOFF': 32, 'INFEASIBLE': 7}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.7871359 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(9) + 1/(2N) < lower bound: margin 0.02184 (holds); zbar(9) = 1.7152829 from cert_N9.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 5 by the energy test (largest Delta E -0.5), 0 with no candidate layer; 34 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 39 records, 0 mismatches; conflicting repeated records: 0
- total nodes 5,229; total thread-seconds 2

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2-2-1 | CUTOFF | 3,691 | 0.8 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-2-3-2 | CUTOFF | 274 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-3-2-2 | CUTOFF | 142 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-2-2-3 | CUTOFF | 75 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-2-3-1-1 | CUTOFF | 140 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-3-2-1-1 | CUTOFF | 153 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-3-3-1 | CUTOFF | 1 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
| 1-2-1-2-2-2 | CUTOFF | 100 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.7871359-C2L-SP1-PDcert_N10_s42_profiles.csv |
