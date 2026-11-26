# Layer-profile certificate, N = 11

- certified: **True**
- feasible value (lower bound): 1.8587795
- upper bound (cutoff): 1.85878960  (relative gap 5.43e-06)
- profiles: 70 / 70 closed; status counts {'CUTOFF': 60, 'INFEASIBLE': 10}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.8587896 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(10) + 1/(2N) < lower bound: margin 0.02619 (holds); zbar(10) = 1.7871359 from cert_N10.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 5 by the energy test (largest Delta E -0.5), 0 with no candidate layer; 65 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 70 records, 0 mismatches; conflicting repeated records: 0
- total nodes 22,191; total thread-seconds 8

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-3-2-1 | CUTOFF | 3,976 | 1.4 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-2-2-3-1 | CUTOFF | 4,955 | 1.2 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-2-2-2-2 | CUTOFF | 4,053 | 1.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-3-2-2-1 | CUTOFF | 2,178 | 0.5 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-2-2-2-1-1 | CUTOFF | 2,534 | 0.5 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-3-3-2 | CUTOFF | 299 | 0.2 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-2-3-3 | CUTOFF | 443 | 0.2 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
| 1-2-1-2-2-2-1 | CUTOFF | 443 | 0.1 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.8587896-C2L-SP1-PDcert_N11_s42_profiles.csv |
