# Layer-profile certificate, N = 12

- certified: **True**
- feasible value (lower bound): 1.9250913
- upper bound (cutoff): 1.92510140  (relative gap 5.25e-06)
- profiles: 126 / 126 closed; status counts {'CUTOFF': 106, 'INFEASIBLE': 20}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.9251014 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(11) + 1/(2N) < lower bound: margin 0.02464 (holds); zbar(11) = 1.8587896 from cert_N11.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 12 by the energy test (largest Delta E -0.25), 0 with no candidate layer; 114 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 126 records, 0 mismatches; conflicting repeated records: 0
- total nodes 125,352; total thread-seconds 37

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2-2-2-1 | CUTOFF | 33,098 | 7.8 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-3-3-2-1 | CUTOFF | 11,069 | 3.5 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-2-2-3-2 | CUTOFF | 7,492 | 2.5 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-3-2-2-2 | CUTOFF | 7,843 | 2.4 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-2-3-2-1-1 | CUTOFF | 9,867 | 2.3 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-2-3-2-2 | CUTOFF | 8,362 | 2.3 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-2-3-3-1 | CUTOFF | 6,458 | 2.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
| 1-2-3-2-3-1 | CUTOFF | 6,282 | 1.9 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.9251014-C2L-SP1-PDcert_N12_s42_profiles.csv |
