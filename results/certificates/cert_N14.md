# Layer-profile certificate, N = 14

- certified: **True**
- feasible value (lower bound): 2.0382169
- upper bound (cutoff): 2.03822700  (relative gap 4.96e-06)
- profiles: 403 / 403 closed; status counts {'CUTOFF': 367, 'INFEASIBLE': 36}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 2.038227 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(13) + 1/(2N) < lower bound: margin 0.01798 (holds); zbar(13) = 1.984519 from cert_N13.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 25 by the energy test (largest Delta E -0.25), 1 with no candidate layer; 377 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 403 records, 0 mismatches; conflicting repeated records: 0
- total nodes 3,866,654; total thread-seconds 1,279

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2-3-3-1 | CUTOFF | 457,900 | 144.4 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-2-2-2-2-2-1 | CUTOFF | 434,147 | 120.8 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-2-3-3-2-1 | CUTOFF | 285,191 | 110.2 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-2-2-2-3-2 | CUTOFF | 307,050 | 80.3 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-3-2-2-3-1 | CUTOFF | 280,573 | 79.5 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-3-3-2-2-1 | CUTOFF | 196,186 | 70.9 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-2-2-3-2-1-1 | CUTOFF | 217,471 | 64.6 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
| 1-2-2-2-3-2-2 | CUTOFF | 194,045 | 52.3 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT2.0382270-C2L-SP1-PDcert_N14_s42_profiles.csv |
