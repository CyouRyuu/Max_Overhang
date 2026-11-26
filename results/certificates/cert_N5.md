# Layer-profile certificate, N = 5

- certified: **True**
- feasible value (lower bound): 1.3045548
- upper bound (cutoff): 1.30456490  (relative gap 7.74e-06)
- profiles: 2 / 2 closed; status counts {'CUTOFF': 2}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: 1.3045649 (<= certificate cutoff)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(4) + 1/(2N) < lower bound: margin 0.03666 (holds); zbar(4) = 1.1678933 from M0cert_N4_s42_result.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 0 by the energy test (largest Delta E None), 0 with no candidate layer; 2 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 2 records, 0 mismatches; conflicting repeated records: 0
- bounds that are NaN, or not finite on a CUTOFF or OPTIMAL record: 0; Gurobi status, objective bound and cutoff of the records not recorded by this run
- total nodes 12; total thread-seconds 0

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2 | CUTOFF | 11 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.3045649-C2L-SP1-PDcert_N5_s42_profiles.csv |
| 1-2-1-1 | CUTOFF | 1 | 0.0 | M1u-TB-S1-S2P-S4-C2-RB-EN-SP-CUT1.3045649-C2L-SP1-PDcert_N5_s42_profiles.csv |
