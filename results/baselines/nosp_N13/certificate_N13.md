# Layer-profile certificate, N = 13

- certified: **True**
- feasible value (lower bound): 1.9845089
- upper bound (cutoff): 1.98550556  (relative gap 5.02e-04)
- profiles: 225 / 225 closed; status counts {'CUTOFF': 220, 'INFEASIBLE': 5}
- missing 0, unexpected 0, repeated records 0 (last record used)
- cutoff of the run: record missing (accepted with --allow-missing-run-record)
- closed records with bound above the cutoff: 0; split profiles with missing/open layers: 0
- S2+ (n_2 = 2) condition zbar(12) + 1/(2N) < lower bound: margin 0.02095 (holds); zbar(12) = 1.9251014 from cert_N12.json, certified: True
- S2+ precondition checked against the audited certificate of N - 1
- closures without a solver run: 0 by the energy test (largest Delta E None), 0 with no candidate layer; 225 closed by the solver
- records without a solver run checked against the recomputed analytic reason: 224 records, 0 mismatches; conflicting repeated records: 0
- total nodes 57,847,794; total thread-seconds 13,266

| slowest profile | status | nodes | seconds | source |
|---|---|---|---|---|
| 1-2-2-2-2-2-1-1 | CUTOFF | 12,994,300 | 2534.6 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-2-2-3-2-1 | CUTOFF | 29,361,715 | 2385.6 | results/baselines/nosp_N13/race/P1-2-2-2-3-2-1_s1_result.json |
| 1-2-2-3-2-2-1 | CUTOFF | 8,580,812 | 1794.6 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-2-3-3-2 | CUTOFF | 3,680,545 | 895.2 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-2-2-2-2-2 | CUTOFF | 1,479,951 | 326.4 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-3-3-2-1-1 | CUTOFF | 253,179 | 72.8 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-2-2-2-3-1 | CUTOFF | 271,659 | 71.1 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
| 1-2-2-3-3-1-1 | CUTOFF | 182,874 | 48.9 | M1u-TB-S1-S2P-S4-C2-RB-C2L-PD6-PDcert_N13_s42_profiles.csv |
