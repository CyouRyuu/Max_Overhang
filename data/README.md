# data/

`rational_stacks.json` — the optimal single-root stacks for N = 4..16 in exact rational arithmetic
(Appendix A of the paper), written by

    python3 src/verify_rational_stacks.py --dump-json data/rational_stacks.json

Conventions: block length 1 (half-length W = 1/2), block weight 1, stability margin eta = 0, table edge at x = 0.
All numbers are exact fractions stored as strings (`"p/q"`).

| field | meaning |
|---|---|
| `N` | number of blocks |
| `k` | parameter of the closed-form construction of Appendix A (`null` for the N = 5 stack, which has no bridge) |
| `profile` | layer profile (n_1, ..., n_H), bottom layer first |
| `overhang_frac`, `overhang_dec12` | exact overhang z_Q(N) and its 12-digit decimal value |
| `blocks` | per block: `id`, `layer`, `role` (spine / counterweight position), `X` = center coordinate |
| `contacts` | per contact (i on j; j = 0 is the table): resultant force `F` and moment `tau` = F * P about x = 0 |

Every stack satisfies force and moment equilibrium with zero residual, contact points inside both contact
intervals, non-overlap within layers, and a single root on the table; `src/check_exact_stacks.py` re-checks
this independently. License: CC BY 4.0 (see `../LICENSE-DATA`).
