# Exact spinal optima: first verifier (Appendix C)

Interval-arithmetic verifier of the exact spinal optima of Appendix C, at stability margin $`\eta = 0`$.
It uses only the Python standard library; `symbolic_check.py` additionally needs SymPy.
The second, independent verifier is `src/check_algebraic_spines.py` together with `src/spine_faces.py`.
`scripts/exact_spines.sh` runs both and compares the output with `results/exact_spines/`.

```bash
python3 src/exact_spines/overhang_spine_check.py                       # N = 4, 6, 8..16 (first version)
python3 src/exact_spines/audit_spines.py --out runs/exact_spines        # N = 4..19, margins, exception faces
python3 src/exact_spines/symbolic_check.py --out runs/exact_spines      # optional: identities in Q(x), polynomials
```

## What is checked

* **Loads and height.** For each height $`K`$, the loads $`t_s`$ are rational functions of $`t_1 = x`$ through the
  recurrence (50), and $`x`$ is the unique positive root of $`t_K(x) = N`$. The root is isolated in an interval of
  width $`3 \times 10^{-26}`$ whose endpoints are checked in exact rational arithmetic, and $`\hat S_K(N)`$ is enclosed
  for every $`K \le N`$. The optimal height is separated from all others, and all added loads are positive there.
* **Ordinary blocks.** The shielded construction of Section 5.2 (equations (60) to (65), the seeds of Table A1,
  and for $`N = 17, 18, 19`$ the seeds of Appendix C) with these exact loads. Every quantity is enclosed in an
  interval with rational endpoints rounded outward to multiples of $`10^{-60}`$; strict inequalities are proved by
  a positive lower end, weak ones by a nonnegative lower end. Equalities follow from the construction and are
  checked separately in $`\mathbb{Q}(x)`$ by `symbolic_check.py`.
* **Exception sizes.** $`N = 5`$ uses equation (66) with $`t_2 = \sqrt{15/2}`$; $`N = 7`$ uses height 3 with loads
  $`(0, x, 2x^2, 7)`$, $`8x^4 - 14x - 7 = 0`$. The faces with a zero added load are covered by rational boxes
  (`exception_faces.csv`).

`overhang_spine_check.py` is kept unchanged (SHA-256
`ad6d1b0a9e5e2ff288b7ae2b94d1774fa0fc53fcac40980aa11c47c7c50002e6`); its output is
`results/exact_spines/original_stdout.txt`. `audit_spines.py` extends it to all $`N = 4, \dots, 19`$, the margins,
and the exception faces.

## Output files

| file | content |
|---|---|
| `heights.csv` | for every $`N`$ and height $`K`$: isolating interval of $`t_1`$, sign margins at its ends, enclosure of $`\hat S_K(N)`$ |
| `optimal_loads.csv`, `loads.csv` | every load $`t_s`$ of the loaded-spine optimum and of the construction with ordinary blocks |
| `blocks.csv`, `contacts.csv` | layer and left edge of every block; force and point of every loaded contact |
| `margins.csv` | every contact, same-layer pair, principal-block comparison, and excluded lower support |
| `summary.csv` | per $`N`$: smallest margin of each kind, the item attaining it, height separation |
| `exception_faces.csv` | rational upper bounds and box counts on the faces with a zero added load, $`N = 5, 7`$ |
| `polynomials.csv`, `identities.csv` | (`symbolic_check.py`) irreducible integer polynomials of the roots; per-block balance residuals in $`\mathbb{Q}(x)`$, all zero |

Block names: B spine, H shield, G bridge, T tower, C counterweight ($`N = 5`$). Fractions are written with
positive denominators; an empty field means that the item does not apply.

## Provenance

The scripts were drafted with the assistance of ChatGPT Astra (OpenAI) for an independent check of the manuscript,
then reproduced and checked by the authors, and are distributed under the MIT license of this repository
(see the declaration on the use of artificial intelligence tools in the paper).
