# Changelog

All versions are archived on Zenodo under the concept DOI [10.5281/zenodo.24518950](https://doi.org/10.5281/zenodo.24518950).

## 1.0.6 (2026-11-27)
The version used by the paper. Exact rational re-check of every analytic screening decision (`src/check_screening.py`,
`results/screening/`); leaves of the box covering at the exception sizes $`N = 5, 7`$ with an independent replay
(`src/spine_faces.py --leaves`, `src/check_face_leaves.py`, `results/exact_spines/face_leaves/`); documentation of
the fields stored by each archive; title of the paper updated; shorter citation section.

## 1.0.5 (2026-11-26), [10.5281/zenodo.24563697](https://doi.org/10.5281/zenodo.24563697)
No code or record changes; refers to the paper by section and table names instead of theorem and table numbers.

## 1.0.4 (2026-11-26), [10.5281/zenodo.24561497](https://doi.org/10.5281/zenodo.24561497)
Raw Gurobi logs and model files of key instances (`results/raw_logs/`); audit of the raw status, ObjBound and cutoff
of every record that stores them; stricter checks of bounds, closure kinds and previous certificates in both
auditors; negative tests of the auditors (`src/test_auditors.py`).

## 1.0.3 (2026-11-26), [10.5281/zenodo.24559614](https://doi.org/10.5281/zenodo.24559614)
Closure kind, raw status and ObjBound in every new record; mechanical audit of all closures without a solver run;
directed interval check of the constants of the spinal cuts.

## 1.0.2, [10.5281/zenodo.24554261](https://doi.org/10.5281/zenodo.24554261)
Two verifiers of the exact spinal optima with their outputs; exactly balanced rational non-spinal stack at $`N = 9`$;
tolerance check of the $`N = 16`$ non-spinal sub-branches; audit of the conditions that link consecutive sizes.

## 1.0.1, [10.5281/zenodo.24546538](https://doi.org/10.5281/zenodo.24546538)
Solver-free auditing tools; records of the matched comparison, the tolerance check and the reverse controls.

## 1.0.0, [10.5281/zenodo.24518951](https://doi.org/10.5281/zenodo.24518951)
First release.
