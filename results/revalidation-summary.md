# Revalidation summary — 2026-11-23

Local source commit: a48377c055ce4929d67c4d996ab16f844458eb7a

| Check | Result |
|---|---|
| Unit tests | 4/4 PASS |
| Symbolic projector identities | PASS |
| Numeric projector identity regression | PASS at 100 / 1,000 / 10,000 samples |
| Maximum residual at 10,000 samples | 1.19e-12 |
| Acceptance threshold | 5e-11 |
| Local cubic W-Goldstone-Higgs benchmark | PASS |
| Off-shell pole-residue ambiguity | PASS |

Interpretation: these checks validate algebraic identities and one local vertex contraction. They do not constitute a complete five-point amplitude calculation, a five-point physical ST validation, or a physical phase-space convergence study.

Native outputs and detailed provenance are in the local snapshot under results/, environment/, and provenance/.
