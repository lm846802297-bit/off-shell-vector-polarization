# Research status — initial theory snapshot

Date: 2026-11-23

## Completed in the local snapshot

1. Fixed the theory target to minimal spontaneously broken SU(2) Yang–Mills theory with one complex scalar doublet and an optional vectorlike fermion doublet.
2. Recorded a tree-level BRST/Slavnov–Taylor benchmark for one off-shell vector leg and four on-shell physical external legs.
3. Derived the off-shell projector contractions in terms of \(\Delta=q^2-m^2\) and verified them symbolically.
4. Verified a local W–Goldstone–Higgs vertex contraction identity against a conserved attached current.
5. Ran deterministic floating-point identity scans at 100, 1,000, and 10,000 points; the maximum residual at 10,000 points was \(1.19\times10^{-12}\), below the \(5\times10^{-11}\) threshold.
6. Preserved an initial mixed-index implementation failure and its source snapshot.
7. Verified the local off-shell ambiguity and pole-residue invariance under \(J^\mu\to J^\mu+(q^2-m^2)K^\mu\).
8. Formulated a conditional five-point recursion target and identified Goldstone-sector factorization closure and local boundary/contact terms as proof obstructions.

## Not completed

- No complete five-point Feynman amplitude has been generated and summed.
- The allowed local \(K^\mu\) basis for the selected five-point current has not been classified.
- No five-point recursion residues have been computed explicitly.
- No five-point phase-space scan or mass-scaling fit has been performed. The 100/1,000/10,000-point scan is an algebraic identity regression, not a physical phase-space calculation.
- No loop-level or complex-pole result has been established.
- The remote GitHub branch currently contains this status document and repository scaffold only; the full local research snapshot still needs to be synchronized.

## Acceptance rule

A result is promoted from `open` to `derived` only after its assumptions and algebra are explicit. It is promoted to `tested` only after an executable test and native output exist. It is promoted to `established` only after an independent physical check (for example, an ST identity or a separate amplitude implementation) agrees.
