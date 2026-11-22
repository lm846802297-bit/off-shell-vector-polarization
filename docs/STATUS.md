# Research status — initial theory snapshot

Date: 2026-11-23

## Completed in the local snapshot

1. Fixed the initial theory target to the minimal SU(2) Yang-Mills theory with one complex scalar doublet and one optional vectorlike fermion doublet.
2. Recorded the tree-level BRST/Slavnov-Taylor benchmark for a current with one off-shell vector leg and four on-shell physical external legs, with the Goldstone phase convention matched to arXiv:2611.13096v1.
3. Derived the off-shell projector contractions in terms of Delta = q^2 - m^2 and verified them symbolically.
4. Verified a local W-Goldstone-Higgs vertex contraction identity against a conserved attached current.
5. Re-ran the deterministic floating-point identity scan at 100, 1,000, and 10,000 points. The maximum residual at 10,000 points is 1.19e-12, below the 5e-11 threshold.
6. Preserved the initial mixed-index implementation failure and its source snapshot; it was an implementation index-placement issue, not a failure of the analytic identities.
7. Verified the local off-shell ambiguity and pole-residue invariance with a separate symbolic regression.
8. Recorded a model-level field-redefinition realization of the ambiguity. The full R_xi/BRST completion remains open.
9. Added a conditional five-point recursion statement and isolated two proof obstructions: Goldstone-sector factorization closure and a possible local boundary polynomial.
10. Re-ran all four unit tests and all four verification scripts in the ChatGPT Linux environment on 2026-11-23.

## Not completed

- No complete five-point Feynman amplitude has been generated or summed.
- The allowed local K^mu basis for the selected five-point current has not been classified.
- No five-point recursion residues or physical phase-space scan have been computed.
- The local W-Goldstone-Higgs check is a cubic-vertex benchmark, not an independent validation of the complete five-point ST identity.
- No loop-level or complex-pole result has been established.
- GitHub Actions has not independently rerun the calculations.
- The public remote repository currently has only the scaffold and review metadata. The complete 48-file source snapshot and original Git history remain in the local downloadable bundle; they are not yet present as individual remote files.

## Revalidation record

- Runtime: environment/revalidation-runtime.json and environment/revalidation-runtime.txt.
- Native output: results/revalidation-2026-11-23.log.
- Structured run provenance: provenance/revalidation-run-manifest.json.
- SHA-256 manifest: environment/sha256-manifest.json.

## Acceptance rule

A result is promoted from open to derived only after assumptions and algebra are explicit. It is promoted to tested only after an executable test and native output exist. It is promoted to established only after an independent physical check (for example, a separate full-amplitude implementation or a complete ST identity check) agrees.
