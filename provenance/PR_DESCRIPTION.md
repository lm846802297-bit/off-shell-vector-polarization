# Draft PR: Research scaffold and reproducibility status

## Scope

This draft establishes the research question, minimal SU(2)-Higgs model, literature map, and acceptance protocol for the study of off-shell vector currents and unstable-vector spin decomposition.

## Local research snapshot

The local snapshot at commit a48377c055ce4929d67c4d996ab16f844458eb7a contains 48 tracked files, including derivations, scripts, tests, runtime records, native outputs, and a corrected SHA-256 manifest.

Local validation on 2026-11-23:
- 4/4 unit tests passed.
- Symbolic projector identities passed.
- Numerical projector identity regression passed at 100, 1,000, and 10,000 samples; threshold 5e-11.
- Local cubic W-Goldstone-Higgs vertex benchmark passed.
- Off-shell pole-residue ambiguity regression passed.

The numeric scan is an algebraic identity regression, not a physical phase-space convergence study. The cubic-vertex benchmark is not a complete five-point ST validation.

## Explicitly open

- Complete five-point Feynman amplitude and all factorization residues.
- Classification of allowed local K^mu shifts.
- Independent on-shell-recursive derivation of the five-point current relation.
- Loop-level and complex-pole extension.
- Independent GitHub Actions rerun.

## Remote synchronization status

This PR is intentionally Draft. The remote branch currently contains the scaffold and review metadata, not the full 48-file local source tree. Do not merge this PR as a complete reproducibility release. The full local Git bundle is available in the originating ChatGPT conversation; its SHA-256 is recorded in provenance/remote-sync.md.
