# PR: Initial theory and reproducibility snapshot

## Summary

- Establish the minimal spontaneously broken SU(2) Yang–Mills–Higgs model and conventions.
- Review Ke, Liu, and Mahbub, arXiv:2611.13096v1.
- Record the tree-level BRST/ST benchmark for an off-shell vector current with four on-shell physical legs.
- Prove algebraically that \(J^\mu\to J^\mu+(q^2-m^2)K^\mu\) preserves the pole residue while shifting regular/contact terms.
- Formulate the five-point on-shell-recursion proof target and identify its main obstructions.
- Add symbolic and independent numerical checks for the off-shell projector identities.
- Preserve the initial mixed-index implementation failure and corrected source snapshot.

## Local validation record

The local snapshot records the following commands as passing:

- `python scripts/verify_projector_identities.py`
- `python scripts/numeric_projector_scan.py` at 100, 1,000, and 10,000 points
- `python scripts/verify_vertex_st_identity.py`
- `python scripts/verify_offshell_ambiguity.py`
- `python -m unittest discover -s tests -v` (4 tests)

These results are reported from the local snapshot. They have not yet been independently rerun on GitHub Actions.

## Scope and limitations

The complete five-point Feynman amplitude and independent on-shell-recursive derivation are not claimed as complete. The BRST/ST relation is a field-theory benchmark, not a new theorem.

## Reproducibility

The local snapshot contains runtime/hardware metadata, source and input hashes, run manifests, native outputs, and preserved failure records under `environment/`, `provenance/`, and `results/`.
