# Off-Shell Vector Polarization

Research repository on gauge-invariant spin decomposition of unstable vector resonances and off-shell currents in spontaneously broken gauge theories.

**Initial scope:** the tree-level five-point current in the minimal spontaneously broken SU(2) Yang-Mills-Higgs model, followed by complex-pole spin residues and polarization density matrices.

## Current research questions

1. Can the three-point off-shell current relation be extended to five-point and higher-point currents directly from on-shell constructibility?
2. What factorization, large-shift, and locality assumptions are needed to exclude a local boundary term?
3. How should spin-resolved residues be defined for mixed vector fields at a complex pole?
4. How should a physical pole density matrix be combined with a non-resonant off-shell current and its interference?

## Current status

- The off-shell projector identities from Ke, Liu, and Mahbub, arXiv:2611.13096v1, have been checked symbolically and numerically in the local research snapshot.
- The local ambiguity J^mu -> J^mu + (q^2-m^2)K^mu is derived: the pole residue is unchanged while regular/contact terms shift.
- A tree-level BRST/Slavnov-Taylor benchmark and a local cubic-vertex check are recorded.
- The complete five-point Feynman amplitude and independent on-shell-recursive derivation are **not yet complete**.
- No loop-level or complex-pole result is claimed.

## Reproducibility and governance

The local snapshot records runtime/hardware, source and input hashes, native outputs, and a corrected SHA-256 manifest. Its numerical scan is an algebraic identity regression, not a physical phase-space convergence study.

The main branch is reserved for reviewed research results. All new work enters main through a pull request. See STATUS, PR_DESCRIPTION, and provenance/remote-sync.md.

**Remote synchronization:** this public GitHub repository currently contains the research scaffold and review metadata. The complete 48-file local snapshot and its original Git history have not yet been imported into the remote tree. Do not merge the current draft PR as a complete source release.
