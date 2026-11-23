# Off-Shell Vector Polarization — Repository Instructions

Shared research governance, compute execution, task handoff, scientific-change tracking, quality gates, result archival, structured records, document review and Git integration are defined by the central research infrastructure repository and pinned in [research/standards-lock.yaml](research/standards-lock.yaml). Start from the [Central Research Standards Index](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/CENTRAL_STANDARDS_INDEX.md).

This repository remains authoritative for the physical definition of off-shell currents, gauge/BRST conventions, Slavnov–Taylor identities, factorization and recursion assumptions, complex-pole residues, spin decomposition, polarization density matrices, and project-specific physical acceptance criteria.

## 1. Preflight and task handoff

Before substantial work:

1. Read this file, README.md, the remote synchronization record, and research/standards-lock.yaml.
2. Read the pinned Central Research Standards Index.
3. For a new derivation, amplitude calculation, symbolic/numerical run, parameter scan or multi-agent handoff, prepare a concise Research Task Brief using the [central Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [template](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/templates/research-task-brief.md).
4. Select the appropriate [central Research Quality Gate](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/RESEARCH_QUALITY_GATE.md) profile and define checks before the calculation.
5. Inspect the exact source revision, local snapshot/bundle, PRs, input manifests, native outputs and previous failure records before repeating work.

## 2. Project-specific theory conventions

- State metric signature, momentum-flow convention, Fourier phases, polarization normalization and all external-state normalization conventions.
- Distinguish on-shell physical polarization states from off-shell vector currents. Do not assign an arbitrary off-shell current a physical three-state spin space without a derived projection.
- State gauge fixing, Goldstone phase conventions, BRST transformations and the precise Slavnov–Taylor identity being tested.
- Distinguish the off-shell contracted current from an on-shell S-matrix element and from the residue of a physical complex pole.
- When using a current redefinition J^mu -> J^mu + (q^2-m^2)K^mu, track separately the pole residue and regular/contact contributions, including their interference in the full process.
- For unstable or mixed vector resonances, derive left/right pole eigenvectors, residue normalization and the dual polarization basis at the complex pole. Do not substitute real-mass on-shell completeness relations without proof.
- Define the polarization density matrix only for a physically specified pole residue or another explicitly defined state space. Establish gauge-parameter and polarization-basis transformation properties before claiming invariance or positivity.
- For on-shell recursion, state the complex shift, factorization channels, large-z behavior, locality/power-counting assumptions and possible boundary polynomial. A zero moving-pole sum does not determine an individual current if a boundary term remains.
- Keep tree-level results distinct from loop corrections, finite widths, complex-pole extensions and non-resonant off-shell contributions.

## 3. Computation, provenance and validation

Follow the pinned central [Research Compute Policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/RESEARCH_COMPUTE_POLICY.md). For every substantive calculation preserve:

1. Exact source commit SHA and clean/dirty tree state.
2. SHA-256 hashes of source, inputs, fixtures and configuration files.
3. Actual runtime, OS, hardware, interpreter/compiler and dependency versions.
4. Exact commands, model/gauge parameters, kinematics, complex shifts and numerical precision.
5. Native symbolic/numerical outputs and logs.
6. Numerical precision, convergence and uncertainty analysis appropriate to the physical claim.
7. Independent physical checks, including the relevant Ward/Slavnov–Taylor identity, factorization residues or an independent amplitude construction.
8. A run manifest connecting the actual executed source, inputs, environment, outputs and checks.

New substantial runs use central RCRPS v1.0.1, including the complete 14-field execution.source_control object. Historical records remain unchanged and must not be retrofitted with invented provenance.

## 4. Remote synchronization and archival boundary

The current public remote repository is not yet a complete reproducibility release. PR #1 records that the authoritative local snapshot has 48 tracked files at local commit `a48377c055ce4929d67c4d996ab16f844458eb7a`, while the remote research branch is only a partial scaffold.

- Treat the local Git bundle and its original history as authoritative until full import and verification are completed.
- Do not infer content identity from filenames, file counts, or a PR description.
- Before claiming synchronization complete, import the complete snapshot, verify every expected path and SHA-256 against the original manifest, fetch the remote tree back, and verify all hashes again.
- Preserve the original local commit identities where the Git bundle permits; do not describe GitHub API-created scaffold commits as equivalent to the original history.
- Keep this governance PR separate from the source-import/reproducibility PR. Do not merge PR #1 as a complete reproducibility release while its own stated synchronization acceptance criteria remain unmet.

## 5. Git and reporting

Follow the pinned central [Git merge policy](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/governance/GIT_MERGE_POLICY.md). Substantive theory and research changes default to Create a merge commit.

Follow the central [Research Task Handoff Standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/RESEARCH_TASK_HANDOFF_STANDARD.md) and [Scientific Diff reporting standard](https://github.com/lm846802297-bit/high-energy-theory-research/blob/142267bc803c90aa681110293b66b73947bbcaa1/docs/scientific-diff-workflow.md#11-research-reporting-and-communication-standard).

Report established algebraic identities, tested numerical regressions, physical validations and pending derivations separately. In particular, a projector identity regression or a cubic-vertex benchmark is not a completed five-point amplitude or a full five-point Slavnov–Taylor validation.
