# Minimal SU(2)–Higgs tree-level baseline

This public derivative preserves the research source and archived results of original local commit `a48377c055ce4929d67c4d996ab16f844458eb7a` (48 tracked files). The original Git tree was `75e0b4fbdf1554dd0828824e092cd8cacf2aaadb`. Its 48 tracked files were independently matched to the pinned source archive, and the original SHA-256 manifest verified all 47 files it covered.

This directory carries **35 unchanged files** from that snapshot and two new files: this README and `sha256-manifest.json`. It is not a byte-for-byte mirror of the full snapshot. The original bundle, its Git history, CI configuration, transient runtime metadata, and superseded repository-level records are not distributed here. The new manifest gives individual exclusions and covers every published file except itself.

## Replay

From this directory, in an environment with Python, NumPy and SymPy:

```sh
python -m unittest discover -s tests -v
```

The test suite rewrites four tracked result JSON files. Use a disposable copy when verifying a pinned source or preserving the archived outputs. An independent local replay with Python 3.12.14, NumPy 2.3.5 and SymPy 1.14.0 passed all four tests. The four generated JSON files matched the archived values, with only Python compiler/version metadata differing. Relevant input files, historical failure records and native outputs are included unchanged.

## What this establishes

- Symbolic off-shell projector identities and a deterministic algebraic numerical regression at 100, 1,000 and 10,000 samples (maximum absolute residual 1.1812773e-12 at 10,000, against the recorded 5e-11 tolerance).
- The local cubic W–Goldstone–Higgs contraction under the stated assumptions, and a separate check of the algebraic pole-residue ambiguity.
- A model and a research target, not a completed five-point physical calculation.

The full five-point Feynman amplitude, full physical Slavnov–Taylor identity, independent on-shell recursion, phase-space convergence and complex-pole extension remain unvalidated in this baseline. No GitHub Actions rerun is claimed.
