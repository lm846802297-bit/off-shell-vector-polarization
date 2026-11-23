# Preserved initial failure

`projector_numeric_scan_mixed_index_convention_failure.json` records the first floating-point implementation attempt. The code represented the projector with a mixed upper/lower index while the contraction target was interpreted as a fully contravariant tensor. The projector-sum identity passed, but the transverse contraction check failed.

The exact source snapshot is preserved in `provenance/source-snapshots/numeric_projector_scan_mixed_index_convention.py`. The implementation was corrected to use fully contravariant `Pi^{mu nu}` components. The corrected symbolic and numerical checks pass. This is a software/index-convention failure, not a counterexample to the paper's analytic identity.
