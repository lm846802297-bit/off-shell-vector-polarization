# Fixed-root Higgs shifts in the charged five-point tree

[Derivation and literature comparison](CHARGED_BOUNDARY.md) ·
[TeX](Charged_Five_Point_Boundary.tex) ·
[Figures](Fixed_Root_Boundary_Figures.pdf).

All 18 trees and the six moving residues are included. The extra cross-arm
boundary vanishes when any two external vector polarizations are orthogonal
to the shift, but the common contribution need not. A stationary scalar pole
prevents full identically vanishing boundary. These statements use physical
vectors, equal vector masses and the fixed-vector, two-Higgs null shift; they
do not assert failure of all-line transverse (ALT) recursion.

Run here using the versions in run_manifest.json:

    sha256sum -c SHA256SUMS
    OPENBLAS_NUM_THREADS=1 python charged.py > check.json
    cmp check.json results.json
    python stationary.py > stationary-check.json
    cmp stationary.json stationary-check.json
    python figures.py ../su2-fixed-root-boundary/full-results.json

The figures use the original neutral signed JSON as a separate comparison.
Nothing adjusts either amplitude to obtain agreement. Values below the plotting
floor only are clipped; the raw outputs are unmodified.
