# Boundary data for fixed-root scalar shifts in broken SU(2)

Read the [complete derivation](OffShell_Polarization_Boundary.md) for the neutral
bridge theorem and the signed result for the complete neutral Yukawa six-point
tree. [TeX source](OffShell_Polarization_Boundary.tex) is provided for review.
The charged non-Abelian extension (including a genuine cubic gauge vertex) is
not part of this verified result.

From this folder, with the surrounding repository checkout and Python, NumPy and
SymPy pinned in run_manifest.json:

    sha256sum -c PARENT.SHA256
    OPENBLAS_NUM_THREADS=1 python boundary.py > replay-conserved.json
    OPENBLAS_NUM_THREADS=1 python full.py .. > replay-yukawa.json
    cmp results.json replay-conserved.json
    cmp full-results.json replay-yukawa.json

The parent inputs already live in sibling research directories; they are not
copied here. The raw output includes signed real and imaginary components,
physical and gauge pole checks, topology-wise boundaries, and reconstruction
residuals. These are tree-level complex-shift tests, not a complex physical-pole
or loop calculation.
