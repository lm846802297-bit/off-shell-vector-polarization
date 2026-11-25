# Finite-gauge check of the charged five-point boundary

All external vectors are on shell in the single-doublet, equal-mass
broken \(SU(2)\) model. The deformation and conventions are those
of [the charged derivation](CHARGED_BOUNDARY.md). For the Cartesian doublet
\(\Phi=(v+h+i\sigma^a\chi^a)(0,1)^T/\sqrt2\) take
\(D_\mu=\partial_\mu-ia\sigma^aW_\mu^a\) and
\(\mathcal L_{\mathrm{gf}}=-(\partial\cdot W^a+\xi m\chi^a)^2/(2\xi)\).
The derivative interactions are
\[
 aW_\mu^a(\chi^a\partial^\mu h-h\partial^\mu\chi^a)
 +a\epsilon^{abc}W_\mu^a\chi^b\partial^\mu\chi^c.
\]
There is no \(WW\chi\) vertex; with no external ghosts, a ghost line
would close and is absent at tree level.

Remove \(-i\) from vector propagators. For \(E(Q)=Q^2-\xi m^2\) the
difference from unitary gauge is \(U_\xi(Q)-U(Q)=QQ/(m^2E(Q))\).
The stationary root \(K_i=k_i+p+r\) contracts with the conserved cubic
current \(J_i\), so its longitudinal difference vanishes. For a same-arm
ordered emission, the difference is
\(g_h^2(L\cdot e_i)(L\cdot J_i)/(m^2E(L)D(K_i))\).
The two \(Wh\chi\) vertices attach a Goldstone with the opposite
contribution: their factors are \(2a L\cdot e_i\) and \(-2a L\cdot J_i\),
with \(g_h/m=2a\).

For a cross-arm assignment put \(L_i+L_j+k_s=0\), \(b_i=L_i\cdot e_i\),
\(b_j=L_j\cdot e_j\). With the vertex momenta held fixed, its contractions are
\[
 V(L_i,e_j,e_s)=-D(L_j)e_j\cdot e_s-b_jL_i\cdot e_s,\quad
 V(e_i,L_j,e_s)=D(L_i)e_i\cdot e_s-b_iL_i\cdot e_s,
 \qquad V(L_i,L_j,e_s)=-m^2L_i\cdot e_s.
\]
They follow by expanding the Yang--Mills vertex and using
\(k_s^2=m^2\) and \(k_s\cdot e_s=0\). The finite-gauge correction to the
two-vector exchange, in units of \(g_h^2/m^2\), is
\[
 -\frac{b_i}{E(L_i)}e_j\cdot e_s
 +\frac{b_j}{E(L_j)}e_i\cdot e_s
 -\frac{b_i b_j}{E(L_i) E(L_j)}L_i\cdot e_s .
\]
The \(W\chi\chi\) exchange cancels the last term for each Higgs assignment.
In the first two terms a given attachment on a given arm occurs twice,
with opposite coefficients and with the same contraction of the other
two polarizations. The complete finite-gauge and unitary amplitudes
therefore agree as rational functions of the complex shift. All their
Laurent coefficients, including the boundary, agree, and the poles
at \(\xi m^2\) are removable. The argument extends to coincident pole
locations at the level of the summed rational function. It makes no
claim about each individual graph, arbitrary off-shell currents,
nondegenerate electroweak masses or loops. For the general BRST relation
see [Chen et al.](https://doi.org/10.1140/epjc/s10052-023-12093-7).

## Independent tree and contour check

Run **OPENBLAS_NUM_THREADS=1 python rxigauge_check.py** beside the pinned
charged.py. The finite-gauge tree generator enumerates labelled diagrams
from the Cartesian vertex catalogue, not from the unitary diagram list,
and evaluates the finite-gauge propagator directly without subtracting a
unitary contribution. It generates 18 no-Goldstone, six one-Goldstone,
and six two-Goldstone trees. Six real phase-space points, both null
shifts and \(\xi=0.31,1.7,4.2\) give 36 deformations; three points use an
explicit longitudinal external polarization. Across 144 finite points,
36 outer contours, 216 physical and 216 candidate unphysical residues,
the largest respective scaled discrepancies are \(2.36\,10^{-14}\),
\(3.22\,10^{-15}\), \(4.26\,10^{-10}\) and \(4.64\,10^{-10}\).
The denominator is the maximum of unity, the reference magnitude,
and the sum of magnitudes of the pieces. Individual spurious residues
reach 3.82 before cancellation. The signed amplitudes and contour
residues are retained in rxi-results.json.
