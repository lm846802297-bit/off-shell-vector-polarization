# Fixed-root fermion-pair shift of the SU(2)–Higgs five-point current

**Scope and pins.** This is a complementary shift, not a repetition of the Higgs-pair shift in [Draft #1 at 6d92e34b](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/6d92e34be9694920dc0c64fd0a1d18747ccb6622/research-baselines/su2-fivepoint-fixedq-boundary). The theory and graph census are unchanged: tree level, massless vectorlike f fbar, no Yukawa, two on-shell external Higgs bosons and an amputated off-shell root W*(q) or χ*(q). The common coupling/color factor and the single stationary K-line propagator are stripped, with K=p_f+p_fbar. Let V be the W-root bare Feynman amplitude contracted with q and divided by i, and G=m times the bare Goldstone-root amplitude. The ST statement for these choices is V=G. This is a Lagrangian calculation, not a claim of constructibility from physical amplitudes.

For spinor conventions and the distinction between off-shell ST and on-shell current constraints see Ke–Liu–Mahbub, [arXiv:2611.13096v1](https://arxiv.org/html/2611.13096v1), Appendix A and §4.2 (Eq. 89 and following paragraph). The shift is a two-line holomorphic spinor deformation, not their ALT shift. For complex two-line recursion in a different Yang–Mills setting see [Britto–Cachazo–Feng–Witten](https://arxiv.org/abs/hep-th/0501052); no assumption about their large-z result is applied to the present off-shell current.

## Deformation and chiral wavefunctions

Use a mostly-minus metric and represent momentum by the matrix P=p⁰ I+ p⃗·σ⃗=λ λ̃ᵀ. The vector reconstructed from a matrix is p^μ=tr(σ^μ P)/2, with σ^μ=(I,σ⃗). At complex z all spinors are holomorphic continuations; do not take z* or conjugate a shifted spinor to obtain its partner. Define

```text
lambda_f(z)       = lambda_f + z lambda_b,
tilde_lambda_b(z) = tilde_lambda_b − z tilde_lambda_f,
tilde_lambda_f, lambda_b unchanged.  Here b denotes fbar.
```

Both shifted momenta remain null by rank one, and their sum K is constant. Choose the two independent chirality choices of the external vector current, with the following fixed normalization and phases (matrix representations):

```text
J_A(z) = 2 (lambda_f + z lambda_b)(tilde_lambda_b − z tilde_lambda_f)^T
       = J_A(0) + 2z(P_b−P_f) − z² J_B,
J_B(z) = 2 lambda_b tilde_lambda_f^T = J_B(0).
```

Their four-vectors j_A and j_B obey K·j=0, or equivalently tr(adj(K) J)=0 where K inside the trace denotes its 2×2 matrix. The identity proves this without assuming a physical polarization for the virtual W line. A separate Weyl-matrix calculation uses the columns u_R=(0,λ_f(z))ᵀ and u_L=(ε λ̃_f,0)ᵀ with the rows vbar_R=(λ̃_b(z)ᵀ,0) and vbar_L=(0,−λ_bᵀ ε), taking ε_12=+1 as a matrix convention. Each column and row satisfies its shifted massless Dirac equation. The phases fix the sign needed to compare the bilinears to J_A,J_B above. Additional z-dependent little-group rescalings would change the quoted large-z grades and are not used.

## Denominators, gauge dependence and boundary

There is one conserved fermion source and all of its bosonic arguments are fixed. In particular q, H=p+r, K and Q_x=q+x for x=p,r are constant. No tree propagator pole moves in z. Each contracted amplitude is exactly a polynomial, generically degree two for A and degree zero for B. More explicitly, writing F_ξ[j]=V_ξ[j]=G_ξ[j] for the linear functional in the chosen R_ξ gauge,

```text
F_A(z) = F[J_A(0)] + 2z F[P_b−P_f] − z² F[J_B],
F_B(z) = F[J_B].
```

The equality of the ST partners follows from the complete bosonic subtree of the pinned baseline for any conserved source, not from supposing that either individual contracted amplitude vanishes at infinity. At generic kinematics both chiral currents need a boundary. For a usual contour formula for F(0) using F(z)/z, the boundary integral is precisely the constant mode F(0); the higher polynomial powers diagnose growth but do not contribute to that integral. The difference V−G vanishes identically for both wavefunctions. The lack of any moving poles under this shift means a zero moving-pole sum cannot determine the individual currents.

The polynomial can depend on the gauge. Set a=g/2, m=a v, M=m_h and d_x=Q_x². Relative to ξ=1, the one-ordering shift in the common functional is

```text
delta F_(xi,1)[j] = 2 a² m² (xi−1) q² (Q_x·j) /
                    [(d_x−m²)(d_x−xi m²)] .
```

Sum it over x=p,r; the h exchange is unchanged. Subtracting the Goldstone-root W and Goldstone exchange terms combines their (q²−M²) and M² factors into q². This is gauge variation of the **off-shell contracted currents**, not of a physical S-matrix. It is generally nonzero in the permitted hard virtual configurations, although it can vanish at special kinematics or in the algebraic q²=0 limit. It is distinct from the earlier Higgs-pair shift, whose two individual amplitudes had a common gauge-independent large-z constant despite gauge-sensitive moving-pole terms.

## Independent checks and replay

The exact SymPy check proves the matrix shift, null momenta, conservation of both currents, the quadratic coefficient and signed gauge variation. The numerical path constructs explicit shifted four-component spinors and barred rows rather than using the outer-product polynomial. It extracts the tree vertices directly from the scalar action with Fourier tags (same vertex extraction mechanism as the prior baseline), contracts complex wavefunctions, and compares the signed amplitudes at finite z and via Fourier coefficients on a circle. Nonzero quadratic coefficient and gauge variation are explicit controls. Unwanted negative-power modes and the purported flat-current growth are checked separately. Across 20 accepted fixtures and three gauge values, the explicit Dirac and two-component currents agreed to 9.94e-16 in the fixture units. The largest scaled ST residual was 7.28e-15, the gauge-variation mismatch 2.31e-15, and the forbidden contour modes 1.55e-14. The nonzero controls are order one with the archived scaling. Fixtures use timelike and crossed signed massless fermion momenta. No loops, finite widths, or general claims for other shifts or wavefunction choices enter.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_fermion_shift.py --seed 82023 --points 8 --output results-primary.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_fermion_shift.py --seed 512039 --points 12 --output results-holdout.json
sha256sum -c SHA256SUMS
```
