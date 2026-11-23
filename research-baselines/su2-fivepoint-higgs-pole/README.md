# The physical Higgs pole of the five-point SU(2)–Higgs current

**Scope:** the `hh | W*(q) f fbar` channel in the tree-level current with four physical external legs `f, fbar, h, h`. It is a different cut from the [previous physical-vector checkpoint](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/69e6f49e089c0e05eefb51a6b6ad071354843390/research-baselines/su2-fivepoint-vector-pole). The model and phase convention are those of the [pinned baseline](https://github.com/lm846802297-bit/off-shell-vector-polarization/blob/69e6f49e089c0e05eefb51a6b6ad071354843390/research-baselines/su2-initial-tree/docs/derivations/five-point-current.md). Ke–Liu–Mahbub, [arXiv:2611.13096v1 §4.2](https://arxiv.org/html/2611.13096v1#S4.SS2), Eq. 89 and the text following it, distinguish their on-shell current constraint from the off-shell Lagrangian Slavnov–Taylor relation. We use Lagrangian vertices, not that relation, as input.

## Amplitudes and residue

Take all momenta incoming, `H=p+r`, `K=p_f+p_fbar`, `q=-H-K`. Use metric `(+---)`, `a=g/2`, `m=a v`, `M=m_h`. The fermion doublet is massless with no scalar Yukawa coupling. On the two-Higgs side of the cut, `Whh` and `chi hh` are absent and only `hhh` contributes. It therefore isolates the physical scalar pole `H²=M²` at generic other invariants. The other subtree is `W*(q) h f fbar` (or `chi*(q) h f fbar)); only a virtual W coupling to the fermions can connect it. No other tree diagram generates this scalar pole. This proof is unaffected by possible accidental equality with an unphysical scalar mass, because its two-Higgs cubic vertex vanishes.

Write `B=3 M²/v`; the cubic h vertex is `-i B` and the h propagator is `i/(H²-M²)`. The W and chi roots see cubic vertices `2 i a m eta_mu,nu` and `a (q-H)_nu` under the baseline convention `D=partial+i a W` and Fourier convention `exp(-ip·x)`. Strip only the common `g T^a_colour` factor from the fermion vertex; its remaining bare factor is `i gamma^nu`. Use

```text
D_K = K²-m²,   C_K = K²-xi m²,
P_mu,nu = eta_mu,nu - (1-xi) K_mu K_nu / C_K.
```

The two bare-phase residue operators on the open fermion line are, with the W root already contracted with q,

```text
R_W   = 2 i a m B q^mu P_mu,nu gamma^nu / D_K,
R_chi =     a B (q-H)^mu P_mu,nu gamma^nu / D_K .
```

With the repository's real-phase ST convention the bare difference is `R_W-i m R_chi`. It is **not zero before projection**:

```text
R_W-i m R_chi
    = i a m B (q+H)^mu P_mu,nu gamma^nu / D_K
    = -i a m B xi slash(K) / C_K .
```

Consequently the residue vanishes between the massless on-shell Dirac wavefunctions, for any `xi` as a rational identity away from coincident denominators in the other invariant. In matrix form, `slash(p_fbar) slash(K) slash(p_f) = 0` because `K=p_f+p_fbar` and both endpoint slashes square to zero. Each *projected* amplitude is separately independent of `xi` here because every gauge-dependent part of the virtual W contracts the conserved external fermion current. Both projected residues are individually nonzero at generic kinematics. The cancellation is not a claim about unprojected off-shell Green functions.

## Independent contour route

At the pole the equal-mass on-shell kinematics of `h -> hh` are not real. Use a complex analytic continuation, without conjugating kinematic components. Let `z=H²-M²`, fix the fermion momenta, and use a fixed Lorentz boost for the scalar pair:

```text
H_rest(z) = (sqrt(M²+z), 0),
p_rest(z) = (sqrt(M²+z)/2,   i n sqrt(M²-(M²+z)/4)),
r_rest(z) = (sqrt(M²+z)/2,  -i n sqrt(M²-(M²+z)/4)),
q(z) = -K-H(z),   n·n (Euclidean) = 1.
```

The numerical method extracts cubic coefficients of the scalar action directly using a signed Fourier-tag sum, multiplies the bare vertices and propagators at finite `z`, and uses explicit Dirac matrices. It extracts the pole by averaging `z A(z)` around a contour, comparing **complex open-line matrices** before any projection to the independent factorized formulas. Two contours use radius `0.12 M²`/12 nodes and `0.055 M²`/20 nodes. All contour momenta obey the external mass shells, and the square-root branch points lie outside both radii. The crossed half of the fixtures uses a past-directed signed null momentum in the algebraic projection, not a positive-norm spin density. This contour deforms `q`; it is **not** the fixed-root large-shift contour for the complete amplitude. Only the diagrams capable of having this pole were included.

The two native JSON runs record projected and unprojected evidence for 16 and 32 accepted fixtures at `xi=0.41, 1, 2.35, 4.8`, including timelike and spacelike external-root virtualities. Across the two runs, the worst scaled contour-to-factorization disagreement was `9.59e-14`, the projected ST residual `1.25e-14`, and the projected gauge variation `3.52e-14`. Wrong Goldstone-phase and skipped-Dirac-projection controls fail at order unity. The contour and threshold are numerical verification, not an estimate of the theory’s accuracy. The first fixture of each run includes the signed matrices for every gauge value. Exact SHA-256 fingerprints and native numerical summaries are included.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_higgs_pole.py --seed 20261123 --points 16 --output higgs-pole-results.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_higgs_pole.py --seed 41037 --points 32 --output higgs-pole-results-holdout.json
sha256sum -c SHA256SUMS
```

**Still open:** the other five-point channels, the regular/fixed-root boundary contribution, the full on-shell-recursive construction, and loop or complex physical pole results. Complex kinematics in this test are a tree-level analytic tool, not a finite-width pole prescription.
