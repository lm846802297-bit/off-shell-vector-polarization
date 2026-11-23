# Fixed-root boundary for the SU(2)–Higgs five-point current

**Pinned inputs:** [Draft #1 at 46cb5f6](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/46cb5f6781b4c1fedf85ff2407d33bad2f41ff5f), the [physical W pole](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/46cb5f6781b4c1fedf85ff2407d33bad2f41ff5f/research-baselines/su2-fivepoint-vector-pole) and [hh Higgs pole](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/46cb5f6781b4c1fedf85ff2407d33bad2f41ff5f/research-baselines/su2-fivepoint-higgs-pole). Compare Ke–Liu–Mahbub [arXiv:2611.13096v1 §4.2](https://arxiv.org/html/2611.13096v1#S4.SS2), Eq. 89 in HTML, and [Appendix B](https://arxiv.org/html/2611.13096v1#A2). This uses their distinction between the Lagrangian current relation and a derivation from on-shell recursion. The shift and the Lagrangian input here are not the paper’s on-shell ALT construction.

## Exact scope and graph census

Take the same minimal SU(2) model, massless vectorlike fermion doublet, no scalar Yukawa, no quadratic W–Goldstone mixing in the R_xi propagators, and all-incoming momenta. Let `K=p_f+p_fbar`, `H=p+r`, `q=−K−H`, `a=g/2`, `m=a v`, `M=m_h`, and `p²=r²=M²`. The current root is off shell and amputated. The external Dirac spinors are on shell, so the current `e^mu` carrying the off-shell exchange momentum K obeys `K·e=0`; e is **not** an on-shell polarization of K. Strip only the common virtual K-line factor, fermion coupling and color matrix; none depends on the shift.

The open fermion line has only W insertions. Two or more insertions on a connected tree would partition the three external bosonic legs into independent nonempty branches. With two insertions, either one branch is a lone h, which cannot attach by W–h mixing, or a lone root is opposite an hh branch, which would require a `Whh` vertex. With three insertions there are again lone h branches. A connection between such branches would be a loop on the graph. Hence only the one-W attachment survives. The remaining W–W (or chi–W)–hh subtree contains two orderings of W exchange, two of chi exchange, a scalar h exchange, and the WWhh contact on the W-root side. The would-be W³-based cubic pairing would require the absent Whh vertex at its other end. This is a full-tree diagram census **for this specified channel/model only**, not a statement about other currents or external states.

The reduced contracted W amplitude `V(z)` is the bare contracted W Feynman amplitude divided by i; write `G(z)` for the bare Goldstone-root amplitude times m. The comparison `V−G` is the Slavnov–Taylor difference in this convention. It is not the uncontracted current.

## Complex shift and boundary

Use `p(z)=p+z n`, `r(z)=r−z n` with `n²=n·p=n·r=0` and `q·n!=0`. Such a complex n can be chosen in the hh rest frame using two orthogonal unit vectors perpendicular to the Higgs three-momenta; take `n=(0,u+i v_perp)` then boost back. Thus q, K, H remain fixed and the Higgs legs stay on shell. Work away from the stationary singularities and generic moving-pole collisions.

In one exchange ordering let `x(z)=x+s z n`, `Q=q+x(z)`, `d=Q²=d0+s L z`, `t=Q·e=t0+s U z`, `A=q·e`, `L=2 q·n`, `U=n·e`, with `s=+1` or `−1`. The reduced Goldstone terms obtained from the cubic vertices are

```text
G_W   = 2a²m² [ 2A − t − (1−xi)(q²−M²)t/(d−xi m²) ]/(d−m²)
G_chi = 2a²M² t/(d−xi m²).
```

Per ordering the exact large-z constant is `2a²(M²−m²)U/L`; there is no positive power left in G. The corresponding W-root terms satisfy `V_W+V_chi−G_W−G_chi=−2a²t` before summing the orderings. Indeed, the two terms of this difference are `2a² xi m² (Q·e)/(d−xi m²)` and `−2a² d (Q·e)/(d−xi m²)` after external-mass and current conservation identities. The two ordering remainders add to `−2a²A` even at finite z, and the WWhh vertex supplies `+2a²A`. The stationary Higgs graph equals `2 a m B A` on both sides, where `B=3M²/[v(H²−M²)]`. Therefore `V(z)=G(z)` for this complete reduced tree, with equal generally **nonzero** boundary:

```text
V_infinity = G_infinity
           = 2 a m B (q·e) + 2 a² (M²−m²) (n·e)/(q·n).
```

The stationary h term in this boundary still carries its (H^2=M^2) pole when viewed as a function of all invariants; it is not a local five-point polynomial. The boundary is shift-dependent in general, although independent of the R_xi parameter. The two separate currents at finite z do **not** have to be gauge independent: at an isolated gauge pole `d(z_g)=xi m²` with `xi!=1`, the residue is `2a² q² t(z_g)/(s L)` in both V and G. It is generically nonzero for an off-shell q, and cancels in their difference. Do not treat every moving pole of these separate amplitudes as a physical-state pole.

## Independent finite-z and Cauchy checks

The numerical path extracts cubic/quartic coefficients directly from labelled Fourier waves in the scalar action, with explicit complex momenta and propagators; it does not call the closed boundary or residue formulas. The Fourier action coefficient method is the same independently implemented method used for the earlier W-pole checkpoint, now evaluated at fixed-root complex-shifted momenta and **generic off-shell K**. For each fixture e is formed using explicit gamma matrices and on-shell endpoint slashes. The computation samples generic timelike and spacelike off-shell configurations, two outer circles to recover the z-regular boundary, and smaller disjoint circles around each expected moving pole to determine its residue. It independently compares the residue sum alone to the full value at z=0, then adds the boundary. The numeric gauge poles at xi=1 are merged rather than counted twice. Over 20 accepted fixtures at three gauge values, the largest scaled finite-shift ST difference was 7.05e-13, boundary-to-formula mismatch 3.41e-13, and full reconstruction mismatch 6.99e-12. Omitting the boundary gave a scaled mismatch of order one; these are numerical checks, not a theory uncertainty. Keeping complex signed data avoids interpreting a cancellation within norms as an amplitude identity.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_fixedq_boundary.py --seed 62023 --points 8 --output fixedq-results.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_fixedq_boundary.py --seed 287431 --points 12 --output fixedq-results-holdout.json
sha256sum -c SHA256SUMS
```

**Limits:** This is a Lagrangian verification in one external-state sector and at one class of shifts. It establishes neither a Lagrangian-independent on-shell recursion nor absence of possible boundary terms for other shifts, current components or theories. There is no loop or complex physical-pole result.
