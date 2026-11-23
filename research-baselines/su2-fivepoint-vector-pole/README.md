# A physical-vector factorization residue of the five-point SU(2)–Higgs current

**Scope:** the `ffbar | W*(q) h h` physical-vector pole in the five-point current with external physical states `f, fbar, h, h`. This checks one residue of the current identity, not the full five-point amplitude or a recursion proof. Baseline target: [five-point-current.md](https://github.com/lm846802297-bit/off-shell-vector-polarization/blob/2f51ced35ea534983ca15b63639e5bfcb9fc5853/research-baselines/su2-initial-tree/docs/derivations/five-point-current.md). The reference is Ke–Liu–Mahbub, [arXiv:2611.13096v1, §4.2, current constraint and following paragraph](https://arxiv.org/html/2611.13096v1#S4.SS2) (Eq. 89 in HTML, labeled Eq. (4.21) using the baseline's PDF notation). The paper discusses the on-shell smooth-limit constraint and identifies the Lagrangian ST benchmark; this calculation does not import that identity as a step.

## Physical pole and conventions

Use `g_mu_nu = diag(1,-1,-1,-1)`, `a=g/2`, `m=a v` and `M=m_h`. The fermion is a massless vectorlike doublet without a Yukawa coupling. In all-incoming notation put `K=p_f+p_fbar`; the right side has `q+K+p+r=0` and `p²=r²=M²`; `q²` is arbitrary. The physical factorization pole is at `K²=m²`. The left subtree is the conserved current `L^mu = g ubar gamma^mu T^b v` (up to an overall crossing phase). The massless Dirac equations imply `K·L=0`, and the auxiliary Goldstone cannot couple directly to the left pair. Consequently the longitudinal gauge-propagator part, including its possible pole coincident at `xi=1`, does not contribute on the left. Ghost number forbids a single ghost on this cut. The physical residue can be computed using the three physical cut polarizations `e_lambda`, each satisfying `K·e_lambda=0`. The overall common pole convention and the left colour matrix are left factored out. The right amplitude is diagonal in the cut/root adjoint index; with two external Higgs singlets the nonzero diagrams are captured by the fixed-colour Abelian subsector.

Take the lower scalar-doublet component `(v+h+i chi)/sqrt(2)` for the third gauge axis, with `D=partial+i a W` on that component, `F=partial·W-xi m chi`, and Fourier convention `exp(-i p·x)`. The gauge fixing cancels quadratic (W)–Goldstone mixing. Denote the amputated bare-phase Feynman amplitudes by `i V^mu` and `G`, respectively. With this convention the repository's real-phase ST relation is `q·V=m G` (equivalently `q·(i V)=i m G`). The cubic rules are

| Rule | Value |
| --- | --- |
| `W W h` | `2 i a m g_mu_nu` |
| `W chi h` | `a (p_chi - p_h)_mu` |
| `h h h` | `-i 3 M²/v` |
| `h chi chi` | `-i M²/v` |
| `W W h h` (quartic) | `2 i a² g_mu_nu` |

No `W chi h h` quartic occurs. Four external legs imply either two cubic vertices or one quartic. For the chosen particles that gives two permutations each of W and chi exchange, one h exchange, and, only on the vector-root side, the seagull. A non-Abelian W³ vertex cannot enter these diagrams: its opposite end in a four-point tree would require a `W h h` vertex, which is absent. Tree diagrams with a closed ghost line are not part of this four-point subtree.

## Residue by diagram algebra

For one Higgs ordering, `Q=q+p=-(K+r)` and let `D=Q²-m²` and `C=Q²-xi m²`. The R_xi W propagator is `-i [g_mu_nu-(1-xi) Q_mu Q_nu/C]/D` and the chi propagator is `i/C`. With `c=(1-xi)/C`, the two W-exchange contributions on the right are

```text
q·V_W = 4 a² m² [ q·e - c (q·Q)(Q·e) ] / D
m G_W = 2 a² m² [ (q-p)·e - c ((q-p)·Q)(Q·e) ] / D .
```

Their difference simplifies to `2 a² xi m² (Q·e)/C`. The chi-exchange terms are

```text
q·V_chi = -a² [q·(Q+p)] (Q-r)·e / C
m G_chi =  a² M² (Q-r)·e / C .
```

Using `q·(Q+p)=Q²-M²` and `(Q-r)·e=2 Q·e`, their difference is `-2 a² Q² (Q·e)/C`. The combined difference for this one ordering is therefore the nonzero *regular* term `-2 a² Q·e`; the same relation holds with the Higgs legs exchanged. The sum is `-2 a² q·e` because `(2q+p+r)·e=(q-K)·e=q·e`. It cancels the vector seagull `+2 a² q·e` exactly. The remaining h exchange cancels separately: with `H=p+r` and `B=3 M²/[v(H²-M²)]`, its difference is `a m B [2q·e-(q-H)·e]=0`.

This gives `q·V-m G=0` for **each** physical cut polarization. Thus

```text
Res_(K²=m²) [q·J_5 - m J_5,chi]
  = sum_lambda L_lambda [q·V(e_lambda) - m G(e_lambda)] = 0
```

in the stated phase convention and for general `xi` as a rational identity away from coincident right-subchannel poles. The scalar-like/Goldstone amplitude is not dropped in the right subtree. In particular, a residue-based argument that kept the right exchange residues but omitted the seagull would leave a regular term.

## Independent computation and replay

The independent route extracts the cubic and quartic vertices directly as coefficients of labelled Fourier waves in the scalar action, using subset polarization rather than the reduced Feynman rules above. It then contracts explicit 4-component tensors with the R_xi propagators. A separate Dirac gamma-matrix trace constructs the left spin-summed tensor, and a Lorentz-boosted three-polarization basis reconstructs the cut. The program compares the two methods **per graph before summation**. Both controls remove a required contribution and produce nonzero residues: omitting the vector seagull and omitting the chi exchange from the Goldstone-root current.

The archived runs include 24 and 72 accepted generic kinematic points, respectively, at `xi=0.37, 1, 2.1, 5.2` and with three polarizations each. These check the stated model and pole, not physical phase-space integration. Across the two runs the worst scaled per-graph discrepancy was `2.45e-15`; the worst contracted ST difference was `3.11e-15`; the scaled spin-summed squared difference was below `4.71e-29`. The two omission controls produced scaled residuals of order unity. The JSON includes the actual masses, four-vectors, and per-point spin-summed quantities.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_vector_pole.py --seed 1209 --points 24 --output vector-pole-results.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_vector_pole.py --seed 36290 --points 72 --output vector-pole-results-holdout.json
sha256sum -c SHA256SUMS
```

**Still open:** the remaining five-point channels, the regular part of the full five-point current, and the fixed-root large-shift argument that could eliminate or determine a five-point boundary term. The seagull fixed here is a **four-point subtree contact term**; cancellation of this particular pole does not fix a local five-point ambiguity or constitute a Lagrangian-independent ST derivation.
