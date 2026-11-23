# Sewing the SU(2)–Higgs current to a second conserved fermion pair

**Scope and pins.** This follows [Draft #1 at 1e73ed55](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/1e73ed552855c064d9f9176d0ce4254f0075fd3b/research-baselines/su2-fivepoint-fermion-shift). It adds a physical external source to the earlier amputated current; it does not redo either complex shift. The field-theory comparison is Ke–Liu–Mahbub [arXiv:2611.13096v1](https://arxiv.org/html/2611.13096v1), introductory comparison of the propagators and §4.2 after Eq. 89. The calculation uses the SU(2) Lagrangian and ordinary R_ξ propagators, not a reconstruction of all multiparticle physical amplitudes from on-shell data.

Use two distinguishable vectorlike massless fermion doublets with zero scalar Yukawa couplings. Use one fixed diagonal isospin component of each doublet (so the fermion line couples to W³ through a nonzero diagonal T³), two external on-shell Higgs bosons, and standard R_ξ gauge with no quadratic W–chi mixing. Both massless external pairs are on shell. Let e carry momentum K and j carry q, all incoming, with e·K=0, j·q=0 and q+K+p+r=0. The color and fermion coupling factors of the two lines are stripped throughout. Write a=g/2, m=a v, M=m_h, and D_x=x²−m². All displayed formulas strip common phase factors. The external fermion source cannot emit a Goldstone directly in this model.

## Removing the sewing pole

At arbitrary q² away from the coincident physical/gauge pole, the root propagator tensor (phase stripped) is

```text
P(q) = [η − (1−ξ_root) q q/(q²−ξ_root m²)]/(q²−m²).
j contracted with P(q) = j/(q²−m²).
```

Thus the gauge-pole coefficient is zero *before* any special property of the vector-valued amputated current is used. The same holds for the other fermion line with K. At q²=m² the remaining pole numerator is exactly the contraction with the massive physical sum −η+q q/m² up to the conventional common propagator sign, because the q q term is killed by the conserved source. The result does not replace the off-shell amplitude by an on-shell polarization sum at arbitrary q², and the physical pole stays present. Algebraic cancellation also explains the ξ=1 degeneracy without interpreting it as two physical poles. As a nonphysical negative control, if j is deformed to violate conservation, at fixed q² the coefficient of 1/(q²−ξ_root m²) is (j·q) q/m² times the remaining subtree and the other propagator. This control is not a gauge completion of a fermion source with a Yukawa interaction.

## Internal cancellation, and why it differs from the amputated diagnostic

For one Higgs ordering x,y with Q=q+x=−K−y, set d=Q², C=d−ξ_inner m². Define E as a generic test vector on the q end and keep e·K=0. The reduced W and chi terms are

```text
T_W(E,e)   = 4 a² m² [E·e − (1−ξ_inner)(E·Q)(e·Q)/C]/(d−m²),
T_chi(E,e) = −2 a² [E·(2Q−q)](e·Q)/C.
```

Here (Q−y)·e=2Q·e follows from conservation and the derivative W–chi–h rule is fixed by the same Goldstone phase as in the prior baseline. The gauge change of their sum relative to ξ_inner=1 is

```text
δξ(T_W+T_chi) =
  2 a² m² (ξ_inner−1) (E·q)(e·Q)/[(d−m²)(d−ξ_inner m²)].
```

This is pure q at the amputated index. It can be nonzero for the previous contraction E=q. With E=j it vanishes for each ordering separately, and the result is

```text
T_W(j,e)+T_chi(j,e) = 4 a² [m²(j·e)−(j·Q)(e·Q)]/(d−m²).
```

The reduced numerator has the massive physical-projector form with the same overall sign convention. Each h exchange or seagull is already ξ independent; with H=p+r they contribute, respectively, (2 a m)(3 M²/v)(j·e)/(H²−M²) and 2 a²(j·e). **Dropping the internal chi exchange does not produce this cancellation.** It is different from not needing an external Goldstone coupling for the conserved source.

## Which graphs are in scope

The verified object is the doubly conserved exchange contribution, with one W attachment from each fermion line and the complete WW hh subtree (two W exchanges, two chi exchanges, scalar exchange, and WWhh contact). There is no other tree with the same external legs and assumptions hidden by this sewing: a connected tree with six external leaves and only cubic/quartic interactions obeys n₃+2 n₄=4. At least two cubic fermion vertices are needed. Three or more W insertions would leave at most one bosonic cubic to host both external Higgs legs, requiring the forbidden Whh or chi hh vertex, or an hhh vertex whose remaining h cannot join a W–fermion vertex without mixing/Yukawas. With distinguishable species there is no alternative pairing of the external fermion lines. This argument does not generalize to massive or Yukawa-coupled fermions, external nonsinglet scalars, or identical-fermion exchange channels. Here no general six-point gauge-invariance theorem or loop result is claimed.

## Verification and replay

The analytic result is checked symbolically. The independent numerical route extracts the scalar-action vertices as labelled Fourier coefficients, builds each uncontracted Lorentz tensor before sewing, and applies the explicit R_ξ propagators to currents constructed from Dirac matrices and on-shell external wavefunctions. It varies the gauge parameters of the two external propagators and those of the internal subtree independently as a diagnostic. The near-root test varies ξ_root through q²/m² at fixed kinematics (a test of the same denominator, not a q² deformation), while keeping the interior finite. Across 20 fixtures the internal direct-to-closed comparison was at most 3.25e-15 scaled and the independent outer-gauge comparison 2.07e-15. At offsets down to 1e-7, the unsimplified tensor evaluation amplified floating-point longitudinal source contamination to up to 8.44e-7. This raw result is retained, alongside an independent contraction of the same tensor after subtracting the measured numerical longitudinal contamination by its explicit propagator coefficient (6.83e-10 maximum scaled discrepancy). The exact Ward identity, not the numerical subtraction, is the root-pole proof. A deliberately nonconserved source produces the predicted near-root residue, and omitting the internal chi exchange produces a nonzero gauge change. The fixtures are signed real on-shell configurations, not a physical cross-section or phase-space integral.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_sewn_sources.py --seed 102023 --points 8 --output results-primary.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_sewn_sources.py --seed 560321 --points 12 --output results-holdout.json
sha256sum -c SHA256SUMS
```
