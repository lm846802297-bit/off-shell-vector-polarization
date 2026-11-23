# Nonconserved on-shell Yukawa source sewn to the SU(2)–Higgs current

**Scope and pins.** This is a new external source for the [conserved-source checkpoint at ecaa3fc6](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/ecaa3fc63158d3f2f0564b171536984c30138ef6/research-baselines/su2-fivepoint-sewn-sources). The same WW hh / chi–W–hh subtree is used with the second fermion pair still massless and conserved. The tree-level exchange subset reported here is not the complete six-point amplitude. For the field-theory distinction between the off-shell ST relation and on-shell recursion see Ke–Liu–Mahbub [arXiv:2611.13096v1](https://arxiv.org/html/2611.13096v1), §4.2 after Eq. 89.

## A physical, nonconserved source with a fixed phase

Introduce a distinguishable massive chiral fermion f whose left component is the lower member of an SU(2) doublet (T³=−1/2) and whose right component is a singlet. The mass comes from the real Yukawa coupling −y Qbar_L Φ f_R + h.c. With a=g/2, m=a v, M=m_h and μ=y v/√2, take Φ_lower=(v+h+iχ)/√2 and the same D=∂+iaW³ convention as the earlier baseline. Spectator doublets can be chosen to cancel global SU(2) anomalies; they do not appear in this tree example. Work with a fixed external diagonal isospin component and no quadratic W–χ mixing in R_ξ gauge.

The neutral mass/Yukawa terms are

```text
−μ fbar f − (μ/v) h fbar f − i(μ/v) χ fbar γ5 f.
```

Consequently the bare vertices are Γ_W^α=−i a γ^α P_L and Γ_χ=(μ/v) γ5. The phases are from the Lagrangian, not chosen after sewing. Write j^α=vbar(p_b)γ^α P_L u(p_f). Both massive momenta are incoming to the source vertex. Their sum q=p_f+p_b flows out of the source into the right subtree, so the propagator momentum entering the source vertex in the all-incoming rule is −q. The momentum and spinor signs are specified concretely by p_f,p_b future-directed and slash(p_f)u=μu, vbar slash(p_b)=−μvbar. Thus

```text
u_div = q·j = μ vbar (P_R−P_L) u = μ vbar γ5 u,
vbar Γ_χ u = a u_div/m.
```

The source Ward relation is not ordinary conservation; u_div is generically nonzero for the massive spinors. Let K be the momentum sum on the other fermion line, e its conserved current, p and r the Higgs momenta, all incoming, q+K+p+r=0. The common far-line W propagator and fermion factors are stripped.

## Root cancellation and inner gauge sensitivity

Write the bare vector-root subtree as iT and the bare scalar-root subtree as G. With V=q·T its known ST relation is V=mG. Use D=q²−m² and C_r=q²−ξ_r m² for the root propagator. Multiplying the bare vertices and propagators, and dividing the full expression by the common factor −i a at the near source, gives

```text
A_root W = [ j·T − (1−ξ_r) u_div V/C_r ] / D,
A_root χ = −u_div G/(m C_r).
```

Thus the sum simplifies exactly to

```text
A = [j·T − u_div V/m²] / (q²−m²).
```

The two gauge-pole residues are ±u_div V/m² at C_r=0, and cancel at fixed finite interior parameters. Keeping only W would give a generally nonzero residue. The physical W propagator denominator remains. This check also applies when the gauge is common to all lines, provided one stays at generic values of the other invariants. The independent gauge choices in verification are a diagnostic that distinguishes the inner and outer factors, not a claim that the model chooses separate physical gauges.

Now use the previously verified tensor variation of the right subtree with e·K=0. Defining Q_x=q+x for x=p,r, its change relative to ξ_i=1 is pure q on the amputated near index:

```text
δ T(E) = c (E·q),
c = 2 a² m² (ξ_i−1) sum_x (e·Q_x)/
                             [(Q_x²−m²)(Q_x²−ξ_i m²)].
```

For the conserved source of the earlier checkpoint this vanished. Here the numerator acting on T is effectively E=j−u_div q/m², with E·q=−u_div(q²−m²)/m². Therefore the exchange subset keeps generally nonzero **inner** variation δA=−u_div c/m² even after the root pole has cancelled. This is not a gauge-dependent physical S-matrix result. The graph census of the no-Yukawa case no longer applies: a Higgs can be emitted through a Yukawa insertion anywhere along the massive fermion line, with the other Higgs emitted from the W/χ bridge or the same fermion line. A virtual Higgs emitted by the massive line can also split via hhh. Such contributions and their Goldstone-related partners are absent here. Their calculation and the full six-point cancellation are not inferred just from the subset residual.

## Independent computation and replay

The code constructs actual on-shell massive spinor columns and rows via the Dirac operators and derives the scalar-source vertex by differentiation of the signed Yukawa mass matrix. The bosonic subtree is computed directly by differentiating the scalar action with labelled Fourier waves, retaining separate W and χ root terms. It compares action-derived values with the exact root formula at separately varied or common gauges. Near C_r=0 it evaluates the unsimplified W and χ terms separately on both sides at fixed physical momenta, and records their residue estimates with signs. The W-only contribution is the explicit negative control. The finite nonzero inner gauge variation is independently checked against the tensor formula. Over 20 accepted fixtures, the source Ward residual was 7.78e-16, the same-gauge root-pair difference 3.15e-16 scaled, and the summed pole-residue residual 5.67e-15 scaled. The direct finite evaluation down to |C_r|=1e-6 agreed with the simplified root formula to 8.38e-9 scaled, as expected from cancellation of large terms. The inner variation formula agreed to 1.67e-16 scaled while the largest relative nonzero variation was 0.201 in the archived scale. No phase is fitted, and there is no massless-current result rerun apart from the far-current conservation check.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_yukawa_source.py --seed 112023 --points 8 --output results-primary.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_yukawa_source.py --seed 163401 --points 12 --output results-holdout.json
sha256sum -c SHA256SUMS
```
