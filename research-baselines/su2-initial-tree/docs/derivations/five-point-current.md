# Five-point off-shell current: baseline identity and on-shell bootstrap target

The off-shell equivalence-class issue is developed separately in [`off-shell-ambiguity.md`](off-shell-ambiguity.md). It is important because a five-point on-shell recursion can fix pole data without uniquely fixing the off-shell current representative.

## 1. Object to be studied

Let `X = {X_1, X_2, X_3, X_4}` denote four on-shell physical external states. Define the amputated tree current with one off-shell vector leg by

`J_5^{a mu}(q; X)`,

and the corresponding amputated Goldstone current by

`J_{5,chi}^a(q; X)`.

All external physical legs are on shell; `q` is not constrained to satisfy `q^2 = m^2`. The current is defined with the external vector propagator removed. The Goldstone phase is chosen to match the convention in which the three-point relation is `p dot J = m A_phi`.

## 2. Established field-theory benchmark

In the spontaneously broken SU(2) theory, the Slavnov-Taylor identity obtained by differentiating the functional ST identity with respect to the four physical external fields gives, after LSZ reduction on those physical legs,

`q_mu J_5^{a mu}(q; X) = m J_{5,chi}^a(q; X)`

at tree level, up to the phase convention just specified. Terms proportional to inverse propagators of the other external legs vanish after their on-shell LSZ reduction and contraction with physical wavefunctions.

This is an established BRST/ST benchmark, not a new result of this project. It holds at arbitrary multiplicity under the same assumptions. Its role here is to provide an independent physical target for an on-shell construction.

For a massless matter sector with no Yukawa coupling, the direct Goldstone-fermion vertex vanishes, but `J_{n,chi}` need not vanish in amplitudes containing additional physical Higgs/vector legs: the Goldstone can couple through the `W-chi-h` interaction and internal gauge/Higgs lines.

## 3. Relation to the four-point condition in Ke-Liu-Mahbub

Ke, Liu, and Mahbub write a longitudinal external-vector amplitude as

`A_{4,L} = (p_1^mu/m_1) A_{4,mu} - (2 eta_1^mu/m_1) A_{4,mu}`

and identify the four-point smooth-limit condition

`p_1^mu A_{4,mu} = O(m_1)`.

They explicitly state that the corresponding off-shell vector-Goldstone relation follows from ST identities in a Lagrangian formulation, while deriving it directly from their on-shell construction is left open (arXiv:2611.13096v1, Eq. (4.21), PDF pp. 19-20).

The five-point project target is the corresponding relation for `J_5`, derived from on-shell factorization/constructibility rather than imported as an ST identity.

## 4. Proposed recursive proof

Define the difference

`Delta_5(q; X) = q_mu J_5^{a mu}(q; X) - m J_{5,chi}^a(q; X)`.

A candidate induction argument is:

1. **Base data:** establish the three-point current relations from the minimal SU(2)-Higgs cubic amplitudes.
2. **Factorization:** enumerate every physical factorization channel of the five-point current. Express each residue as products of lower-point on-shell amplitudes and lower-point currents.
3. **Inductive cancellation:** use the lower-point current identities on the side containing the off-shell root. Check explicitly that each residue of `Delta_5` vanishes.
4. **Local ambiguity:** after all factorization residues vanish, `Delta_5` is a local polynomial/contact ambiguity. Determine its allowed Lorentz, group, mass-dimension, and high-energy structures.
5. **Constructibility:** choose a shift of the four on-shell legs that preserves their mass shells and the total momentum. Prove the required large-z falloff of `Delta_5`; this must eliminate the local boundary term rather than assume it away.
6. **Independent check:** compute the same current from tree-level Feynman rules in R_xi gauge and verify the ST relation at generic off-shell `q^2`.

The proof is not complete until steps 2-5 are explicit. In particular, vanishing residues alone do not exclude a local contact term.

## 5. Conditional on-shell recursion statement and the main obstruction

Let `z` deform the four on-shell external legs while keeping the off-shell root momentum `q` fixed. Define

`Delta_5(z) = q_mu J_5^{a mu}(z) - m J_{5,chi}^a(z)`.

If (i) `Delta_5(z)` is meromorphic, (ii) every factorization residue vanishes after the *complete* intermediate-state sum, and (iii) `Delta_5(z) = O(1/z)` at large `z`, then a contour argument gives `Delta_5(0)=0`. This is a conditional statement, not yet a proof that the assumptions hold in the SU(2)-Higgs model.

The main obstruction is that the Goldstone insertion is an auxiliary, gauge-dependent object rather than a physical asymptotic state. Its factorization cannot automatically be represented by a sum over only the three physical massive-vector polarizations. A valid recursive proof must either derive a closed factorization rule for the scalar-like/Goldstone amplitudes or show that the required unphysical-sector contributions cancel in the current difference.

There is also a potential circularity in the large-`z` estimate. Appendix B.3 of arXiv:2611.13096v1, Eqs. (B.17)-(B.19), explicitly assumes that current constraints, Ward identities, or mass factors remove inverse-mass enhancements from internal longitudinal propagators. Therefore the five-point falloff cannot be imported as an unconditional theorem while the relevant lower-point current constraints remain unproved. The induction must list exactly which lower-point constraints are used to establish the bound.

For the on-shell five-point amplitude with two external fermions and no external transverse vector, the dimensional estimate in Eq. (B.19) would give `gamma <= 4 - 5 - d_g - 2/2 = -2` for dimensionless coupling products (`d_g=0`), *provided* the internal longitudinal enhancements have already been controlled. This is a useful conditional target, not yet a proof for the off-shell current with arbitrary `q^2`.

There is a second obstruction: vanishing factorization residues leave a possible local polynomial. The large-`z` bound must be demonstrated for the *difference* `Delta_5`, not inferred from the behavior of the full physical amplitude. If a nonzero local polynomial survives, on-shell factorization alone is insufficient and additional local input is required.

This is the central test for M1: establish the complete cut-state content, prove the required lower-point bounds without circular assumptions, classify the local-current equivalence class, determine the boundary term, and only then claim a recursive result.

## 6. First concrete five-point channel

Use the minimal SU(2)-Higgs model with a vectorlike fermion doublet and consider the current with external states `f, fbar, h, h`:

`J_5^{a mu}(q; f(p_1), fbar(p_2), h(p_3), h(p_4))`.

This choice is useful because the Goldstone current has nontrivial contributions through the derivative `W-chi-h` interaction, the `W-f-f` interaction, and the Higgs self-interaction. The vector current receives the corresponding `W-W-h`, `W-W-h-h`, and exchange contributions. The complete set of tree topologies must be enumerated before claiming an explicit amplitude identity.

Momentum convention: all momenta incoming, `q + p_1 + p_2 + p_3 + p_4 = 0`. Fermion current conservation is used only where the relevant internal fermion pair is massless and on shell.

## 7. What is currently verified

- The off-shell polarization-tensor contractions are symbolically verified in `scripts/verify_projector_identities.py`.
- The local three-point W-Goldstone-Higgs kinematic identity attached to a conserved current is verified in `scripts/verify_vertex_st_identity.py`.
- Neither check is a calculation of the complete five-point current. The full five-point amplitude and its all-channel recursion proof remain open.

## 8. Falsifiable outcomes

The on-shell bootstrap programme will be considered successful only if:

- every factorization residue of `Delta_5` vanishes analytically;
- the candidate local polynomial is excluded by a demonstrated large-z bound or fixed by a separately stated physical input;
- the result agrees with the tree-level ST identity for generic off-shell `q^2`;
- the result is independent of the R_xi parameter after the physical external-state reduction;
- the smooth simultaneous massless limit has the expected scaling.
