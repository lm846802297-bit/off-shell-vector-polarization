# On-shell data and the ambiguity of off-shell currents

## Statement

Let `D(q) = q^2 - m^2` and let `J^mu(q;X)` be a current whose vector leg becomes a physical external vector at `D=0`. For any local analytic vector `K^mu(q;X)`, define

`J'^mu = J^mu + D K^mu`.

At the physical pole, `J'` and `J` give identical helicity amplitudes because `D=0`. Thus on-shell S-matrix data do not determine the off-shell current uniquely.

## Effect on an internal exchange

For a channel exchange

`A = J_L^mu N_{mu nu} J_R^nu / D`,

make the local redefinitions

`J_L -> J_L + D K_L`,

`J_R -> J_R + D K_R`.

Then

`A' - A = K_L^mu N_{mu nu} J_R^nu + J_L^mu N_{mu nu} K_R^nu + D K_L^mu N_{mu nu} K_R^nu`.

The difference is regular at `D=0`; the pole residue is unchanged. In a fixed channel this is precisely the freedom to redistribute regular/contact contributions between the exchange numerator and the contact part of the amplitude.

## Consequence for a current constraint

If one completion satisfies

`q_mu J^mu - m J_chi = 0`,

then the redefined currents satisfy

`q_mu J'^mu - m J'_chi = D (q dot K - m K_chi)`

when `J_chi -> J_chi + D K_chi`.

Therefore an off-shell vector-Goldstone relation at arbitrary `q^2` is not determined by on-shell pole data alone. It requires additional structure, such as a specified BRST/gauge-fixed field basis, a prescription for contact terms, or a stronger constructibility condition that excludes the allowed local ambiguity.

## A model-level realization by field redefinition

In the SU(2)-Higgs theory with a vectorlike fermion doublet, consider the local perturbative field redefinition

`W_mu^a -> W_mu^a + alpha/Lambda^4 * (Phi^dagger Phi) * (bar(Psi) gamma_mu T^a Psi)`.

The bilinear `bar(Psi) gamma_mu T^a Psi` transforms in the adjoint, while `Phi^dagger Phi` is a gauge singlet, so the shift has the correct homogeneous transformation to redefine the gauge connection. Its canonical dimension is one. The induced change in the action is proportional to the gauge-field equations of motion, together with the corresponding gauge-fixing/BRST completion. By the equivalence theorem for local invertible field redefinitions, on-shell S-matrix elements are unchanged.

In unitary-gauge variables, `Phi^dagger Phi = (v+h)^2/2`, so this redefinition contains a term proportional to `h^2 bar(Psi) gamma_mu T^a Psi`. At the linearized level, the gauge-field equation of motion supplies the inverse transverse propagator factor `(q^2-m^2)`, realizing the `D K^mu` ambiguity in the five-point `W^* f fbar h h` current. In a general R_xi gauge the associated Goldstone and gauge-fixing terms must be transformed consistently; the transverse-sector argument alone is not the full BRST completion.

This makes the ambiguity concrete in the chosen model: it is not merely an arbitrary algebraic modification, but can arise from a change of field coordinates that leaves physical scattering unchanged. The coefficient `alpha` and the operator basis must therefore be fixed before comparing off-shell currents.

## Research implication

The target should not be phrased as finding a unique off-shell current from the S-matrix alone. A better target is:

1. characterize the equivalence class of off-shell currents modulo `(q^2-m^2)K^mu`;
2. identify the BRST-compatible representative selected by the minimal SU(2)-Higgs theory;
3. determine whether on-shell recursion plus a demonstrable large-z condition selects that representative, or only its pole-residue class;
4. keep the invariant pole residue separate from the regular contact/interference contribution.

This is a structural constraint on the research programme, not a claim that the equivalence-class observation itself is new.

## Verification

The algebraic pole-residue invariance is checked in `scripts/verify_offshell_ambiguity.py` and `results/off_shell_ambiguity.json`.
