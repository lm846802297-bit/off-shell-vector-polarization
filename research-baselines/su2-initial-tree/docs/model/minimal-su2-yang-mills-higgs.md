# Minimal spontaneously broken SU(2) Yang-Mills-Higgs model

## 1. Field content and conventions

Use a four-dimensional SU(2) gauge theory with gauge field `W_mu^a`, coupling `g`, and one complex scalar doublet `Phi`. The vacuum expectation value is chosen as

`<Phi> = (0, v/sqrt(2))^T`.

The three gauge bosons acquire a common tree-level mass

`m = g v / 2`.

The physical radial Higgs mode has

`m_h^2 = 2 lambda v^2`.

For optional matter-current tests, include a vectorlike Dirac doublet `Psi` with no Yukawa coupling in the massless-fermion benchmark. This avoids introducing chiral anomaly issues and makes the Goldstone-fermion coupling vanish in that benchmark.

## 2. Lagrangian

`L = -1/4 W^a_{mu nu} W^{a mu nu} + (D_mu Phi)^dagger (D^mu Phi) - lambda (Phi^dagger Phi - v^2/2)^2 + L_f + L_gf + L_gh`.

`D_mu = partial_mu - i g T^a W_mu^a`, with `T^a = sigma^a/2`.

The field strength is

`W^a_{mu nu} = partial_mu W^a_nu - partial_nu W^a_mu + g epsilon^{abc} W^b_mu W^c_nu`.

A convenient R_xi gauge-fixing function is

`F^a = partial_mu W^{a mu} - xi m chi^a`,

with `L_gf = - F^a F^a/(2 xi)`. The Goldstone fields `chi^a` are the three would-be Goldstone modes.

## 3. Relevant tree-level couplings

With the standard normalization above:

- `W W h`: coefficient proportional to `g m g_{mu nu} delta^{ab}`.
- `W W h h`: coefficient proportional to `(g^2/2) g_{mu nu} delta^{ab}`.
- `h h h`: coefficient `-i 3 m_h^2/v` in the usual all-incoming Feynman-rule convention.
- `W chi h`: derivative coupling with magnitude `g/2`; its overall phase is fixed together with the Goldstone-field convention.
- `W f f`: `i g gamma_mu T^a` for a vectorlike doublet.

The exact signs of the W-chi-h vertex depend on the definition of `chi^a` and Fourier convention. In the derivation notes the Goldstone amplitude is phased so that the ST identity reads `q_mu J^{a mu} = m J_chi^a`.

## 4. Scope

The first benchmark is tree level and uses generic hard kinematics. The simultaneous massless scaling is

`m_i = tau * mbar_i`, `tau -> 0`,

with mass ratios held fixed. This follows the mass-singularity convention of arXiv:2611.13096v1. Soft/collinear singularities and loop infrared logarithms are outside the first benchmark.
