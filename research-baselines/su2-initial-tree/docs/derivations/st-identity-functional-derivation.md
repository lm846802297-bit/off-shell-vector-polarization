# BRST/ST benchmark for an off-shell vector current

## 1. BRST setup

For SU(2), with structure constants `epsilon^{abc}`, the gauge-field BRST transformation is

`s W_mu^a = (D_mu c)^a = partial_mu c^a + g epsilon^{abc} W_mu^b c^c`.

The scalar doublet transforms as

`s Phi = i g c^a T^a Phi`,

and the ghost sector has

`s c^a = -(g/2) epsilon^{abc} c^b c^c`,

`s bar(c)^a = B^a`, `s B^a = 0`.

Choose the R_xi gauge-fixing fermion so that the gauge-fixing function is

`F^a = partial_mu W^{a mu} - xi*m*chi^a`.

The gauge-fixed action is invariant under the nilpotent BRST differential. Its 1PI effective action obeys the Slavnov-Taylor functional identity.

## 2. Differentiated identity

Differentiate the functional ST identity with respect to four physical external fields `X_1,...,X_4`, then set all fields to zero. After amputation, the resulting identity contains:

- the contraction of the vector leg with its momentum;
- the corresponding Goldstone insertion;
- terms proportional to inverse two-point functions of the other external legs, or to their BRST variations.

After LSZ reduction on the four on-shell physical external states, the inverse-propagator terms vanish. In the Goldstone phase convention used in this repository, the tree-level identity is

`q_mu J_5^{a mu}(q; X_1,X_2,X_3,X_4) = m J_{5,chi}^a(q; X_1,X_2,X_3,X_4)`.

The off-shell momentum `q` is not required to satisfy `q^2=m^2`. This is the field-theory benchmark stated in the discussion following Eq. (4.21) of arXiv:2611.13096v1, PDF p. 20.

## 3. Consequence for the longitudinal mass limit

For an on-shell external massive vector, write

`epsilon_L^mu = q^mu/m - 2 eta^mu/m`, with `eta^mu = O(m^2)` at fixed `q dot n`.

Then

`A_{n,L} = (q dot J_n)/m - 2 eta dot J_n/m = J_{n,chi} - 2 eta dot J_n/m`.

If the Goldstone current and the residual `eta` contribution have a smooth simultaneous massless limit, then the scalar-like part is finite and

`q dot J_n = O(m)`.

The boundedness of the residual term must be checked for the chosen amplitude; it is not inferred from `eta=O(m^2)` alone if the current itself has inverse powers of the mass.

## 4. What is not proved by this benchmark

The ST identity is an established field-theory statement. It does not by itself provide the desired Lagrangian-independent on-shell-recursive derivation. The latter must reconstruct the same identity from factorization residues and constructibility, and must explicitly rule out a local boundary/contact ambiguity.
