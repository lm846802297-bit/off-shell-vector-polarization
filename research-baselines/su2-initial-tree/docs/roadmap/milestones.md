# Milestones and acceptance criteria

## M0 — Reproducible baseline (this snapshot)

- [x] Freeze the model and conventions.
- [x] Record the exact Ke-Liu-Mahbub equations relevant to the project.
- [x] Symbolically verify the off-shell projector identities.
- [x] Verify one local W-Goldstone-Higgs vertex contraction identity.
- [x] Verify the algebraic off-shell ambiguity and pole-residue invariance.
- [x] Run the 10,000-point numerical identity regression.
- [ ] Complete a full five-point current calculation.

## M1 — Off-shell equivalence class and five-point current

- [x] Establish the generic local ambiguity `J^mu -> J^mu + (q^2-m^2)K^mu` and its effect on the regular contact part.
- [ ] Classify allowed `K^mu` structures for the five-point SU(2)-Higgs current.
- [ ] Identify the BRST-compatible representative in a fixed R_xi field basis.



- [ ] Enumerate all tree topologies for `J_5^{a mu}(q; f fbar h h)` and `J_{5,chi}^a(q; f fbar h h)`.
- [ ] Derive all factorization residues with fixed momentum-flow and group-index conventions.
- [ ] Compute the large-z behavior under a specified shift of the four on-shell legs.
- [ ] Classify and eliminate/fix any local polynomial ambiguity.
- [ ] Compare the recursive result to an independent R_xi Feynman-rule calculation.
- [ ] Scan generic off-shell `q^2` and verify the ST identity numerically.

## M2 — General n-point induction

- [ ] State the induction hypothesis for lower-point currents.
- [ ] Prove residue cancellation in every channel.
- [ ] Establish sufficient large-z falloff and identify non-constructible sectors.
- [ ] Identify the minimum contact-term data required when constructibility fails.

## M3 — Complex-pole spin residue

- [ ] Construct the mixed vector/Goldstone inverse propagator matrix.
- [ ] Determine complex poles and left/right null eigenvectors.
- [ ] Factorize the pole residue and separate Lorentz spin from field-space mixing.
- [ ] Define a Hermitian density matrix from physical amplitudes, not from a complex projector alone.
- [ ] Check gauge-parameter independence using Nielsen identities.

## M4 — Off-shell current and collider observable

- [ ] Keep the non-resonant current tensor distinct from the pole spin density matrix.
- [ ] Retain pole-background interference explicitly.
- [ ] Compare against `H -> ZZ* -> 4l` angular moments and tomography literature.
- [ ] Quantify positivity, basis dependence, and finite-width effects.

## Promotion rules

A milestone is complete only when the derivation, tests, native outputs, environment snapshot, source/input hashes, and an independent physical validation are committed together. A failed test is preserved as a result and investigated; it is not silently overwritten.
