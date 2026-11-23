# Equation-level review: Ke, Liu, and Mahbub (2026)

**Reference:** W. Ke, Z. Liu, and I. Mahbub, *Off-shell continuation from on-shell recursion*, arXiv:2611.13096v1, submitted 9 Nov 2026.

## Relevant equations and PDF pages

| Location | Content | Project assessment |
|---|---|---|
| Eqs. (2.1)-(2.9), PDF pp. 5-6 | Massive null decomposition, helicity projectors, transverse/longitudinal sum | Standard on-shell spin-1 basis; reference-vector dependent at the level of individual helicities |
| Eqs. (2.12)-(2.13), PDF p. 6 | Three-point off-shell current constraints | Input identities used by their recursion construction |
| Eqs. (3.10)-(3.16), PDF pp. 10-11 | Polynomial continuation of transverse and longitudinal tensors under ALT shift | Valid continuation of on-shell tensors; not physical projectors away from the pole |
| Eqs. (3.18)-(3.23), PDF pp. 12-13 | Four-fermion example | Internal longitudinal contribution rewritten using scalar-like current pieces |
| Eqs. (3.24)-(3.41), PDF pp. 13-15 | `VVff` and `VVVV` examples | Root sums cancel shift-linear terms; quadratic terms generate contact contributions |
| Eqs. (4.20)-(4.21), PDF pp. 19-20 | Four-point longitudinal split and `p dot A_4 = O(m)` condition | Four-point smooth-limit condition; direct off-shell-current derivation from their on-shell construction is left open |
| Appendix B.3, Eqs. (B.17)-(B.19), PDF pp. 27-28 | Large-z constructibility estimate | The bound assumes current constraints/Ward identities remove internal longitudinal inverse-mass enhancements; it cannot be used unconditionally to prove those same constraints |
| Appendix C.3, Eq. (C.17), PDF p. 30 | Off-shell contractions of continued tensors | Explicitly shows loss of transversality/normalization away from the massive pole |

## Direct algebraic check of Appendix C.3

Set `Delta = q^2 - m^2`, `d = q dot n`, and

`ell^mu = q^mu/m - m n^mu/d`.

Then

`q dot ell = Delta/m`,

`ell^2 = -1 + Delta/m^2`,

`q_mu Pi_T^{mu nu} = Delta n^nu/d`.

Thus the continued longitudinal vector is transverse and unit-normalized only at `q^2=m^2`. The continued transverse tensor is likewise not transverse away from the pole. The sum remains the unitary-gauge numerator `-g^{mu nu}+q^mu q^nu/m^2` for all virtualities.

These identities are checked by the repository's symbolic script.

## What the paper does not establish

1. A general n-point off-shell current constraint derived directly from on-shell constructibility.
2. A proof that all possible local boundary/contact ambiguities vanish for a generic five-point current.
3. A classification of the off-shell current ambiguity `J^mu -> J^mu + (q^2-m^2)K^mu` and a BRST prescription selecting a representative.
4. A loop-level complex-pole extension of its tree-level ALT construction.
5. A Hermitian, positive semidefinite spin density matrix for a general off-shell vector current.

The first two points define the initial research programme. The latter two are later-stage directions and must be compared against existing pole-residue and mixing literature before being claimed as novel.
