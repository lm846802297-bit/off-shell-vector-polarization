# Exact off-shell contractions of the continued spin-1 tensors

## Definitions

Use the mostly-minus metric and let `n^2 = 0`, `d = q dot n != 0`, and

`Delta = q^2 - m^2`,

`eta^mu = m^2 n^mu/(2 d)`,

`ell^mu = (q^mu - 2 eta^mu)/m = q^mu/m - m n^mu/d`.

The continued tensors are

`Pi_T^{mu nu} = -g^{mu nu} + (q^mu n^nu + q^nu n^mu)/d - m^2 n^mu n^nu/d^2`,

`Pi_L^{mu nu} = ell^mu ell^nu`.

## Derivation

Using `q dot n = d` and `n^2 = 0`,

`q_mu Pi_T^{mu nu} = -q^nu + (q^2 n^nu + d q^nu)/d - m^2 n^nu/d = Delta*n^nu/d`.

For the longitudinal tensor,

`q dot ell = q^2/m - m = Delta/m`,

so

`q_mu Pi_L^{mu nu} = (q dot ell) ell^nu = Delta/m^2 * (q^nu - m^2 n^nu/d)`.

The norm is

`ell^2 = q^2/m^2 - 2 + m^2 n^2/d^2 = -1 + Delta/m^2`.

Finally, direct expansion gives

`Pi_T^{mu nu} + Pi_L^{mu nu} = -g^{mu nu} + q^mu q^nu/m^2`.

## Interpretation

At the physical pole `Delta=0`,

- `q_mu Pi_T^{mu nu}=0`;
- `q_mu Pi_L^{mu nu}=0`;
- `ell^2=-1`.

Away from the pole, these properties fail by terms proportional to `Delta`. The continued tensors still sum to the unitary-gauge numerator, but the individual terms are not physical polarization projectors away from the pole.

This is the precise distinction required by Appendix C.3 of arXiv:2611.13096v1, PDF p. 30, Eq. (C.17)-(C.18).

## Verification

- Exact symbolic check: `results/projector_identities.json`.
- Independent floating-point scan: `results/projector_numeric_scan.json`.
- The initial mixed-index implementation failure is retained in `results/failures/` and documented in `docs/STATUS.md`.
