# Verification summary — initial snapshot

Date: 2026-11-23

| Check | Result | Detail |
|---|---|---|
| Symbolic off-shell projector identities | PASS | Six exact algebraic checks in `projector_identities.json` |
| Floating-point projector regression | PASS | 10,000 deterministic random kinematic points; max residual `1.19e-12`; threshold `5e-11` |
| Local W-Goldstone-Higgs vertex contraction | PASS | Residual vanishes after momentum conservation and attached-current conservation |
| Off-shell current ambiguity | PASS | Pole residue unchanged under `J -> J + (q^2-m^2)K`; exchange difference is regular |
| Unit tests | PASS | 4 tests |

The local vertex check is not a complete five-point amplitude. The complete five-point current and its on-shell-recursive derivation remain open.
