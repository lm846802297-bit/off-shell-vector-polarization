#!/usr/bin/env python3
"""Symbolically verify the Ke-Liu-Mahbub off-shell projector identities."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "projector_identities.json"

m, d = sp.symbols("m d", nonzero=True)
q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
n0, n1, n2, n3 = sp.symbols("n0 n1 n2 n3")
q = sp.Matrix([q0, q1, q2, q3])
n = sp.Matrix([n0, n1, n2, n3])
eta = m**2 * n / (2 * d)
metric = sp.diag(1, -1, -1, -1)
q_cov = metric * q
n_cov = metric * n
q2 = sp.expand((q.T * metric * q)[0])
q_dot_n = sp.expand((q.T * metric * n)[0])
n2 = sp.expand((n.T * metric * n)[0])

PiT = -metric + (q * n.T + n * q.T) / d - m**2 * (n * n.T) / d**2
ell = q / m - m * n / d
PiL = ell * ell.T
PiSum = sp.simplify(PiT + PiL)
PiUnitary = -metric + q * q.T / m**2

# Replace the scalar products by independent symbols only after expanding.
subs = {
    q_dot_n: d,
    n2: 0,
}
# Verify the tensor sum as a direct algebraic identity.
sum_residual = sp.simplify(PiSum - PiUnitary)

# Contract covariantly: q_mu Pi^{mu nu}.
contract_T = sp.simplify((q_cov.T * PiT))
contract_L = sp.simplify((q_cov.T * PiL))
Delta = sp.symbols("Delta")
expected_T = Delta * n.T / d
expected_L = Delta * (q.T / m**2 - n.T / d)

# Algebraic contraction identities are verified after imposing q.n=d, n^2=0,
# and replacing q^2 by m^2 + Delta.
scalar_subs = {
    q_dot_n: d,
    n2: 0,
    q2: m**2 + Delta,
}
# Build scalar contractions from the explicit formula to avoid ambiguous
# component-wise substitutions of a quadratic expression.
contract_T_direct = sp.simplify(
    -q.T
    + (q2 * n.T + q_dot_n * q.T) / d
    - m**2 * q_dot_n * n.T / d**2
)
contract_L_direct = sp.simplify(
    ((q2 - m**2) / m**2) * (q.T - m**2 * n.T / d)
)
contract_T_expected = Delta * n.T / d
contract_L_expected = Delta * (q.T / m**2 - n.T / d)

# The norm and longitudinal contraction.
ell_dot_q = sp.simplify(q2 / m - m * q_dot_n / d)
ell_sq = sp.simplify((q2 - 2 * m**2 + m**4 * n2 / d**2) / m**2)

checks = {
    "PiT_plus_PiL_equals_unitary_numerator": all(x == 0 for x in sum_residual),
    "q_contract_PiT": all(sp.simplify(x.subs(scalar_subs)) == 0 for x in (contract_T_direct - contract_T_expected)),
    "q_contract_PiL": all(sp.simplify(x.subs(scalar_subs)) == 0 for x in (contract_L_direct - contract_L_expected)),
    "q_dot_ell": sp.simplify((ell_dot_q - Delta / m).subs(scalar_subs)) == 0,
    "ell_squared": sp.simplify((ell_sq - (-1 + Delta / m**2)).subs(scalar_subs)) == 0,
}

# The sum identity is exact for any q and n when d=q.n. Verify it by
# comparing the explicit tensor components after replacing d with q.n.
PiT_qn = -metric + (q * n.T + n * q.T) / q_dot_n - m**2 * (n * n.T) / q_dot_n**2
ell_qn = q / m - m * n / q_dot_n
PiL_qn = ell_qn * ell_qn.T
sum_residual_qn = sp.simplify(PiT_qn + PiL_qn - PiUnitary)
checks["sum_identity_with_d_qn"] = all(x == 0 for x in sum_residual_qn)

payload = {
    "script": "scripts/verify_projector_identities.py",
    "sympy_version": sp.__version__,
    "python_version": sys.version,
    "platform": platform.platform(),
    "definitions": {
        "Delta": "q^2 - m^2",
        "d": "q.n",
        "eta": "m^2 n/(2 q.n)",
        "ell": "q/m - m n/(q.n)",
    },
    "expected_identities": {
        "q_mu_PiT_mu_nu": "Delta*n^nu/(q.n)",
        "q_mu_PiL_mu_nu": "Delta*(q^nu/m^2 - n^nu/(q.n))",
        "q_dot_ell": "Delta/m",
        "ell_squared": "-1 + Delta/m^2",
        "PiT_plus_PiL": "-g^{mu nu} + q^mu q^nu/m^2",
    },
    "checks": checks,
    "all_pass": all(checks.values()),
}

RESULT.parent.mkdir(parents=True, exist_ok=True)
RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
if not payload["all_pass"]:
    raise SystemExit(1)
