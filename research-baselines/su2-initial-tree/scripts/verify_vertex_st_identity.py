#!/usr/bin/env python3
"""Verify the local W-chi-h kinematic contraction against a conserved current.

This is a vertex-level check only. It is not the complete five-point current.
"""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "vertex_st_identity.json"

g, m = sp.symbols("g m", nonzero=True)
q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
k0, k1, k2, k3 = sp.symbols("k0 k1 k2 k3")
r0, r1, r2, r3 = sp.symbols("r0 r1 r2 r3")
j0, j1, j2, j3 = sp.symbols("j0 j1 j2 j3")
q = sp.Matrix([q0, q1, q2, q3])
k = sp.Matrix([k0, k1, k2, k3])
r = sp.Matrix([r0, r1, r2, r3])
j = sp.Matrix([j0, j1, j2, j3])
metric = sp.diag(1, -1, -1, -1)

def dot(a, b):
    return sp.expand((a.T * metric * b)[0])

momentum_conservation = [sp.Eq(q[i] + k[i] + r[i], 0) for i in range(4)]
current_conservation = sp.Eq(dot(k, j), 0)

# With the Goldstone phase fixed consistently, the derivative W-chi-h
# vertex is proportional to (g/2)*(q-r)_mu. The WWh vertex is proportional
# to g*m*g_{mu nu}. For a conserved attached current, q+k+r=0 and k.j=0
# imply (q-r).j = 2 q.j.
left = sp.Rational(1, 2) * g * m * dot(q - r, j)
right = g * m * dot(q, j)

# Substitute q = -k-r and use k.j=0.
left_reduced = sp.simplify(left.subs({q[i]: -k[i] - r[i] for i in range(4)}))
right_reduced = sp.simplify(right.subs({q[i]: -k[i] - r[i] for i in range(4)}))
residual_before_conservation = sp.factor(left_reduced - right_reduced)
current_conservation_solution = sp.solve(sp.Eq(dot(k, j), 0), j0)[0]
residual = sp.simplify(residual_before_conservation.subs(j0, current_conservation_solution))
residual_ratio = sp.simplify(residual_before_conservation / dot(k, j))

payload = {
    "script": "scripts/verify_vertex_st_identity.py",
    "python_version": sys.version,
    "platform": platform.platform(),
    "assumptions": ["q+k+r=0", "k.j=0", "Goldstone phase chosen to match q.J=m*J_chi"],
    "identity": "m*(g/2)*(q-r).j = g*m*q.j",
    "residual_before_current_conservation": str(residual_before_conservation),
    "coefficient_of_k_dot_j": str(residual_ratio),
    "residual_after_assumptions": str(residual),
    "pass": residual == 0 and residual_ratio == g * m / 2,
    "scope_warning": "This checks one local cubic-vertex contraction only; it is not a complete five-point amplitude or a proof of the on-shell-recursive current identity.",
}
RESULT.parent.mkdir(parents=True, exist_ok=True)
RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
if not payload["pass"]:
    raise SystemExit(1)
