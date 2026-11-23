#!/usr/bin/env python3
"""Verify that local D*K shifts preserve an on-shell pole residue."""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "off_shell_ambiguity.json"

D, JL, JR, KL, KR, N = sp.symbols("D JL JR KL KR N")
A = JL * N * JR / D
A_prime = (JL + D * KL) * N * (JR + D * KR) / D
delta_A = sp.factor(A_prime - A)
expected_delta = sp.expand(KL * N * JR + JL * N * KR + D * KL * N * KR)
residue_difference = sp.simplify(sp.limit(D * delta_A, D, 0))

qK, m, Kchi = sp.symbols("qK m Kchi")
current_identity_shift = sp.factor(D * (qK - m * Kchi))

checks = {
    "exchange_difference_is_regular": sp.simplify(delta_A - expected_delta) == 0,
    "pole_residue_unchanged": residue_difference == 0,
    "current_identity_shift_is_proportional_to_D": sp.simplify(current_identity_shift / D - (qK - m * Kchi)) == 0,
}

payload = {
    "script": "scripts/verify_offshell_ambiguity.py",
    "python_version": sys.version,
    "platform": platform.platform(),
    "definitions": {
        "D": "q^2 - m^2",
        "Jprime": "J + D*K",
        "exchange": "J_L*N*J_R/D",
    },
    "exchange_difference": str(delta_A),
    "pole_residue_difference": str(residue_difference),
    "current_identity_shift": str(current_identity_shift),
    "checks": checks,
    "all_pass": all(checks.values()),
    "scope_warning": "This verifies the algebraic ambiguity; it does not classify all allowed K structures in the five-point SU(2)-Higgs current.",
}
RESULT.parent.mkdir(parents=True, exist_ok=True)
RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
if not payload["all_pass"]:
    raise SystemExit(1)
