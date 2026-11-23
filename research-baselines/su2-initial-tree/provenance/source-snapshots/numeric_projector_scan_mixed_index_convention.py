#!/usr/bin/env python3
"""Independent floating-point regression for the off-shell projector identities."""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "projector_numeric_scan.json"
SEED = 20261123
SAMPLES = 1000
rng = np.random.default_rng(SEED)
metric = np.diag([1.0, -1.0, -1.0, -1.0])
max_sum = 0.0
max_t_contract = 0.0
max_l_contract = 0.0
max_norm = 0.0
min_abs_qn = float("inf")

for _ in range(SAMPLES):
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    n = np.concatenate(([1.0], direction))
    q = rng.normal(size=4)
    m = float(rng.uniform(0.2, 3.0))
    qn = float(q @ metric @ n)
    if abs(qn) < 0.05:
        q[0] += 1.0
        qn = float(q @ metric @ n)
    min_abs_qn = min(min_abs_qn, abs(qn))
    q2 = float(q @ metric @ q)
    delta = q2 - m * m
    q_cov = metric @ q
    n_cov = metric @ n

    PiT = -metric + (np.outer(q, n_cov) + np.outer(n, q_cov)) / qn - (m * m) * np.outer(n, n_cov) / (qn * qn)
    ell = q / m - m * n / qn
    PiL = np.outer(ell, metric @ ell)
    PiU = -metric + np.outer(q, q_cov) / (m * m)

    sum_res = np.max(np.abs(PiT + PiL - PiU))
    t_res = np.max(np.abs(q_cov @ PiT - delta * n_cov / qn))
    l_res = np.max(np.abs(q_cov @ PiL - delta * (q_cov / (m * m) - n_cov / qn)))
    norm_res = abs(float(q @ metric @ ell) - delta / m) + abs(float(ell @ metric @ ell) - (-1.0 + delta / (m * m)))

    max_sum = max(max_sum, float(sum_res))
    max_t_contract = max(max_t_contract, float(t_res))
    max_l_contract = max(max_l_contract, float(l_res))
    max_norm = max(max_norm, float(norm_res))

threshold = 5e-11
checks = {
    "projector_sum": max_sum < threshold,
    "transverse_contraction": max_t_contract < threshold,
    "longitudinal_contraction": max_l_contract < threshold,
    "longitudinal_vector_norm_and_contraction": max_norm < threshold,
}
payload = {
    "script": "scripts/numeric_projector_scan.py",
    "seed": SEED,
    "samples": SAMPLES,
    "threshold": threshold,
    "min_abs_q_dot_n": min_abs_qn,
    "max_abs_residuals": {
        "projector_sum": max_sum,
        "transverse_contraction": max_t_contract,
        "longitudinal_contraction": max_l_contract,
        "longitudinal_vector_norm_and_contraction": max_norm,
    },
    "checks": checks,
    "all_pass": all(checks.values()),
    "numpy_version": np.__version__,
    "python_version": sys.version,
    "platform": platform.platform(),
}
RESULT.parent.mkdir(parents=True, exist_ok=True)
RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
if not payload["all_pass"]:
    raise SystemExit(1)
