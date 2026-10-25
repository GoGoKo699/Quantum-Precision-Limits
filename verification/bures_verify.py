"""Adversarial small-matrix checks for the Bures tomography proof.

Scope: deterministic matrix geometry and structural rank preservation.
This is NOT an implementation of the quantum Frobenius tomography primitive,
amplitude estimation, a fault-tolerant circuit, or an empirical query benchmark.

Run: python bures_verify.py --output path/to/report.json
Only NumPy is required. Without --output, the full report is printed and no
file is written. With --output, a compact summary is printed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


RNG = np.random.default_rng(20261025)
C = 0.01


def hermitian(x):
    return (x + x.conj().T) / 2


def haar_columns(rows, columns):
    x = RNG.normal(size=(rows, columns)) + 1j * RNG.normal(size=(rows, columns))
    return np.linalg.qr(x, mode="reduced")[0]


def state_factor(d, weights):
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    return haar_columns(d, len(weights)) * np.sqrt(weights)


def bures_from_factors(a, b):
    """Both factors have Frobenius norm one; avoids ill-conditioned PSD roots."""
    fidelity = np.linalg.svd(a.conj().T @ b, compute_uv=False).sum()
    return float(np.sqrt(max(0.0, 2 - 2 * fidelity)))


def truncate_rank_psd(x, r):
    values, vectors = np.linalg.eigh(hermitian(x))
    values = np.maximum(values, 0.0)
    if len(values) > r:
        values[:-r] = 0.0
    return values, vectors, (vectors * values) @ vectors.conj().T


def weighted_error(x, weights):
    diff = x - np.diag(weights)
    denominators = weights[:, None] + weights[None, :]
    return float(np.sum(2 * np.abs(diff) ** 2 / denominators))


def verify_peeling(d, weights, h, perturbation_kind):
    """Stress the geometry with a valid rank-r Frobenius oracle surrogate."""
    r = len(weights)
    a = state_factor(d, weights)
    rho = a @ a.conj().T
    residual_basis = np.eye(d, dtype=complex)
    kept_vectors, kept_values, rounds = [], [], []
    weighted_budget = 0.0
    min_progress = 1.0
    for _ in range(d + 1):
        residual = hermitian(residual_basis.conj().T @ rho @ residual_basis)
        tau = float(np.trace(residual).real)
        if residual_basis.shape[1] == 0:
            break
        # Both signs and the full allowed mass-error scale are exercised.
        mass_error = (C / 4) * np.sqrt(max(tau, 0) * h / r)
        tau_hat = np.clip(tau + RNG.choice([-1.0, 1.0]) * mass_error, 0, 1)
        if tau_hat <= 2 * h:
            assert tau <= 3 * h + 2e-13
            break
        assert tau >= h - 2e-13
        target = C * np.sqrt(tau * h / r)
        n = residual.shape[0]
        if perturbation_kind == "dense":
            z = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
            perturbation = hermitian(z)
        elif perturbation_kind == "diagonal":
            perturbation = np.diag(RNG.choice([-1.0, 1.0], size=n))
        else:
            # A rank-two coherent perturbation mixes a large-eigenvalue
            # direction and an independent direction, stressing cross blocks.
            v = haar_columns(n, min(2, n))
            perturbation = v[:, :1] @ v[:, :1].conj().T
            if n > 1:
                perturbation = (
                    v[:, :1] @ v[:, 1:2].conj().T
                    + v[:, 1:2] @ v[:, :1].conj().T
                )
        perturbation *= (0.499 * target) / np.linalg.norm(perturbation, "fro")
        values, vectors, estimate = truncate_rank_psd(residual + perturbation, r)
        error = float(np.linalg.norm(estimate - residual, "fro"))
        assert error <= target + 3e-13
        keep = values >= tau_hat / (4 * r)
        assert 1 <= keep.sum() <= r
        v_keep = residual_basis @ vectors[:, keep]
        captured = float(np.trace(v_keep.conj().T @ rho @ v_keep).real)
        progress = captured / tau
        assert progress >= 0.5 - 1e-10
        assert values[keep].min() >= tau / (5 * r) - 2e-13
        min_progress = min(min_progress, progress)
        weighted_budget += 10 * r * error * error / tau
        kept_vectors.append(v_keep)
        kept_values.extend(values[keep])
        residual_basis = residual_basis @ vectors[:, ~keep]
        rounds.append(
            {
                "mass": tau,
                "mass_estimate": float(tau_hat),
                "captured_fraction": progress,
                "retained_rank": int(keep.sum()),
                "frobenius_error": error,
                "allowed_error": float(target),
            }
        )
    else:
        raise AssertionError("The residual did not terminate within d rounds")

    v = np.concatenate(kept_vectors, axis=1)
    s = np.asarray(kept_values)
    x = hermitian(v.conj().T @ rho @ v)
    p, q = float(np.trace(x).real), float(s.sum())
    g = weighted_error(x, s)
    assert g <= weighted_budget + 3e-13
    chi = weighted_error(x / p, s / q)
    identity_rhs = q / p**2 * (g - (p - q) ** 2 / q)
    identity_residual = abs(chi - identity_rhs)
    assert identity_residual <= 3e-12
    assert (p - q) ** 2 <= q * g + 3e-13
    b = v * np.sqrt(s / q)
    distance = bures_from_factors(a, b)
    tail_bound = np.sqrt(max(0.0, 2 * (1 - np.sqrt(min(p, 1.0)))))
    triangle_bound = tail_bound + np.sqrt(max(0.0, chi))
    # Factor-based fidelity has O(machine epsilon) error in squared distance.
    assert distance * distance <= triangle_bound * triangle_bound + 2e-12

    # Verify finite factor rounding and exact normalization algebra.
    kappa = 2e-5
    step = kappa / (2 * np.sqrt(2 * b.size))
    br = (np.round(b.real / step) + 1j * np.round(b.imag / step)) * step
    factor_error = float(np.linalg.norm(br - b, "fro"))
    assert factor_error <= kappa
    br_normalized = br / np.linalg.norm(br, "fro")
    normalized_error = float(np.linalg.norm(br_normalized - b, "fro"))
    assert normalized_error <= 2 * factor_error + 2e-14
    rounded_distance = bures_from_factors(b, br_normalized)
    assert rounded_distance**2 <= normalized_error**2 + 2e-12
    assert abs(np.linalg.norm(br_normalized, "fro") - 1) <= 1e-14
    return {
        "dimension": d,
        "input_rank": r,
        "smallest_weight": float(min(weights) / sum(weights)),
        "h": h,
        "perturbation": perturbation_kind,
        "rounds": rounds,
        "minimum_captured_fraction": min_progress,
        "discarded_mass": max(0.0, 1 - p),
        "weighted_error": g,
        "weighted_upper_bound": weighted_budget,
        "normalization_identity_residual": identity_residual,
        "bures_distance": distance,
        "bures_triangle_bound": float(triangle_bound),
        "factor_rounding_error": factor_error,
        "normalized_factor_error": normalized_error,
    }


def verify_amplification_rank(d, r, env_dim, phase_count):
    """Arbitrary generalized phases must preserve the claimed rank and span."""
    weights = np.geomspace(1.0, 1e-10, r)
    a = state_factor(d, weights)
    env_isometry = haar_columns(env_dim, r)
    psi = a @ env_isometry.conj().T
    assert abs(np.linalg.norm(psi, "fro") - 1) < 1e-14
    p_basis = haar_columns(d, max(1, d // 2))
    projector = p_basis @ p_basis.conj().T
    good = projector @ psi
    bad = psi - good
    span = np.column_stack([good.reshape(-1), bad.reshape(-1)])
    z = psi.copy()
    max_span_residual = 0.0
    max_rank_tail = 0.0
    for _ in range(phase_count):
        # Dyadic rounded angles, not ideal continuously sampled angles.
        angle_p, angle_psi = np.round(RNG.uniform(-np.pi, np.pi, 2) * 2**18) / 2**18
        z += np.expm1(1j * angle_p) * (projector @ z)
        z += np.expm1(1j * angle_psi) * psi * np.vdot(psi, z)
        fit = span @ np.linalg.lstsq(span, z.reshape(-1), rcond=None)[0]
        span_residual = np.linalg.norm(z.reshape(-1) - fit)
        singular = np.linalg.svd(z, compute_uv=False)
        rank_tail = np.linalg.norm(singular[r:])
        max_span_residual = max(max_span_residual, float(span_residual))
        max_rank_tail = max(max_rank_tail, float(rank_tail))
        assert span_residual < 2e-11
        assert rank_tail < 2e-11
        assert abs(np.linalg.norm(z, "fro") - 1) < 2e-11
    return {
        "system_dimension": d,
        "input_rank": r,
        "purifier_dimension": env_dim,
        "phase_pairs": phase_count,
        "max_two_dimensional_span_residual": max_span_residual,
        "max_singular_value_tail_after_rank_r": max_rank_tail,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the full JSON report to this path.")
    args = parser.parse_args()
    profiles = [
        (8, [1]),
        (8, [1, 1]),
        (8, [.9, .099999, .000001]),
        (8, np.geomspace(1, 1e-14, 8)),
        (12, np.ones(12)),
        # Fifteen eigenvalues just below tau/(4r): retained mass is near 3/4.
        (24, [1 - 15 * .249 / 16] + [.249 / 16] * 15),
    ]
    peeling = [
        verify_peeling(d, weights, h, kind)
        for h in (1e-3, 1e-7)
        for d, weights in profiles
        for kind in ("dense", "diagonal", "coherent")
    ]
    amplification = [
        verify_amplification_rank(d, r, env, 60)
        for d, r, env in ((8, 1, 17), (8, 3, 31), (12, 6, 64), (16, 16, 127))
    ]
    result = {
        "status": "passed",
        "seed": 20261025,
        "scope": (
            "Small-matrix structural and deterministic checks only; "
            "no quantum tomography circuit or query benchmark is simulated."
        ),
        "peeling_cases": peeling,
        "amplification_rank_cases": amplification,
        "summary": {
            "peeling_case_count": len(peeling),
            "amplification_case_count": len(amplification),
            "maximum_peeling_rounds": max(len(c["rounds"]) for c in peeling),
            "minimum_captured_fraction": min(c["minimum_captured_fraction"] for c in peeling),
            "maximum_normalization_identity_residual": max(
                c["normalization_identity_residual"] for c in peeling
            ),
            "maximum_rank_tail": max(
                c["max_singular_value_tail_after_rank_r"] for c in amplification
            ),
        },
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({"status": result["status"], **result["summary"]}, indent=2))
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
