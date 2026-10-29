"""Numerical Floquet analysis of one ideal Mathieu coordinate."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class FloquetResult:
    monodromy: np.ndarray
    multipliers: np.ndarray
    discriminant: float
    stable: bool
    beta: float | None


def floquet_analysis(
    a: float,
    q: float,
    *,
    method: str = "DOP853",
    rtol: float = 2e-12,
    atol: float = 2e-14,
    stability_tolerance: float = 2e-9,
) -> FloquetResult:
    """Integrate the fundamental matrix across τ∈[0,π].

    Mathieu convention: u'' + (a - 2 q cos(2τ))u = 0.
    For stable motion, β=acos(Tr(M)/2)/π and ωsec=β Ωrf/2.
    """
    if not np.isfinite(a) or not np.isfinite(q):
        raise ValueError("a and q must be finite")

    def rhs(tau: float, flat: np.ndarray) -> np.ndarray:
        fundamental = flat.reshape(2, 2)
        matrix = np.array([[0.0, 1.0], [-(a - 2.0 * q * np.cos(2.0 * tau)), 0.0]])
        return (matrix @ fundamental).ravel()

    solution = solve_ivp(
        rhs,
        (0.0, np.pi),
        np.eye(2).ravel(),
        method=method,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Floquet integration failed: {solution.message}")
    monodromy = solution.y[:, -1].reshape(2, 2)
    discriminant = float(np.trace(monodromy))
    stable = abs(discriminant) < 2.0 - stability_tolerance
    beta = float(np.arccos(np.clip(discriminant / 2.0, -1.0, 1.0)) / np.pi) if stable else None
    multipliers = np.linalg.eigvals(monodromy)
    return FloquetResult(monodromy, multipliers, discriminant, stable, beta)
