"""Ideal four-branch hyperbolic RF electrode geometry in the XY plane."""

from __future__ import annotations

from enum import Enum

import numpy as np
from numpy.typing import NDArray


class HyperbolaBranch(str, Enum):
    X_PLUS = "x_plus"
    X_MINUS = "x_minus"
    Y_PLUS = "y_plus"
    Y_MINUS = "y_minus"


def branch_points(
    r0_m: float, u: NDArray[np.float64] | float, branch: HyperbolaBranch | str
) -> NDArray[np.float64]:
    """Return branch coordinates with final dimension (x, y).

    X branches obey x²-y²=+r0²; Y branches obey x²-y²=-r0².
    """
    if not np.isfinite(r0_m) or r0_m <= 0:
        raise ValueError("r0_m must be finite and positive")
    branch = HyperbolaBranch(branch)
    u_arr = np.asarray(u, dtype=float)
    c, s = np.cosh(u_arr), np.sinh(u_arr)
    if branch is HyperbolaBranch.X_PLUS:
        x, y = r0_m * c, r0_m * s
    elif branch is HyperbolaBranch.X_MINUS:
        x, y = -r0_m * c, r0_m * s
    elif branch is HyperbolaBranch.Y_PLUS:
        x, y = r0_m * s, r0_m * c
    else:
        x, y = r0_m * s, -r0_m * c
    return np.stack(np.broadcast_arrays(x, y), axis=-1)


def branch_outward_normals(
    u: NDArray[np.float64] | float, branch: HyperbolaBranch | str
) -> NDArray[np.float64]:
    """Return unit normals from the central vacuum toward each electrode."""
    branch = HyperbolaBranch(branch)
    u_arr = np.asarray(u, dtype=float)
    c, s = np.cosh(u_arr), np.sinh(u_arr)
    norm = np.sqrt(np.cosh(2.0 * u_arr))
    signs = {
        HyperbolaBranch.X_PLUS: (c, -s),
        HyperbolaBranch.X_MINUS: (-c, -s),
        HyperbolaBranch.Y_PLUS: (-s, c),
        HyperbolaBranch.Y_MINUS: (-s, -c),
    }
    nx, ny = signs[branch]
    return np.stack(np.broadcast_arrays(nx / norm, ny / norm), axis=-1)


def hyperbola_residual(
    points_xy_m: NDArray[np.float64], branch: HyperbolaBranch | str, r0_m: float
) -> NDArray[np.float64]:
    """Return x²-y²-(+/-r0²) in square metres."""
    branch = HyperbolaBranch(branch)
    points = np.asarray(points_xy_m, dtype=float)
    if points.shape[-1] != 2:
        raise ValueError("points_xy_m must have final dimension 2")
    x, y = points[..., 0], points[..., 1]
    target_sign = 1.0 if branch in (HyperbolaBranch.X_PLUS, HyperbolaBranch.X_MINUS) else -1.0
    return x * x - y * y - target_sign * r0_m * r0_m
