"""Mesh-independent residual estimator and Dörfler marking helpers for P2 FEM.

The solver adapter must supply cellwise Laplacian residuals and face flux jumps
from the actual FEM solution. This module does not solve the PDE or certify it.
"""
from __future__ import annotations
import numpy as np


def p2_residual_indicators(cell_diameters_m, cell_volumes_m3,
                           laplacian_residual_v_m2, interior_face_areas_m2,
                           face_jumps_v_m, face_left, face_right,
                           boundary_neumann_residual_v_m, boundary_areas_m2,
                           boundary_cells):
    h = np.asarray(cell_diameters_m, float)
    vol = np.asarray(cell_volumes_m3, float)
    lap = np.asarray(laplacian_residual_v_m2, float)
    area = np.asarray(interior_face_areas_m2, float)
    jump = np.asarray(face_jumps_v_m, float)
    left, right = np.asarray(face_left, int), np.asarray(face_right, int)
    bn = np.asarray(boundary_neumann_residual_v_m, float)
    ba = np.asarray(boundary_areas_m2, float)
    bc = np.asarray(boundary_cells, int)
    n = len(h)
    if any(len(x) != n for x in (vol, lap)) or len(area) != len(jump) or len(left) != len(area) or len(right) != len(area):
        raise ValueError("inconsistent cell/face arrays")
    if np.any(h <= 0) or np.any(vol <= 0) or np.any(area < 0) or np.any(ba < 0):
        raise ValueError("lengths, volumes, and areas must be non-negative")
    eta2 = (h*h)*(lap*lap)*vol
    face_h = 0.5*(h[left]+h[right])
    contribution = 0.5*face_h*area*jump*jump
    np.add.at(eta2, left, contribution)
    np.add.at(eta2, right, contribution)
    if len(bn):
        if len(bn) != len(ba) or len(bn) != len(bc):
            raise ValueError("boundary residual and area lengths differ")
        np.add.at(eta2, bc, h[bc]*ba*bn*bn)
    return np.sqrt(np.maximum(eta2, 0.0))


def dorfler_mark(indicators, theta: float = 0.5) -> np.ndarray:
    """Return minimal descending-indicator prefix accounting for theta of eta²."""
    eta = np.asarray(indicators, float)
    if eta.ndim != 1 or np.any(eta < 0) or not 0 < theta <= 1:
        raise ValueError("indicators must be nonnegative 1D; theta in (0,1]")
    total = float(np.dot(eta, eta))
    if total == 0:
        return np.zeros(eta.shape, dtype=bool)
    order = np.argsort(eta*eta)[::-1]
    count = int(np.searchsorted(np.cumsum((eta[order])**2), theta*total, side="left"))+1
    marked = np.zeros(eta.shape, dtype=bool)
    marked[order[:count]] = True
    return marked
