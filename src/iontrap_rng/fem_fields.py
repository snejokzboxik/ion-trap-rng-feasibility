"""Point evaluation of saved scalar P2 tetrahedral fields."""
from __future__ import annotations

from pathlib import Path
import numpy as np
from skfem import Basis, ElementTetP2, MeshTet


class P2ScalarField:
    """Evaluate a saved scikit-fem P2 potential and its elementwise gradient.

    Coordinates in the NPZ are dimensionless multiples of r0. The returned
    gradient is therefore d(phi)/d(x/r0), in volts per normalized coordinate.
    P2 potential is continuous; its gradient is generally discontinuous across
    element faces. The trajectory integrator must monitor face-crossing effects.
    """

    def __init__(self, npz_path: str | Path, field_name: str):
        d = np.load(npz_path)
        self.mesh = MeshTet(d["points"].T, d["tetrahedra"].T)
        self.basis = Basis(self.mesh, ElementTetP2())
        self.values = np.asarray(d[field_name], dtype=float)
        if self.values.shape != (self.basis.N,):
            raise ValueError(f"{field_name}: expected {(self.basis.N,)}, got {self.values.shape}")
        self.find_element = self.mesh.element_finder(mapping=self.basis.mapping)

    def _coordinates(self, points):
        x = np.asarray(points, dtype=float)
        if x.ndim == 1:
            x = x[:, None]
        if x.shape[0] != 3:
            raise ValueError("points must have shape (3, n)")
        cells = np.asarray(self.find_element(*x), dtype=int)
        if np.any(cells < 0):
            raise ValueError("query point lies outside the saved FEM mesh")
        X = self.basis.mapping.invF(x[:, :, None], tind=cells)
        return x, cells, X

    def potential(self, points):
        x, cells, X = self._coordinates(points)
        out = np.zeros(x.shape[1])
        edofs = self.basis.element_dofs[:, cells]
        coeff = self.values[edofs]
        for k in range(self.basis.Nbfun):
            shape = self.basis.elem.lbasis(X[:, :, 0], k)[0]
            out += coeff[k] * shape
        return out

    def gradient(self, points):
        x, cells, X = self._coordinates(points)
        out = np.zeros((3, x.shape[1]))
        edofs = self.basis.element_dofs[:, cells]
        coeff = self.values[edofs]
        for k in range(self.basis.Nbfun):
            grad = self.basis.elem.gbasis(self.basis.mapping, X, k, tind=cells)[0].grad[:, :, 0]
            out += grad * coeff[k][None, :]
        return out

