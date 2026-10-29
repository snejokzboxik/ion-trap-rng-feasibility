import numpy as np
from skfem import Basis, ElementTetP2, MeshTet
from iontrap_rng.fem_fields import P2ScalarField


def test_p2_point_gradient_matches_exact_quadratic(tmp_path):
    axis = np.linspace(-1., 1., 5)
    mesh = MeshTet.init_tensor(axis, axis, axis)
    basis = Basis(mesh, ElementTetP2())
    H = np.array([[1.2, .2, -.1], [.2, -.7, .05], [-.1, .05, -.5]])
    points_dof = basis.doflocs
    phi = .5 * np.einsum('in,ij,jn->n', points_dof, H, points_dof)
    path = tmp_path / 'quadratic.npz'
    np.savez(path, points=mesh.p.T, tetrahedra=mesh.t.T, potential=phi)
    field = P2ScalarField(path, 'potential')
    x = np.array([[.13, .21, .3], [.1, -.2, .2], [.07, .13, -.11]])
    np.testing.assert_allclose(field.potential(x), .5*np.einsum('in,ij,jn->n', x, H, x), atol=1e-12)
    np.testing.assert_allclose(field.gradient(x), H @ x, atol=1e-11)
