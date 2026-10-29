"""P1/P2 manufactured-solution validation for the ideal quadrupole Laplace field."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from skfem import Basis, MeshTri, ElementTriP1, ElementTriP2, asm, condense, solve
from skfem.models.poisson import laplace


ROOT = Path(__file__).resolve().parents[1]


def exact(points):
    x, y = points
    return 0.5 * (x*x-y*y)


def exact_grad(points):
    x, y = points
    return np.array((x, -y))


def solve_case(n: int, order: int) -> dict:
    # The square and exact boundary data define a clean manufactured reference.
    axis = np.linspace(-6.0, 6.0, n+1)
    mesh = MeshTri.init_tensor(axis, axis)
    element = ElementTriP1() if order == 1 else ElementTriP2()
    basis = Basis(mesh, element, intorder=8)
    A = asm(laplace, basis)
    D = basis.get_dofs().all()
    g = exact(basis.doflocs)
    uh = solve(*condense(A, np.zeros(basis.N), x=g, D=D))
    X = basis.global_coordinates()
    w = basis.dx
    grad = basis.interpolate(uh).grad
    uhq = basis.interpolate(uh)
    phi = exact(X)
    phi_num=float(np.sum((uhq-phi)**2*w)); phi_den=float(np.sum(phi**2*w))
    phi_local=float(np.sum((uhq-phi)**2*w*(np.sum(X*X,axis=0)<1.0)))
    phi_local_den=float(np.sum(phi**2*w*(np.sum(X*X,axis=0)<1.0)))
    eg = exact_grad(X)
    e2 = np.sum((grad-eg)**2, axis=0)
    d2 = np.sum(eg**2, axis=0)
    err2 = float(np.sum(e2*w)); den2 = float(np.sum(d2*w))
    mask = np.sum(X*X, axis=0) < 1.0
    local_err2 = float(np.sum(e2[mask]*w[mask])); local_den2 = float(np.sum(d2[mask]*w[mask]))
    # Independent check for P1: exact element moment of the linear reference gradient.
    # With affine triangles grad(uh) is constant; quadrature above must equal this integral.
    p = mesh.p[:, mesh.t].transpose(2, 1, 0)
    centroid = p.mean(axis=1)
    a = p[:, 1]-p[:, 0]; b = p[:, 2]-p[:, 0]
    area = np.abs(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0])/2
    gu = np.zeros((len(mesh.t.T), 2))
    # Recover constant element gradient from local nodal values (valid only P1).
    p1 = Basis(mesh, ElementTriP1())
    u1 = exact(p1.doflocs)
    loc = mesh.p[:, mesh.t].transpose(2, 1, 0)
    M = np.concatenate((loc, np.ones((len(loc),3,1))),axis=2)
    coeff = np.linalg.solve(M, u1[mesh.t.T][...,None])[:, :2, 0]
    gu[:] = coeff
    variance = np.sum((p-centroid[:,None,:])**2, axis=(1,2))/12
    check_num = np.sum(area*(np.sum((gu-np.column_stack((centroid[:,0],-centroid[:,1])))**2,axis=1)+variance))
    check_den = np.sum(area*(np.sum(np.column_stack((centroid[:,0],-centroid[:,1]))**2,axis=1)+variance))
    hessian = []
    for e, dofs in enumerate(basis.element_dofs.T):
        xy = basis.doflocs[:, dofs].T
        vals = uh[dofs]
        dx = xy[:, 0]-centroid[e, 0]; dy = xy[:, 1]-centroid[e, 1]
        if order == 1:
            H = np.zeros((2, 2))
        else:
            design = np.column_stack((np.ones(len(dx)), dx, dy, 0.5*dx*dx, dx*dy, 0.5*dy*dy))
            c = np.linalg.lstsq(design, vals, rcond=None)[0]
            H = np.array([[c[3], c[4]], [c[4], c[5]]])
        hessian.append(H)
    hessian = np.asarray(hessian)
    hessian_err = np.sqrt(np.sum(area*np.sum((hessian-np.diag([1.,-1.]))**2,axis=(1,2)))/np.sum(area)) / np.sqrt(2.)
    return {
        "n_cells_per_side": n, "order": order, "dofs": int(basis.N),
        "h_max_over_r0": float(np.sqrt(2)*12/n),
        "global_potential_relative_l2": float(np.sqrt(phi_num/phi_den)),
        "local_potential_relative_l2_r_lt_1": float(np.sqrt(phi_local/phi_local_den)),
        "global_field_relative_l2": float(np.sqrt(err2/den2)),
        "local_field_relative_l2_r_lt_1": float(np.sqrt(local_err2/local_den2)) if local_den2 else None,
        "independent_p1_moment_relative_l2": float(np.sqrt(check_num/check_den)) if order == 1 else None,
        "p1_moment_vs_quadrature_abs": float(abs(np.sqrt(check_num/check_den)-np.sqrt(err2/den2))) if order == 1 else None,
        "elementwise_hessian_relative_rms": float(hessian_err),
        "boundary_geometry": "straight square boundary; no curved conductor in manufactured benchmark",
    }


def main():
    results = [solve_case(n, p) for n in (12, 24, 48) for p in (1, 2)]
    payload = {
        "scope": "M2b phase 1 only; exact harmonic quadrupole on square with exact Dirichlet boundary",
        "phi": "(x^2-y^2)/2", "r0_normalized": 1.0,
        "limitations": ["uniform affine triangles", "no isoparametric curved surfaces", "no 3D device conclusion"],
        "results": results,
    }
    (ROOT / "results" / "m2b_manufactured_validation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    fig, ax = plt.subplots(figsize=(6.8, 4.5), constrained_layout=True)
    for order, marker in ((1, "o"), (2, "s")):
        rr=[r for r in results if r["order"]==order]
        ax.loglog([r["h_max_over_r0"] for r in rr], [r["global_field_relative_l2"] for r in rr], marker+"-", label=f"P{order}: global")
        ax.loglog([r["h_max_over_r0"] for r in rr], [r["local_field_relative_l2_r_lt_1"] for r in rr], marker+"--", label=f"P{order}: $r<r_0$")
    ax.set(xlabel=r"$h_{max}/r_0$", ylabel="relative electric-field L2 error", title="Manufactured quadrupole: P1/P2 comparison")
    ax.grid(True, which="both", alpha=.3); ax.legend()
    fig.savefig(ROOT/"results"/"figures"/"m2b_p1_p2_convergence.png",dpi=180); plt.close(fig)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
