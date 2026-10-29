"""Targeted M2a.1 mesh/enclosure audit; no new 3-D or stochastic model."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from iontrap_rng.readout_fem import (  # noqa: E402
    _assemble_stiffness,
    _fit_center_polynomial,
    _make_rect_reference,
    _relative_l2,
    _solve_dirichlet,
    ideal_potential_v,
    solve_mesh_level,
)


def regional_field_error(points, triangles, values, grads, radius):
    p = points[triangles]
    uh_grad = np.einsum("eia,ei->ea", grads, values[triangles])
    ctr = p.mean(axis=1)
    exact_grad = np.column_stack((ctr[:, 0], -ctr[:, 1]))
    e1, e2 = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    area = np.abs(e1[:, 0]*e2[:, 1] - e1[:, 1]*e2[:, 0]) / 2
    use = np.linalg.norm(ctr, axis=1) <= radius
    centered=p-ctr[:,None,:]
    trace_cov=np.sum(centered**2,axis=(1,2))/12.0
    return float(np.sqrt(np.sum(area[use] * (np.sum((uh_grad[use]-exact_grad[use])**2, axis=1)+trace_cov[use])) /
                         np.sum(area[use] * (np.sum(exact_grad[use]**2, axis=1)+trace_cov[use]))))


def main():
    ideal_rows = []
    # Existing locally refined meshes: separate the center-sensitive volume
    # from the full square, where the far-field element size stayed fixed.
    for h in (1/20, 1/40, 1/80, 1/160, 1/320):
        p, tri, b = _make_rect_reference(h, 6.0)
        K, area, grads = _assemble_stiffness(p, tri)
        bc = {int(i): float(ideal_potential_v(p[i, 0], p[i, 1], 1.0, 1.0)) for i in b["outer"]}
        u = _solve_dirichlet(K, bc)
        ideal_rows.append({
            "mesh":"local", "hmin_over_r0":h, "nodes":len(p), "triangles":len(tri),
            "global_field_l2":_field_error(p, tri, u, grads, None),
            "field_l2_r_lt_1r0":regional_field_error(p, tri, u, grads, 1.0),
            "field_l2_r_lt_2r0":regional_field_error(p, tri, u, grads, 2.0),
            "center_hessian_relative":_hessian_error(p, u),
        })

    # A uniformly sized mesh over the same 12 r0-wide validation domain.
    # The Gmsh background field forces global refinement at each h.
    for h in (0.4, 0.2, 0.1):
        p, tri, b = make_uniform_square(h, 6.0)
        K, area, grads = _assemble_stiffness(p, tri)
        bc = {int(i): float(ideal_potential_v(p[i, 0], p[i, 1], 1.0, 1.0)) for i in b["outer"]}
        u = _solve_dirichlet(K, bc)
        ideal_rows.append({
            "mesh":"global", "h_over_r0":h, "nodes":len(p), "triangles":len(tri),
            "global_field_l2":_field_error(p, tri, u, grads, None),
            "field_l2_r_lt_1r0":regional_field_error(p, tri, u, grads, 1.0),
            "field_l2_r_lt_2r0":regional_field_error(p, tri, u, grads, 2.0),
            "center_hessian_relative":_hessian_error(p, u),
        })

    enclosure_rows = []
    for outer in (4.0, 6.0, 10.0):
        result, _ = solve_mesh_level(40 + int(outer), 1/80, outer_over_r0=outer, n_curve=100)
        enclosure_rows.append({
            "outer_over_r0":outer,
            "drive_hessian_V_per_r0_2":result.drive_quadrupole_hessian_v_per_r0_2,
            "center_gradient_fraction":result.drive_center_gradient_fraction,
            "differential_weighting_gradient_per_r0":result.differential_weighting_gradient_per_r0,
            "per_electrode_Cprime_F_per_m":result.per_electrode_capacitance_f_per_m,
            "differential_Cprime_F_per_m":result.differential_port_capacitance_f_per_m,
        })
    out = {"scope":"targeted P1 ideal-reference mesh comparison and finite-device enclosure sensitivity",
           "ideal_reference":ideal_rows, "grounded_enclosure_sensitivity_at_h_over_r0_1_80":enclosure_rows,
           "notes":["The ideal field domain has exact quadrupole Dirichlet data on the outer square.",
                    "The finite-device runs retain the grounded enclosure; no ideal quadrupole is prescribed there.",
                    "This remains P1 with straight-sided sampled boundaries; no P2 claim is made."]}
    (ROOT / "results" / "m2a1_numeric_audit.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    # Correct the saved M2a ideal-field metric with exact element integration.
    prior_path=ROOT/"results"/"m2a_readout.json"
    prior=json.loads(prior_path.read_text(encoding="utf-8"))
    local_by_h={round(row["hmin_over_r0"],12):row["global_field_l2"]
                for row in ideal_rows if row["mesh"]=="local"}
    for row in prior["mesh_levels"]:
        key=round(row["h_min_over_r0"],12)
        if key in local_by_h:
            row["ideal_relative_field_l2"]=local_by_h[key]
    prior["ideal_field_error_note"]=("Corrected M2a.1: exact element integration includes the within-triangle variance of the linear exact quadrupole gradient; prior centroid-only metric underestimated the P1 global L2 field error.")
    prior_path.write_text(json.dumps(prior,indent=2),encoding="utf-8")
    print(json.dumps(out, indent=2))


def _field_error(p, tri, u, grads, radius):
    if radius is None:
        ctr = p[tri].mean(axis=1)
        use = None
    else:
        ctr = p[tri].mean(axis=1)
        use = np.linalg.norm(ctr, axis=1) <= radius
    uh_grad = np.einsum("eia,ei->ea", grads, u[tri])
    ex = np.column_stack((ctr[:, 0], -ctr[:, 1]))
    edges = p[tri]
    e1, e2 = edges[:, 1]-edges[:, 0], edges[:, 2]-edges[:, 0]
    area = np.abs(e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0])/2
    if use is not None:
        uh_grad, ex, area = uh_grad[use], ex[use], area[use]
        edges=edges[use]
    centroid=edges.mean(axis=1)
    trace_cov=np.sum((edges-centroid[:,None,:])**2,axis=(1,2))/12.0
    return float(np.sqrt(np.sum(area*(np.sum((uh_grad-ex)**2, axis=1)+trace_cov)) /
                         np.sum(area*(np.sum(ex**2, axis=1)+trace_cov))))


def _hessian_error(points, values):
    coef = _fit_center_polynomial(points, values, rfit=0.25)
    hess = np.array([coef[3], coef[4], coef[5]])
    return float(np.linalg.norm(hess-np.array([1.0, 0.0, -1.0]))/np.sqrt(2.0))


def make_uniform_square(h, outer):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.add("uniform_ideal_audit")
    geo = gmsh.model.geo
    corners = [geo.addPoint(-outer,-outer,0,h), geo.addPoint(outer,-outer,0,h),
               geo.addPoint(outer,outer,0,h), geo.addPoint(-outer,outer,0,h)]
    curves=[geo.addLine(corners[i],corners[(i+1)%4]) for i in range(4)]
    loop=geo.addCurveLoop(curves); surface=geo.addPlaneSurface([loop]); geo.synchronize()
    outer_tag=gmsh.model.addPhysicalGroup(1,curves,tag=1)
    gmsh.model.addPhysicalGroup(2,[surface],tag=1)
    f=gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f,"F",str(h))
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.MeshSizeMin",h)
    gmsh.option.setNumber("Mesh.MeshSizeMax",h)
    gmsh.model.mesh.generate(2)
    nt,xyz,_=gmsh.model.mesh.getNodes(); points=np.asarray(xyz).reshape(-1,3)[:,:2]
    mapping={int(tag):i for i,tag in enumerate(nt)}
    types,_,node_blocks=gmsh.model.mesh.getElements(2)
    tri_blocks=[nodes.reshape(-1,3) for typ,nodes in zip(types,node_blocks) if typ==2]
    triangles=np.array([[mapping[int(t)] for t in row] for row in np.vstack(tri_blocks)],dtype=int)
    btags=gmsh.model.mesh.getNodesForPhysicalGroup(1,outer_tag)[0]
    boundary=np.array([mapping[int(t)] for t in btags],dtype=int)
    gmsh.finalize()
    used=np.unique(triangles); remap=np.full(len(points),-1,dtype=int); remap[used]=np.arange(len(used))
    return points[used],remap[triangles],{"outer":remap[boundary]}


if __name__ == "__main__":
    main()
