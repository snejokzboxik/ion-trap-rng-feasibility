"""Bounded 2-D P1 FEM electrostatics for finite hyperbolic trap electrodes.

All geometry is dimensionless (r0 = 1).  Gmsh supplies a conforming, locally
refined triangular mesh; SciPy assembles and solves Laplace's equation.  The
outer grounded shield belongs only to the finite-device model.  The ideal
solver validation instead applies the exact quadrupole Dirichlet data there.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asinh, cosh, pi, sinh, sqrt

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from .constants import EPSILON_0_F_PER_M, K_B_J_PER_K
from .ions import IonSpecies


def _gmsh_module():
    """Load Gmsh only when a mesh is requested (keeps M1 algebra lightweight)."""
    try:
        import gmsh
    except ImportError as exc:
        raise RuntimeError("install the optional FEM dependency: pip install .[fem]") from exc
    return gmsh


@dataclass(frozen=True)
class MeshResult:
    level: int
    h_min_over_r0: float
    nodes: int
    triangles: int
    ideal_relative_potential_l2: float
    ideal_relative_field_l2: float
    ideal_relative_hessian: float
    drive_quadrupole_hessian_v_per_r0_2: float
    drive_center_gradient_fraction: float
    drive_quadrupole_residual_fraction: float
    drive_higher_multipole_fraction_at_0_4r0: float
    single_electrode_weighting_gradient_per_r0: float
    paired_rf_weighting_gradient_per_r0: float
    differential_weighting_gradient_per_r0: float
    per_electrode_capacitance_f_per_m: float
    differential_port_capacitance_f_per_m: float


def hyperbola_branch(parameter: np.ndarray | float, branch: str, r0: float = 1.0) -> np.ndarray:
    """Return one of the four required hyperbola branches, parameterized by t."""
    t = np.asarray(parameter, dtype=float)
    if branch == "right":
        return np.vstack((r0 * np.cosh(t), r0 * np.sinh(t)))
    if branch == "left":
        return np.vstack((-r0 * np.cosh(t), r0 * np.sinh(t)))
    if branch == "top":
        return np.vstack((r0 * np.sinh(t), r0 * np.cosh(t)))
    if branch == "bottom":
        return np.vstack((r0 * np.sinh(t), -r0 * np.cosh(t)))
    raise ValueError(f"unknown branch {branch!r}")


def _normal_outward(points: np.ndarray, sign: float) -> np.ndarray:
    x, y = points
    raw = sign * np.vstack((x, -y))
    return raw / np.linalg.norm(raw, axis=0)


def truncated_electrode_polygon(
    branch: str, r0: float = 1.0, thickness_over_r0: float = 0.1,
    extent_over_r0: float = 1.5, n_curve: int = 100,
) -> np.ndarray:
    """A finite strip following a hyperbola, with semicircular end caps.

    The central-facing boundary is exactly the required hyperbola in the
    parameterization.  Backing thickness is measured along the outward normal.
    The cap radius is half the strip thickness. It regularizes the artificial
    in-plane termination; it does not represent 3-D electrode end rounding.
    """
    if r0 <= 0 or thickness_over_r0 <= 0 or extent_over_r0 <= 1:
        raise ValueError("r0/thickness must be positive and extent_over_r0 > 1")
    if n_curve < 12:
        raise ValueError("n_curve must be >= 12")
    T = asinh(extent_over_r0)
    t = np.linspace(-T, T, n_curve)
    face = hyperbola_branch(t, branch, r0)
    # +grad(x²-y²) points outward from the central aperture on left/right;
    # -grad(x²-y²) does so on top/bottom.
    outward_sign = +1.0 if branch in {"right", "left"} else -1.0
    normals = _normal_outward(face, outward_sign)
    thick = thickness_over_r0 * r0
    back = face + thick * normals

    # Rounded end caps: a semicircle in the local tangent/normal frame.
    def unit(v):
        return v / np.linalg.norm(v)
    tangent_start = unit(face[:, 1] - face[:, 0])
    tangent_end = unit(face[:, -1] - face[:, -2])
    n_start = unit(back[:, 0] - face[:, 0])
    n_end = unit(back[:, -1] - face[:, -1])
    radius = thick / 2.0
    mid_start = 0.5 * (face[:, 0] + back[:, 0])
    mid_end = 0.5 * (face[:, -1] + back[:, -1])
    theta_end = np.linspace(pi, 0.0, 13)
    end_arc = mid_end[:, None] + radius * (
        n_end[:, None] * np.cos(theta_end) + tangent_end[:, None] * np.sin(theta_end)
    )
    theta_start = np.linspace(0.0, -pi, 13)
    start_arc = mid_start[:, None] + radius * (
        n_start[:, None] * np.cos(theta_start) + tangent_start[:, None] * np.sin(theta_start)
    )
    polygon = np.hstack((face, end_arc[:, 1:-1], back[:, ::-1], start_arc[:, 1:-1]))
    area2 = np.sum(polygon[0] * np.roll(polygon[1], -1) - polygon[1] * np.roll(polygon[0], -1))
    if area2 < 0:
        polygon = polygon[:, ::-1]
    return polygon


def ideal_potential_v(x: np.ndarray, y: np.ndarray, v_differential_peak_v: float = 1.0,
                      r0: float = 1.0) -> np.ndarray:
    """Ideal field for electrode potentials +/-Vdiff/2."""
    return v_differential_peak_v * (x*x - y*y) / (2.0*r0*r0)


def _add_polygon_loop(gmsh_model, polygon: np.ndarray, point_size: float) -> tuple[int, list[int]]:
    tags = [gmsh_model.geo.addPoint(float(x), float(y), 0.0, point_size) for x, y in polygon.T]
    curves = [gmsh_model.geo.addLine(tags[i], tags[(i+1) % len(tags)]) for i in range(len(tags))]
    # A conductor is a hole in the vacuum domain, so its loop orientation is
    # opposite to the counter-clockwise outer boundary.
    return gmsh_model.geo.addCurveLoop([-c for c in curves[::-1]]), curves


def _make_mesh(hmin: float, r0: float, finite_electrodes: bool,
               thickness_over_r0: float = 0.1, extent_over_r0: float = 1.5,
               outer_over_r0: float = 6.0, n_curve: int = 100):
    gmsh = _gmsh_module()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.add("hyperbolic_trap_2d")
    geo = gmsh.model.geo
    outer = outer_over_r0*r0
    corners = [geo.addPoint(-outer, -outer, 0, 0.5*r0), geo.addPoint(outer, -outer, 0, 0.5*r0),
               geo.addPoint(outer, outer, 0, 0.5*r0), geo.addPoint(-outer, outer, 0, 0.5*r0)]
    outer_curves = [geo.addLine(corners[i], corners[(i+1) % 4]) for i in range(4)]
    outer_loop = geo.addCurveLoop(outer_curves)
    electrode_curves: dict[str, list[int]] = {}
    holes = []
    if finite_electrodes:
        for branch in ("right", "left", "top", "bottom"):
            polygon = truncated_electrode_polygon(branch, r0, thickness_over_r0,
                                                   extent_over_r0, n_curve)
            loop, curves = _add_polygon_loop(gmsh.model, polygon, hmin*r0)
            holes.append(loop)
            electrode_curves[branch] = curves
    surface = geo.addPlaneSurface([outer_loop, *holes])
    geo.synchronize()
    phys_outer = gmsh.model.addPhysicalGroup(1, outer_curves, tag=1)
    gmsh.model.setPhysicalName(1, phys_outer, "outer_ground")
    phys_electrodes: dict[str, int] = {}
    if finite_electrodes:
        for k, branch in enumerate(("right", "left", "top", "bottom"), start=2):
            tag = gmsh.model.addPhysicalGroup(1, electrode_curves[branch], tag=k)
            gmsh.model.setPhysicalName(1, tag, f"electrode_{branch}")
            phys_electrodes[branch] = tag
    gmsh.model.addPhysicalGroup(2, [surface], tag=1)

    # Two overlapping threshold fields give adaptive refinement near the
    # electrode-facing surfaces and the central ion-sensitive volume.
    fields = []
    electrode_entities = [c for cs in electrode_curves.values() for c in cs]
    if electrode_entities:
        fdist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fdist, "CurvesList", electrode_entities)
        fth = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(fth, "InField", fdist)
        gmsh.model.mesh.field.setNumber(fth, "SizeMin", hmin*r0)
        gmsh.model.mesh.field.setNumber(fth, "SizeMax", 0.25*r0)
        gmsh.model.mesh.field.setNumber(fth, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(fth, "DistMax", 0.75*r0)
        fields.append(fth)
    center_point = geo.addPoint(0, 0, 0, hmin*r0)
    geo.synchronize()
    fdistc = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fdistc, "PointsList", [center_point])
    fthc = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fthc, "InField", fdistc)
    gmsh.model.mesh.field.setNumber(fthc, "SizeMin", hmin*r0)
    gmsh.model.mesh.field.setNumber(fthc, "SizeMax", 0.30*r0)
    gmsh.model.mesh.field.setNumber(fthc, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(fthc, "DistMax", 1.5*r0)
    fields.append(fthc)
    if fields:
        fmin = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hmin*r0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4*r0)
    gmsh.model.mesh.generate(2)

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    points = np.asarray(coords).reshape(-1, 3)[:, :2]
    tag_to_idx = {int(tag): i for i, tag in enumerate(node_tags)}
    element_types, element_tags, element_nodes_by_type = gmsh.model.mesh.getElements(2)
    triangle_blocks = [nodes.reshape(-1, 3) for typ, nodes in zip(element_types, element_nodes_by_type) if typ == 2]
    if not triangle_blocks:
        raise RuntimeError(f"No 3-node triangles in mesh; Gmsh element types={element_types}")
    gmsh_nodes = np.vstack(triangle_blocks)
    triangles = np.array([[tag_to_idx[int(t)] for t in row] for row in gmsh_nodes], dtype=int)
    boundary_nodes: dict[str, np.ndarray] = {}
    outer_tags = gmsh.model.mesh.getNodesForPhysicalGroup(1, phys_outer)[0]
    boundary_nodes["outer"] = np.array([tag_to_idx[int(t)] for t in outer_tags], dtype=int)
    for branch, tag in phys_electrodes.items():
        tags = gmsh.model.mesh.getNodesForPhysicalGroup(1, tag)[0]
        boundary_nodes[branch] = np.array([tag_to_idx[int(t)] for t in tags], dtype=int)
    gmsh.finalize()
    # Standalone points used only to drive Gmsh's distance field are not part
    # of the FE space; drop any unreferenced mesh nodes to avoid null modes.
    used = np.unique(triangles)
    remap = np.full(len(points), -1, dtype=int)
    remap[used] = np.arange(len(used))
    points = points[used]
    triangles = remap[triangles]
    boundary_nodes = {key: remap[val] for key, val in boundary_nodes.items()}
    return points, triangles, boundary_nodes


def _assemble_stiffness(points: np.ndarray, triangles: np.ndarray):
    p = points[triangles]
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]
    det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    if np.any(det <= 0):
        bad = np.count_nonzero(det <= 0)
        raise RuntimeError(f"Gmsh returned {bad} non-positive-orientation triangles")
    area = det/2
    grad = np.empty((len(triangles), 3, 2))
    grad[:, 0, :] = np.column_stack((y2-y3, x3-x2)) / det[:, None]
    grad[:, 1, :] = np.column_stack((y3-y1, x1-x3)) / det[:, None]
    grad[:, 2, :] = np.column_stack((y1-y2, x2-x1)) / det[:, None]
    local = area[:, None, None] * np.einsum("eia,eja->eij", grad, grad)
    rows = np.repeat(triangles, 3, axis=1).ravel()
    cols = np.tile(triangles, (1, 3)).ravel()
    matrix = coo_matrix((local.ravel(), (rows, cols)), shape=(len(points), len(points))).tocsr()
    return matrix, area, grad


def _solve_dirichlet(matrix, boundary: dict[int, float]) -> np.ndarray:
    ids = np.fromiter(boundary.keys(), dtype=int)
    vals = np.fromiter(boundary.values(), dtype=float)
    u = np.zeros(matrix.shape[0])
    u[ids] = vals
    mask = np.ones(matrix.shape[0], dtype=bool)
    mask[ids] = False
    interior = np.flatnonzero(mask)
    u[interior] = spsolve(matrix[interior][:, interior], -matrix[interior][:, ids] @ vals)
    return u


def _fit_center_polynomial(points: np.ndarray, values: np.ndarray, rfit: float = 0.25):
    radius = np.linalg.norm(points, axis=1)
    selected = radius < rfit
    x, y, z = points[selected, 0], points[selected, 1], values[selected]
    A = np.column_stack((np.ones_like(x), x, y, 0.5*x*x, x*y, 0.5*y*y))
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    # potential: c + gx*x + gy*y + 0.5*hxx*x² + hxy*x*y + 0.5*hyy*y²
    return coef


def _fit_multipoles(points: np.ndarray, values: np.ndarray, rfit: float = 0.4):
    """Fit a local Taylor polynomial through fourth order in dimensionless XY."""
    radius = np.linalg.norm(points, axis=1)
    selected = radius < rfit
    x, y, z = points[selected, 0]/rfit, points[selected, 1]/rfit, values[selected]
    powers = [(i, n-i) for n in range(5) for i in range(n, -1, -1)]
    A = np.column_stack([x**i * y**j for i,j in powers])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    return powers, coef


def _relative_l2(points: np.ndarray, triangles: np.ndarray, values: np.ndarray,
                 exact_fn) -> float:
    p = points[triangles]
    centroids = p.mean(axis=1)
    interp = values[triangles].mean(axis=1)
    exact = exact_fn(centroids[:, 0], centroids[:, 1])
    d1, d2 = p[:, 1]-p[:, 0], p[:, 2]-p[:, 0]
    w = np.abs(d1[:, 0]*d2[:, 1] - d1[:, 1]*d2[:, 0])/2
    return float(np.sqrt(np.sum(w*(interp-exact)**2)/np.sum(w*exact**2)))


def _relative_field_l2(points: np.ndarray, triangles: np.ndarray, values: np.ndarray,
                       basis_gradients: np.ndarray) -> float:
    """Relative L2 electric-field error for the exact unit quadrupole reference."""
    p = points[triangles]
    grad_u = np.einsum("eia,ei->ea", basis_gradients, values[triangles])
    # E=-grad(phi); exact phi=(x²-y²)/2 so grad(phi)=(x,-y).
    # Integrate |grad_h-grad_exact|² exactly over each triangle.  A centroid-
    # only rule misses the within-element variance of the linear exact field.
    centroids = p.mean(axis=1)
    exact_grad = np.column_stack((centroids[:,0], -centroids[:,1]))
    d1, d2 = p[:,1]-p[:,0], p[:,2]-p[:,0]
    area = np.abs(d1[:,0]*d2[:,1]-d1[:,1]*d2[:,0])/2
    centered = p-centroids[:,None,:]
    trace_cov = np.sum(centered**2,axis=(1,2))/12.0
    err = np.sum(area*(np.sum((grad_u-exact_grad)**2,axis=1)+trace_cov))
    norm = np.sum(area*(np.sum(exact_grad**2,axis=1)+trace_cov))
    return float(np.sqrt(err/norm))


def solve_mesh_level(level: int, hmin_over_r0: float,
                     vdiff_v: float = 10.0, thickness_over_r0: float = 0.1,
                     extent_over_r0: float = 1.5, outer_over_r0: float = 6.0,
                     n_curve: int = 100) -> tuple[MeshResult, dict[str, np.ndarray]]:
    """Run ideal reference and finite-device drive/weighting solves at one h level."""
    hdim = hmin_over_r0
    points, triangles, boundaries = _make_mesh(hdim, 1.0, True, thickness_over_r0,
                                               extent_over_r0, outer_over_r0, n_curve)
    K, areas, grads = _assemble_stiffness(points, triangles)
    allb = np.concatenate(list(boundaries.values()))
    # Ideal reference: exact potential on an outer square, no internal conductors.
    # Re-use finite mesh domain but ideal quadrupole data cannot coexist with holes;
    # a separate rectangular mesh is generated below for a clean Laplace benchmark.
    p0, t0, b0 = _make_rect_reference(hdim, outer_over_r0)
    K0, area0, grad0 = _assemble_stiffness(p0, t0)
    ref_b = {int(i): float(ideal_potential_v(p0[i,0], p0[i,1], 1.0, 1.0)) for i in b0["outer"]}
    uref = _solve_dirichlet(K0, ref_b)
    relpot = _relative_l2(p0, t0, uref, lambda x,y: ideal_potential_v(x,y,1.0,1.0))
    relfield = _relative_field_l2(p0,t0,uref,grad0)
    coef_ref = _fit_center_polynomial(p0, uref)
    hessian_ref = np.array([coef_ref[3], coef_ref[4], coef_ref[5]])
    exact_hess = np.array([1.0, 0.0, -1.0])
    relhess = float(np.linalg.norm(hessian_ref-exact_hess)/np.linalg.norm(exact_hess))

    drive_bc: dict[int,float] = {int(i): 0.0 for i in boundaries["outer"]}
    for branch in ("right", "left"):
        drive_bc.update({int(i): +0.5*vdiff_v for i in boundaries[branch]})
    for branch in ("top", "bottom"):
        drive_bc.update({int(i): -0.5*vdiff_v for i in boundaries[branch]})
    drive = _solve_dirichlet(K, drive_bc)
    cd = _fit_center_polynomial(points, drive)
    hxx, hxy, hyy = cd[3:6]
    qscale = 0.5*(hxx-hyy)
    residual = float(np.sqrt((0.5*(hxx+hyy))**2 + hxy**2) / max(abs(qscale), 1e-30))
    center_gradient_fraction = float(np.linalg.norm(cd[1:3]) / max(abs(qscale)*0.25, 1e-30))
    powers, mc = _fit_multipoles(points, drive)
    high = np.array([c for (i,j),c in zip(powers,mc) if i+j >= 3])
    multipole_fraction = float(np.linalg.norm(high) / max(abs(qscale)*0.4**2,1e-30))

    weights = {}
    wgrads = {}
    weight_element_gradients = {}
    caps = {}
    for branch in ("right", "left", "top", "bottom"):
        bc = {int(i): 0.0 for i in boundaries["outer"]}
        for other in ("right", "left", "top", "bottom"):
            val = 1.0 if other == branch else 0.0
            bc.update({int(i): val for i in boundaries[other]})
        w = _solve_dirichlet(K, bc)
        weights[branch] = w
        fit = _fit_center_polynomial(points, w)
        wgrads[branch] = np.array(fit[1:3])
        # C' = eps0 * integral |grad Vw|² dA for a 1-V terminal.
        grad_u = np.einsum("eia,ei->ea", grads, w[triangles])
        weight_element_gradients[branch] = grad_u
        caps[branch] = EPSILON_0_F_PER_M * float(np.sum(areas*np.sum(grad_u**2, axis=1)))
    capacitance_matrix = {
        (a,b): EPSILON_0_F_PER_M*float(np.sum(
            areas*np.sum(weight_element_gradients[a]*weight_element_gradients[b],axis=1)))
        for a in weights for b in weights
    }
    c_diff = 0.25*(capacitance_matrix[("right","right")]
                   +capacitance_matrix[("left","left")]
                   -2*capacitance_matrix[("right","left")])
    pair_gradient = wgrads["right"] + wgrads["left"]
    # One volt across the differential pair means +/-0.5 V on its terminals.
    differential_gradient = 0.5*(wgrads["right"] - wgrads["left"])
    single_grad = float(np.linalg.norm(wgrads["right"]))
    paired_grad = float(np.linalg.norm(pair_gradient))
    diff_grad = float(np.linalg.norm(differential_gradient))
    result = MeshResult(
        level=level, h_min_over_r0=hdim, nodes=len(points), triangles=len(triangles),
        ideal_relative_potential_l2=relpot, ideal_relative_field_l2=relfield,
        ideal_relative_hessian=relhess,
        drive_quadrupole_hessian_v_per_r0_2=float(qscale),
        drive_center_gradient_fraction=center_gradient_fraction,
        drive_quadrupole_residual_fraction=residual,
        drive_higher_multipole_fraction_at_0_4r0=multipole_fraction,
        single_electrode_weighting_gradient_per_r0=single_grad,
        paired_rf_weighting_gradient_per_r0=paired_grad,
        differential_weighting_gradient_per_r0=diff_grad,
        per_electrode_capacitance_f_per_m=float(np.mean(list(caps.values()))),
        differential_port_capacitance_f_per_m=float(c_diff),
    )
    fields = {"points": points, "triangles": triangles, "drive": drive, **{f"weight_{k}": v for k,v in weights.items()}}
    return result, fields


def _make_rect_reference(hmin: float, outer_over_r0: float):
    """Gmsh rectangular ideal-validation domain with locally refined center."""
    gmsh = _gmsh_module()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.add("ideal_quadrupole_validation")
    geo = gmsh.model.geo
    o = outer_over_r0
    p = [geo.addPoint(-o,-o,0,0.5), geo.addPoint(o,-o,0,0.5), geo.addPoint(o,o,0,0.5), geo.addPoint(-o,o,0,0.5)]
    lines = [geo.addLine(p[i],p[(i+1)%4]) for i in range(4)]
    loop=geo.addCurveLoop(lines); surf=geo.addPlaneSurface([loop]); geo.synchronize()
    outer=gmsh.model.addPhysicalGroup(1,lines,tag=1)
    gmsh.model.addPhysicalGroup(2,[surf],tag=1)
    center=geo.addPoint(0,0,0,hmin); geo.synchronize()
    d=gmsh.model.mesh.field.add("Distance"); gmsh.model.mesh.field.setNumbers(d,"PointsList",[center])
    th=gmsh.model.mesh.field.add("Threshold"); gmsh.model.mesh.field.setNumber(th,"InField",d)
    gmsh.model.mesh.field.setNumber(th,"SizeMin",hmin); gmsh.model.mesh.field.setNumber(th,"SizeMax",0.3)
    gmsh.model.mesh.field.setNumber(th,"DistMin",0); gmsh.model.mesh.field.setNumber(th,"DistMax",1.5)
    gmsh.model.mesh.field.setAsBackgroundMesh(th)
    gmsh.option.setNumber("Mesh.MeshSizeMin",hmin); gmsh.option.setNumber("Mesh.MeshSizeMax",0.4)
    gmsh.model.mesh.generate(2)
    nt,xyz,_=gmsh.model.mesh.getNodes(); pts=np.asarray(xyz).reshape(-1,3)[:,:2]
    mp={int(tag):i for i,tag in enumerate(nt)}
    element_types, element_tags, node_blocks = gmsh.model.mesh.getElements(2)
    triangle_blocks = [nodes.reshape(-1, 3) for typ, nodes in zip(element_types,node_blocks) if typ == 2]
    if not triangle_blocks:
        raise RuntimeError(f"No ideal-reference triangles; Gmsh element types={element_types}")
    tri=np.array([[mp[int(t)] for t in row] for row in np.vstack(triangle_blocks)],int)
    btags=gmsh.model.mesh.getNodesForPhysicalGroup(1,outer)[0]
    bn=np.array([mp[int(t)] for t in btags],int)
    gmsh.finalize()
    used=np.unique(tri); remap=np.full(len(pts),-1,dtype=int); remap[used]=np.arange(len(used))
    pts=pts[used]; tri=remap[tri]; bn=remap[bn]
    return pts,tri,{"outer":bn}


def image_current_peak_a(species: IonSpecies, omega_rad_s: float, amplitude_m: float,
                         weighting_gradient_per_m: float) -> float:
    """Shockley-Ramo peak current magnitude q*omega*A*|grad(phi_w)|."""
    if omega_rad_s < 0 or amplitude_m < 0 or weighting_gradient_per_m < 0:
        raise ValueError("omega, amplitude, and weighting gradient must be nonnegative")
    return abs(species.charge_c) * omega_rad_s * amplitude_m * weighting_gradient_per_m


def parallel_rc_response(signal_current_a: float, frequency_hz: float, resistance_ohm: float,
                         capacitance_f: float, temperature_k: float, enbw_hz: float,
                         amp_voltage_noise_v_rt_hz: float = 0.0,
                         amp_current_noise_a_rt_hz: float = 0.0) -> dict[str, float]:
    """Unmatched Norton-current readout into R || C with input-referred white noise."""
    if min(frequency_hz, resistance_ohm, capacitance_f, temperature_k, enbw_hz) <= 0:
        raise ValueError("frequency, R, C, T, and ENBW must be positive")
    w=2*pi*frequency_hz
    z=1.0/complex(1.0/resistance_ohm,w*capacitance_f)
    signal_v=abs(signal_current_a*z)
    # Resistor current noise flows through exactly the same impedance;
    # integrate equivalent white source noise across the declared ENBW.
    thermal_i=sqrt(4*K_B_J_PER_K*temperature_k/resistance_ohm*enbw_hz)
    amp_i=amp_current_noise_a_rt_hz*sqrt(enbw_hz)
    amp_v=amp_voltage_noise_v_rt_hz*sqrt(enbw_hz)
    noise_v=sqrt((thermal_i*abs(z))**2+(amp_i*abs(z))**2+amp_v**2)
    return {"z_signal_ohm":abs(z),"signal_v_rms_or_peak_as_input":signal_v,
            "thermal_current_noise_a_rms":thermal_i,"noise_v_rms":noise_v,
            "amplitude_ratio":signal_v/noise_v if noise_v else float("inf"),
            "pole_hz":1/(2*pi*resistance_ohm*capacitance_f)}


def parallel_rp_resonator_response(signal_current_a: float, frequency_hz: float,
                                   parallel_resistance_ohm: float, capacitance_f: float,
                                   temperature_k: float, enbw_hz: float,
                                   amp_voltage_noise_v_rt_hz: float = 0.0,
                                   amp_current_noise_a_rt_hz: float = 0.0) -> dict[str, float]:
    """Parallel resonator sensitivity at resonance; loss represented by Rp.

    Requires a tuning inductance L=1/(omega0² C). Parasitic C is included in
    C and therefore changes both required L and Q.  This remains an equivalent
    lumped model, not an impedance-matched or fabricated receiver.
    """
    if min(frequency_hz, parallel_resistance_ohm, capacitance_f, temperature_k, enbw_hz)<=0:
        raise ValueError("frequency, Rp, C, T, and ENBW must be positive")
    w=2*pi*frequency_hz
    q=parallel_resistance_ohm*w*capacitance_f
    if enbw_hz > frequency_hz/max(q,1e-30):
        raise ValueError("declared ENBW exceeds resonator 3-dB bandwidth in this model")
    z=parallel_resistance_ohm
    signal_v=abs(signal_current_a*z)
    # Rectangular post-detection passband centered on resonance. Integrate the
    # local Lorentzian |Z(f)|^2 instead of treating the resonator as flat.
    x=q*enbw_hz/frequency_hz
    transfer_factor=float(np.arctan(x)/x) if x else 1.0
    noise_bandwidth=enbw_hz*transfer_factor
    thermal_i=sqrt(4*K_B_J_PER_K*temperature_k/parallel_resistance_ohm*noise_bandwidth)
    amp_i=amp_current_noise_a_rt_hz*sqrt(noise_bandwidth)
    amp_v=amp_voltage_noise_v_rt_hz*sqrt(enbw_hz)
    noise=sqrt((thermal_i*z)**2+(amp_i*z)**2+amp_v**2)
    return {"q_loaded_lumped":q,"required_inductance_h":1/(w*w*capacitance_f),
            "resonator_bandwidth_hz":frequency_hz/max(q,1e-30),
            "postfilter_lorentzian_noise_factor":transfer_factor,
            "signal_v_peak":signal_v,"noise_v_rms":noise,
            "amplitude_ratio":signal_v/noise if noise else float("inf")}


def event_bin_probability(rate_per_s: float, bin_width_s: float) -> tuple[float,float]:
    """Exact Poisson probability of >=1 event and ideal binary Hmin per bin."""
    if rate_per_s < 0 or bin_width_s <= 0:
        raise ValueError("rate must be nonnegative and interval positive")
    p_event=-np.expm1(-rate_per_s*bin_width_s)
    p0=1.0-p_event
    pmax=max(p0,p_event)
    hmin=-np.log2(pmax) if pmax>0 else 0.0
    return float(p_event),float(hmin)
