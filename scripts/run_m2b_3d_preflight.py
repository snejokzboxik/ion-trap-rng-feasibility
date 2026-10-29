"""Single ideal-conductor 3-D preflight for the four hyperbolic branches.

This is deliberately not a material-stack or process-qualified model. Silicon
electrodes are treated as equipotential perfect conductors; all surrounding
volume is vacuum. Results are normalized to r0 and reported as a preflight.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import gmsh
import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis, MeshTet, ElementTetP1, asm
from skfem.models.poisson import laplace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iontrap_rng.readout_fem import truncated_electrode_polygon

R0 = 1.0  # dimensionless; physical reference r0=1 um
RAIL_HALF_LENGTH = 5.0
BOX_HALF_XY = 6.0
BOX_HALF_Z = 7.0
RING_Z = 5.5
RING_THICKNESS = 0.2
RING_INNER = 0.5
RING_OUTER = 1.5


def make_mesh(hmax: float, hcenter: float, *, dc_surface_h: float | None = None,
              rf_surface_h: float | None = None, axial_end_h: float | None = None,
              return_metadata: bool = False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    requested=[hcenter]+[v for v in (dc_surface_h,rf_surface_h,axial_end_h) if v is not None]
    # The global floor must not silently suppress any requested local size.
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(requested))
    gmsh.option.setNumber("Mesh.MeshSizeMax", hmax)
    # The background fields below are the explicit sizing controls; disable
    # competing point/curvature/extension heuristics that otherwise obscure
    # actual local h. CAD curvature is still geometrically represented by OCC.
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.add("m2b_hyperbolic_vacuum_preflight")
    occ = gmsh.model.occ
    box = occ.addBox(-BOX_HALF_XY, -BOX_HALF_XY, -BOX_HALF_Z,
                     2*BOX_HALF_XY, 2*BOX_HALF_XY, 2*BOX_HALF_Z)
    solids = []
    for branch in ("right", "left", "top", "bottom"):
        poly = truncated_electrode_polygon(branch, 1.0, 0.10, 1.5, 80)
        pts = [occ.addPoint(float(x), float(y), -RAIL_HALF_LENGTH) for x, y in poly.T]
        lines = [occ.addLine(pts[i], pts[(i+1) % len(pts)]) for i in range(len(pts))]
        loop = occ.addCurveLoop(lines)
        face = occ.addPlaneSurface([loop])
        ext = occ.extrude([(2, face)], 0, 0, 2*RAIL_HALF_LENGTH)
        solids.append((3, ext[1][1]))
    # Axial DC ring electrodes at both ends, annular to retain an open axis.
    for z0 in (-RING_Z-RING_THICKNESS/2, RING_Z-RING_THICKNESS/2):
        outer = occ.addCylinder(0, 0, z0, 0, 0, RING_THICKNESS, RING_OUTER)
        inner = occ.addCylinder(0, 0, z0-0.01, 0, 0, RING_THICKNESS+0.02, RING_INNER)
        annulus, _ = occ.cut([(3, outer)], [(3, inner)], removeObject=True, removeTool=True)
        solids.extend(annulus)
    occ.cut([(3, box)], solids, removeObject=True, removeTool=True)
    occ.synchronize()
    dc_surface_entities=[];rf_surface_entities=[]
    if dc_surface_h is not None or axial_end_h is not None:
        # Select only the CAD boundaries of the two annular DC electrodes.
        # The rail ends lie at |z|=5 r0 and the DC rings at 5.4..5.6 r0.
        for dim,tag in gmsh.model.getEntities(2):
            if dim != 2: continue
            xmin,ymin,zmin,xmax,ymax,zmax=gmsh.model.getBoundingBox(dim,tag)
            near_ring=((5.35 <= zmin and zmax <= 5.65) or
                       (-5.65 <= zmin and zmax <= -5.35))
            bounded_xy=(max(abs(xmin),abs(ymin),abs(xmax),abs(ymax)) <= 1.7)
            if near_ring and bounded_xy: dc_surface_entities.append(tag)
    if rf_surface_h is not None:
        for dim,tag in gmsh.model.getEntities(2):
            if dim != 2: continue
            xmin,ymin,zmin,xmax,ymax,zmax=gmsh.model.getBoundingBox(dim,tag)
            rail_side=(zmin >= -5.05 and zmax <= 5.05 and zmax-zmin > 4.0 and
                       max(abs(xmin),abs(ymin),abs(xmax),abs(ymax)) <= 2.0)
            rail_end=((abs(zmin-5.)<.05 or abs(zmax+5.)<.05) and
                      zmax-zmin < .05 and max(abs(xmin),abs(ymin),abs(xmax),abs(ymax)) <= 2.0)
            if rail_side or rail_end: rf_surface_entities.append(tag)
    # Targeted size fields at the ion center and conductor faces.
    point = gmsh.model.occ.addPoint(0, 0, 0, hcenter)
    gmsh.model.occ.synchronize()
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "PointsList", [point])
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", dist)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", hcenter)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", hmax)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", 1.5)
    fields=[threshold]
    if dc_surface_h is not None:
        if not dc_surface_entities:
            gmsh.finalize(); raise RuntimeError("DC ring surfaces not found for targeted mesh field")
        distance_surface=gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance_surface,"SurfacesList",dc_surface_entities)
        surface_threshold=gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(surface_threshold,"InField",distance_surface)
        gmsh.model.mesh.field.setNumber(surface_threshold,"SizeMin",dc_surface_h)
        gmsh.model.mesh.field.setNumber(surface_threshold,"SizeMax",hmax)
        gmsh.model.mesh.field.setNumber(surface_threshold,"DistMin",0.)
        gmsh.model.mesh.field.setNumber(surface_threshold,"DistMax",.6)
        fields.append(surface_threshold)
    if rf_surface_h is not None:
        if not rf_surface_entities:
            gmsh.finalize(); raise RuntimeError("RF rail surfaces not found for targeted mesh field")
        dr=gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dr,"SurfacesList",rf_surface_entities)
        tr=gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(tr,"InField",dr)
        gmsh.model.mesh.field.setNumber(tr,"SizeMin",rf_surface_h)
        gmsh.model.mesh.field.setNumber(tr,"SizeMax",hmax)
        gmsh.model.mesh.field.setNumber(tr,"DistMin",0.)
        gmsh.model.mesh.field.setNumber(tr,"DistMax",.5)
        fields.append(tr)
    if axial_end_h is not None:
        ends=[]
        for zc in (-RING_Z,RING_Z):
            ends.append(gmsh.model.occ.addPoint(0,0,zc))
        gmsh.model.occ.synchronize()
        de=gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(de,"PointsList",ends)
        te=gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(te,"InField",de)
        gmsh.model.mesh.field.setNumber(te,"SizeMin",axial_end_h)
        gmsh.model.mesh.field.setNumber(te,"SizeMax",hmax)
        gmsh.model.mesh.field.setNumber(te,"DistMin",0.)
        gmsh.model.mesh.field.setNumber(te,"DistMax",1.0)
        fields.append(te)
    if len(fields)==1:
        gmsh.model.mesh.field.setAsBackgroundMesh(fields[0])
    else:
        combine=gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(combine,"FieldsList",fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(combine)
    gmsh.model.mesh.generate(3)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    p = np.asarray(coords).reshape(-1, 3)
    tets = []
    types, _, blocks = gmsh.model.mesh.getElements(3)
    tagmap = {int(t): i for i, t in enumerate(ntags)}
    for etype, block in zip(types, blocks):
        if etype == 4:
            tets.extend([[tagmap[int(t)] for t in row] for row in np.asarray(block).reshape(-1, 4)])
    gmsh.finalize()
    tets = np.asarray(tets, dtype=int)
    used = np.unique(tets)
    remap = np.full(len(p), -1, dtype=int)
    remap[used] = np.arange(len(used))
    result=(p[used], remap[tets])
    if return_metadata:
        return (*result,{"dc_surface_h_r0":dc_surface_h,"dc_surface_entity_count":len(dc_surface_entities),
                         "dc_surface_entities":dc_surface_entities,
                         "rf_surface_h_r0":rf_surface_h,"rf_surface_entity_count":len(rf_surface_entities),
                         "axial_end_h_r0":axial_end_h,
                         "mesh_size_min_effective_r0":min(requested),"mesh_size_max_r0":hmax})
    return result


def polygon_boundary_distance(xy: np.ndarray, poly: np.ndarray) -> np.ndarray:
    a = poly.T
    b = np.roll(a, -1, axis=0)
    v = b-a
    w = xy[:, None, :]-a[None, :, :]
    t = np.clip(np.sum(w*v[None,:,:], axis=2)/np.sum(v*v,axis=1)[None,:],0,1)
    d = xy[:,None,:]-(a[None,:,:]+t[:,:,None]*v[None,:,:])
    return np.sqrt(np.min(np.sum(d*d,axis=2),axis=1))


def boundary_values(p: np.ndarray, *, p2_surface_tolerance: float = 2e-6):
    tol = 1e-7
    x, y, z = p.T
    ids: dict[str, list[int]] = {k: [] for k in ("ground", "right", "left", "top", "bottom", "dc_plus", "dc_minus")}
    outer = (np.isclose(np.abs(x), BOX_HALF_XY, atol=tol) |
             np.isclose(np.abs(y), BOX_HALF_XY, atol=tol) |
             np.isclose(np.abs(z), BOX_HALF_Z, atol=tol))
    for i in np.flatnonzero(outer): ids["ground"].append(int(i))
    for branch in ("right", "left", "top", "bottom"):
        poly = truncated_electrode_polygon(branch, 1.0, 0.10, 1.5, 80)
        dist = polygon_boundary_distance(p[:,:2], poly)
        side = (dist < 2e-6) & (np.abs(z) <= RAIL_HALF_LENGTH+tol)
        # Rail end faces lie at z=+/-5 and inside the footprint polygon.
        from matplotlib.path import Path as MPath
        end = np.isclose(np.abs(z), RAIL_HALF_LENGTH, atol=tol) & MPath(poly.T).contains_points(p[:,:2])
        for i in np.flatnonzero(side | end): ids[branch].append(int(i))
    rho = np.sqrt(x*x+y*y)
    for sign, key in ((1,"dc_plus"),(-1,"dc_minus")):
        zc = sign*RING_Z
        # Annular cylindrical walls and annular end faces.
        # P2 edge nodes lie on chords of the polygonized cylindrical OCC face,
        # not exactly on its analytic radius. The caller must intersect this
        # expanded candidate set with actual boundary DOFs.
        radial_wall = (np.abs(rho-RING_INNER)<p2_surface_tolerance) | (np.abs(rho-RING_OUTER)<p2_surface_tolerance)
        face = np.isclose(np.abs(z-zc),RING_THICKNESS/2,atol=tol) & (rho>=RING_INNER-tol) & (rho<=RING_OUTER+tol)
        sel = (radial_wall & (np.abs(z-zc)<=RING_THICKNESS/2+tol)) | face
        for i in np.flatnonzero(sel): ids[key].append(int(i))
    return ids


def terminal_fields(p, tets, bc_ids):
    mesh = MeshTet(p.T, tets.T)
    basis = Basis(mesh, ElementTetP1(), intorder=4)
    A = asm(laplace, basis).tocsc()
    fixed = np.unique(np.concatenate([np.asarray(bc_ids[k],dtype=int) for k in bc_ids]))
    # Ensure every exterior boundary node is assigned a conductor or shield.
    boundary = mesh.boundary_nodes()
    missing = np.setdiff1d(boundary, fixed)
    if len(missing):
        raise RuntimeError(f"unassigned boundary nodes: {len(missing)}")
    free = np.setdiff1d(np.arange(basis.N), fixed)
    lu = splu(A[free][:,free])
    fields={}
    for key in ("right","left","top","bottom","dc_plus","dc_minus"):
        val=np.zeros(basis.N)
        val[bc_ids[key]]=1.0
        u=val.copy()
        u[free]=lu.solve(-(A[free][:,fixed]@val[fixed]))
        fields[key]=u
    return basis, fields


def center_taylor(p, u, radius=0.5):
    r = np.linalg.norm(p,axis=1)
    q = p[r<radius]
    v = u[r<radius]
    x,y,z=q.T
    M=np.column_stack((np.ones(len(x)),x,y,z,.5*x*x,.5*y*y,.5*z*z,x*y,x*z,y*z))
    c,*_=np.linalg.lstsq(M,v,rcond=None)
    return c


def run_level(hmax,hcenter,save_fields=False):
    p,t=make_mesh(hmax,hcenter)
    ids=boundary_values(p)
    terminal_count={k:len(v) for k,v in ids.items()}
    for key in ("right","left","top","bottom","dc_plus","dc_minus"):
        if terminal_count[key] < 4:
            raise RuntimeError(f"terminal {key} has only {terminal_count[key]} boundary nodes; mesh/geometry tagging failed")
    basis,fields=terminal_fields(p,t,ids)
    rf=.5*(fields["right"]+fields["left"]-fields["top"]-fields["bottom"])
    crf=center_taylor(p,rf)
    dc=fields["dc_plus"]+fields["dc_minus"]
    cdc=center_taylor(p,dc)
    wr=.5*(fields["right"]-fields["left"])
    cw=center_taylor(p,wr)
    names=("right","left","top","bottom","dc_plus","dc_minus")
    grads={k:basis.interpolate(fields[k]).grad for k in names}
    cap=np.empty((len(names),len(names)))
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            cap[i,j]=8.8541878128e-12*1e-6*float(np.sum(np.sum(grads[a]*grads[b],axis=0)*basis.dx))
    vdiff=np.array([.5,-.5,0,0,0,0])
    cdiff=float(vdiff@cap@vdiff)
    result={"hmax_r0":hmax,"hcenter_r0":hcenter,"nodes":len(p),"tetrahedra":len(t),
            "terminal_boundary_node_counts":terminal_count,
            "rf_hessian_xyz":crf[4:7].tolist(),"rf_cross_hessian_xyz":crf[7:10].tolist(),
            "axial_dc_hessian_xyz_per_V_r0_2":cdc[4:7].tolist(),
            "diff_weighting_gradient_per_r0_per_V":cw[1:4].tolist(),
            "maxwell_capacitance_matrix_F":cap.tolist(),
            "differential_capacitance_F":cdiff,
            "rf_hessian_trace_residual_fraction":float(abs(np.sum(crf[4:7]))/max(np.linalg.norm(crf[4:7]),1e-30)),
            "axial_rf_gradient_over_quadrupole_scale":float(np.linalg.norm(crf[1:4])/max(np.linalg.norm(crf[4:6]),1e-30))}
    if save_fields:
        np.savez_compressed(ROOT/"results"/"m2b_3d_fields.npz",points=p,tetrahedra=t,rf=rf,dc=dc,weighting=wr)
    return result


def main():
    rows=[run_level(h,hc,save_fields=(i==2)) for i,(h,hc) in enumerate(((0.9,0.25),(0.75,0.18),(0.65,0.13)))]
    out={"model":"ideal-conductor vacuum preflight; not a technology-qualified Si/dielectric stack",
         "geometry":{"r0_physical_m":1e-6,"hyperbolic_profile":"x^2-y^2=+/-r0^2","rail_length_r0":10,
                     "rail_thickness_r0":0.1,"branch_extent_r0":1.5,"axial_endcap":"annular rings at z=+/-5.5r0, thickness=0.2r0, radii 0.5r0 to 1.5r0",
                     "outer_ground_box_r0":[[-6,6],[-6,6],[-7,7]]},
         "assumptions":["vacuum domain only", "perfectly equipotential silicon conductors", "no oxide/support dielectric", "no RF losses or receiver coupling", "P1 tetrahedra and polygonized hyperbola, not curved high-order geometry"],
         "mesh_levels":rows}
    (ROOT/"results"/"m2b_3d_preflight.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__": main()
