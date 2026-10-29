"""Boundary audit and P2 electrostatics for the symmetric octant reference.

No weighting-potential problem is included: the RF/DC drive patterns are even
under all three coordinate reflections, while individual terminal excitations
do not generally share that parity.
"""
from __future__ import annotations
import json, sys, time, hashlib, os
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis, ElementTetP2, MeshTet, asm
from skfem.models.poisson import laplace
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
from m3a3_octant_mesh import truncated_electrode_polygon, RAIL_L, BOX, RING_Z, RING_T, RIN, ROUT, _stats
from matplotlib.path import Path as MPath

def classify_boundary_facets(p,t,plane_tol=2e-8):
    mesh=MeshTet(p.T,t.T)
    bf=mesh.boundary_facets(); facets=mesh.facets[:,bf].T
    xyz=p[facets]; c=xyz.mean(axis=1); x,y,z=c.T; rho=np.hypot(x,y)
    tags=np.full(len(bf),'',dtype='U16')
    # Physical conductor surfaces, including the conductor cross-sections on
    # x=0/y=0 planes, take precedence over artificial Neumann symmetry faces.
    for branch,v in [('right',.5),('top',-.5)]:
        poly=truncated_electrode_polygon(branch,1.,.10,1.5,80).T
        d=MPath(poly).contains_points(c[:,:2],radius=2e-8)
        # Extruded hyperbolic side faces, truncated rail end at z=+5, and
        # cut faces on x=0 or y=0 lying inside the conductor cross-section.
        side=np.zeros(len(c),bool)
        for a,b in zip(poly,np.roll(poly,-1,axis=0)):
            ab=b-a; tpar=np.clip(((c[:,:2]-a)*ab).sum(axis=1)/(ab@ab),0,1)
            side |= np.linalg.norm(c[:,:2]-(a+tpar[:,None]*ab),axis=1)<2e-7
        rail=(side & (z>=-plane_tol) & (z<=RAIL_L+plane_tol))
        cap=(np.abs(z-RAIL_L)<plane_tol)&d
        cut=((np.abs(x)<plane_tol)|(np.abs(y)<plane_tol))&d&(z>=-plane_tol)&(z<=RAIL_L+plane_tol)
        sel=(rail|cap|cut)&(tags=='')
        tags[sel]='rf_'+branch
    # The positive annular DC ring occupies z=5.4..5.6 and rho=.5..1.5.
    ringz=(z>=RING_Z-RING_T/2-2e-7)&(z<=RING_Z+RING_T/2+2e-7)
    dc_side=ringz & ((np.abs(rho-RIN)<.012)|(np.abs(rho-ROUT)<.012))
    dc_cap=((np.abs(z-(RING_Z-RING_T/2))<2e-7)|(np.abs(z-(RING_Z+RING_T/2))<2e-7))&(rho>=RIN-2e-7)&(rho<=ROUT+2e-7)
    dc_cut=((np.abs(x)<plane_tol)|(np.abs(y)<plane_tol))&ringz&(rho>=RIN-2e-7)&(rho<=ROUT+2e-7)
    tags[(dc_side|dc_cap|dc_cut)&(tags=='')]='dc_plus'
    # Outer grounded enclosure.
    outer=(np.abs(x-BOX[0])<plane_tol)|(np.abs(y-BOX[1])<plane_tol)|(np.abs(z-BOX[2])<plane_tol)
    tags[outer&(tags=='')]='ground'
    # Remaining exposed artificial planes are natural homogeneous Neumann.
    sym=(np.abs(x)<plane_tol)|(np.abs(y)<plane_tol)|(np.abs(z)<plane_tol)
    tags[sym&(tags=='')]='symmetry'
    unclassified=np.flatnonzero(tags=='')
    return mesh,bf,facets,tags,unclassified

def run(tag='plateau018', solve=False):
    dat=np.load(ROOT/'results'/f'm3a3_octant_mesh_{tag}.npz');p=dat['points'];t=dat['tetrahedra']
    mesh,bf,facets,tags,missing=classify_boundary_facets(p,t)
    counts={str(k):int(np.count_nonzero(tags==k)) for k in np.unique(tags)}
    audit={'mesh_tag':tag,'vertices':len(p),'tetrahedra':len(t),'boundary_facets':len(bf),'facet_counts':counts,'unclassified_facets':len(missing)}
    surface_stats={}
    for name in ('dc_plus','rf_right','rf_top','ground'):
        tri=facets[tags==name]
        if len(tri):
            es=np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]])
            es=np.unique(np.sort(es,axis=1),axis=0)
            lengths=np.linalg.norm(p[es[:,0]]-p[es[:,1]],axis=1)
            surface_stats[name]={'triangles':int(len(tri)),'unique_surface_edges_r0':_stats(lengths)}
    audit['surface_mesh_statistics']=surface_stats
    if len(missing):
        audit['unclassified_examples']=np.round(p[facets[missing[:12]]].mean(axis=1),8).tolist()
        return audit
    if not solve: return audit
    basis=Basis(mesh,ElementTetP2(),intorder=4); tic=time.perf_counter(); A=asm(laplace,basis).tocsc(); assembly=time.perf_counter()-tic
    # Interior plus open symmetry-plane DOFs are free; conductor/enclosure DOFs fixed.
    groups={k:np.flatnonzero(np.isin(bf,np.asarray(bf)[tags==k])) for k in []}
    fdict={name:bf[tags==name] for name in ('ground','rf_right','rf_top','dc_plus')}
    fixed_by={name:np.unique(basis.get_dofs(facets=ids).all()) for name,ids in fdict.items()}
    fixed=np.unique(np.concatenate(list(fixed_by.values())))
    free=np.setdiff1d(np.arange(basis.N),fixed)
    tic=time.perf_counter();lu=splu(A[free][:,free]);factor=time.perf_counter()-tic
    outs={}
    for field,values in [('rf',{'rf_right':.5,'rf_top':-.5}),('dc',{'dc_plus':1.})]:
        val=np.zeros(basis.N)
        for name,v in values.items(): val[fixed_by[name]]=v
        u=val.copy();u[free]=lu.solve(-(A[free][:,fixed]@val[fixed]))
        outs[field]=u
    payload={'points':p,'tetrahedra':t,'phi_rf':outs['rf'],'phi_dc':outs['dc']}
    np.savez_compressed(ROOT/'results'/f'm3a3_octant_fields_{tag}.npz',**payload)
    audit.update({'p2_dofs':int(basis.N),'free_dofs':int(len(free)),'fixed_dofs':int(len(fixed)),'assembly_s':float(assembly),'factorization_s':float(factor),
      'solve_s':None,'factor_nnz_LU':int(lu.L.nnz+lu.U.nnz),'solver':'scipy.sparse.linalg.splu; P2 Laplace, vacuum, natural Neumann on open octant planes',
      'RF_basis':'right=+0.5 V, top=-0.5 V; left/bottom inferred by even scalar reflection',
      'DC_basis':'positive and negative rings each +1 V; positive octant ring explicitly +1 V'})
    return audit

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--tag',default='plateau018');ap.add_argument('--solve',action='store_true');a=ap.parse_args()
    res=run(a.tag,a.solve)
    mesh_path=ROOT/'results'/f'm3a3_octant_mesh_{a.tag}.npz'
    res['input_mesh_sha256']=hashlib.sha256(mesh_path.read_bytes()).hexdigest()
    if a.solve:
        fields_path=ROOT/'results'/f'm3a3_octant_fields_{a.tag}.npz'
        res.update({'solver_completed':True,'solver_exit_code':0,
                    'config_sha256':os.environ.get('IONTRAP_RUN_CONFIG_SHA256'),
                    'solution_sha256':hashlib.sha256(fields_path.read_bytes()).hexdigest(),
                    'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
    # Keep the pre-existing mesh-only boundary audit immutable; solved-field
    # audit records are written to a distinct output name.
    audit_name=(f'm3a3_octant_solve_boundary_{a.tag}.json' if a.solve
                else f'm3a3_octant_boundary_{a.tag}.json')
    audit_path=ROOT/'results'/audit_name
    audit_path.write_text(json.dumps(res,indent=2),encoding='utf-8')
    print(json.dumps(res,indent=2))
