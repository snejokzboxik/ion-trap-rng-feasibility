"""Targeted P2 solves for already-probed M3a.2 center-refined meshes."""
from __future__ import annotations
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis,MeshTet,ElementTetP2,asm
from skfem.models.poisson import laplace
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from run_m2b_3d_preflight import boundary_values
from iontrap_rng.hessian import fit_local_hessian

def jsonable(x):
    if isinstance(x,dict): return {str(k):jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [jsonable(v) for v in x]
    if isinstance(x,np.generic): return x.item()
    return x

def solve(hc):
    tag=str(hc).replace('.','p')
    meshfile=ROOT/'results'/f'm3a2_center_mesh_hc{tag}.npz'
    saved=np.load(meshfile); p=saved['points']; t=saved['tetrahedra']
    mesh=MeshTet(p.T,t.T); basis=Basis(mesh,ElementTetP2(),intorder=6)
    start=time.perf_counter(); A=asm(laplace,basis).tocsc(); assemble_s=time.perf_counter()-start
    boundary=np.unique(basis.get_dofs(facets=mesh.boundary_facets()).all())
    mask=np.zeros(basis.N,dtype=bool); mask[boundary]=True
    bc0=boundary_values(basis.doflocs.T,p2_surface_tolerance=.02)
    bc={k:np.asarray(v,int)[mask[np.asarray(v,int)]] for k,v in bc0.items()}
    terms=('right','left','top','bottom','dc_plus','dc_minus')
    # The grounded enclosure is a fixed zero-potential boundary too.
    fixed=np.unique(np.concatenate([np.asarray(v,int) for v in bc.values()])); missing=np.setdiff1d(boundary,fixed)
    if len(missing): raise RuntimeError(f'unassigned boundary dofs: {len(missing)}')
    overlaps=[]
    for i,a in enumerate(terms):
        for b in terms[i+1:]:
            n=len(np.intersect1d(bc[a],bc[b]))
            if n: overlaps.append({'terminals':[a,b],'count':int(n)})
    free=np.setdiff1d(np.arange(basis.N),fixed)
    start=time.perf_counter(); lu=splu(A[free][:,free]); factor_s=time.perf_counter()-start
    fields={}; solve_start=time.perf_counter()
    for key in terms:
        v=np.zeros(basis.N); v[bc[key]]=1.; v[free]=lu.solve(-(A[free][:,fixed]@v[fixed])); fields[key]=v
    solve_s=time.perf_counter()-solve_start
    rf=.5*(fields['right']+fields['left']-fields['top']-fields['bottom'])
    dc=fields['dc_plus']+fields['dc_minus']
    fitrows={}
    for radius in (.35,.50):
        fitrows[str(radius)]={}
        for name,u in (('rf',rf),('dc',dc)):
            f=fit_local_hessian(basis.doflocs.T,u,radius=radius,order=2,harmonic_quadratic=False)
            fitrows[str(radius)][name]={'samples':f.samples,'rank':f.rank,'condition':f.condition_number,
                'relative_potential_residual':f.relative_residual,'hessian_V_per_r0_squared':f.hessian.tolist(),
                'trace_V_per_r0_squared':float(np.trace(f.hessian)),'gradient_at_origin_V_per_r0':f.gradient.tolist()}
    output={'hcenter_target_r0':float(hc),'mesh_file':str(meshfile.relative_to(ROOT)),
      'mesh_probe':json.loads((ROOT/'results'/f'm3a2_center_mesh_hc{tag}.json').read_text()),
      'geometry':'unchanged M2b.2 four hyperbolic RF branches, finite length 10r0, annular DC rings, empty vacuum, ideal conductors, polygonal boundary',
      'normalization':{'rf':'1 V differential peak pattern','dc':'1 V on each ring','coordinates':'r/r0'},
      'P2_dofs':int(basis.N),'fixed_dofs':int(len(fixed)),'missing_boundary_dofs':int(len(missing)),'terminal_boundary_dof_counts':{k:int(len(bc[k])) for k in terms},
      'terminal_tag_overlaps':overlaps,'matrix_nnz':int(A.nnz),'assembly_s':assemble_s,'factorization_s':factor_s,'six_rhs_s':solve_s,
      'local_nodal_potential_hessian_fits':fitrows,
      'limitations':['linear polygonal geometry; perfect conductors; vacuum; no dielectric stack','P2 gradients are elementwise discontinuous','mesh-only central refinement; RF/DC terminal geometry unchanged']}
    np.savez_compressed(ROOT/'results'/f'm3a2_center_fields_hc{tag}.npz',points=p,tetrahedra=t,
      p2_dof_coordinates=basis.doflocs.T,rf_potential=rf,dc_potential=dc)
    (ROOT/'results'/f'm3a2_center_fields_hc{tag}.json').write_text(json.dumps(jsonable(output),indent=2),encoding='utf-8')
    print(json.dumps(jsonable({'hc':hc,'dofs':basis.N,'nnz':A.nnz,'assembly_s':assemble_s,'factorization_s':factor_s,'rhs_s':solve_s,'hessian_fits':fitrows}),indent=2),flush=True)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--hcenter',type=float,required=True);args=ap.parse_args();solve(args.hcenter)
