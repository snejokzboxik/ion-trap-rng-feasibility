"""Bounded P2 solve on one archived M2b geometry scale; no trajectories."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis,MeshTet,ElementTetP2,asm
from skfem.models.poisson import laplace
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from run_m2b_3d_preflight import make_mesh,boundary_values

def solve_level(hmax=.9,hcenter=.25,*,mesh_npz=None,dc_surface_h=None,rf_surface_h=None,axial_end_h=None):
    mesh_meta={}
    if mesh_npz is None:
        generated=make_mesh(hmax,hcenter,dc_surface_h=dc_surface_h,rf_surface_h=rf_surface_h,
                            axial_end_h=axial_end_h,return_metadata=True)
        p,t,mesh_meta=generated
    else:
        saved=np.load(mesh_npz);p,t=saved['points'],saved['tetrahedra']
        mesh_meta={'mesh_source':str(mesh_npz),'dc_surface_h_r0':dc_surface_h,
                   'rf_surface_h_r0':rf_surface_h,'axial_end_h_r0':axial_end_h}
    mesh=MeshTet(p.T,t.T)
    basis=Basis(mesh,ElementTetP2(),intorder=6); A=asm(laplace,basis).tocsc()
    boundary=np.unique(basis.get_dofs(facets=mesh.boundary_facets()).all())
    boundary_mask=np.zeros(basis.N,dtype=bool); boundary_mask[boundary]=True
    # The initial axial rings are polygonized cylinders; chord sagitta at this
    # coarse size can approach 0.1 r0. Candidates are safely restricted to
    # boundary DOFs before application.
    bc0=boundary_values(basis.doflocs.T,p2_surface_tolerance=.02 if dc_surface_h else .1)
    bc={k:np.asarray(v,int)[boundary_mask[np.asarray(v,int)]] for k,v in bc0.items()}
    fixed=np.unique(np.concatenate([np.asarray(v,int) for v in bc.values()]))
    overlaps=[]
    keys=list(bc)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            n_overlap=len(np.intersect1d(bc[a],bc[b]))
            if n_overlap:overlaps.append({'terminals':[a,b],'count':int(n_overlap)})
    missing=np.setdiff1d(boundary,fixed)
    if len(missing):
        print('unassigned boundary examples',basis.doflocs[:,missing[:20]].T.tolist())
        print('unassigned axis/range',basis.doflocs[:,missing].min(axis=1),basis.doflocs[:,missing].max(axis=1))
        raise RuntimeError(f'{len(missing)} P2 boundary dofs unassigned')
    free=np.setdiff1d(np.arange(basis.N),fixed); lu=splu(A[free][:,free])
    terms=('right','left','top','bottom','dc_plus','dc_minus'); fields={}
    for k in terms:
        v=np.zeros(basis.N); v[bc[k]]=1.; v[free]=lu.solve(-(A[free][:,fixed]@v[fixed])); fields[k]=v
    rf=.5*(fields['right']+fields['left']-fields['top']-fields['bottom'])
    dc=fields['dc_plus']+fields['dc_minus']; weighting=.5*(fields['right']-fields['left'])
    # Quadrature-point gradient recovery: fit E/grad(phi) over centered spheres.
    X=basis.global_coordinates(); dx=basis.dx
    def recover(u,r):
        gv=basis.interpolate(u).grad
        xyz=np.moveaxis(X,0,-1).reshape(-1,3); vals=np.moveaxis(gv,0,-1).reshape(-1,3); w=dx.ravel()
        use=np.linalg.norm(xyz,axis=1)<r; q=xyz[use]/r; ww=np.sqrt(w[use]); y=vals[use]
        M=np.column_stack((np.ones(len(q)),q)); Mw=M*ww[:,None]; Yw=y*ww[:,None]
        c,_,rank,s=np.linalg.lstsq(Mw,Yw,rcond=None); H=c[1:4,:].T/r
        cond=float(s[0]/s[-1]); residual=float(np.linalg.norm(Mw@c-Yw)/np.linalg.norm(Yw))
        return {'radius_r0':r,'quadrature_samples':int(use.sum()),'rank':int(rank),'condition':cond,'fit_residual':residual,
                'gradient_at_center':c[0,:].tolist(),'Hessian':H.tolist(),'trace':float(np.trace(H))}
    # Gaussian flux balance via FEM weak form and smooth-sphere sampling is deferred;
    # this reports global PDE residual and Hessian trace recovery independently.
    rows={name:[recover(u,r) for r in (.35,.50,.75)] for name,u in [('rf',rf),('dc',dc),('weighting',weighting)]}
    names=terms; G={k:basis.interpolate(fields[k]).grad for k in names}
    C=np.empty((6,6))
    for i,a in enumerate(names):
      for j,b in enumerate(names): C[i,j]=8.8541878128e-12*1e-6*np.sum(np.sum(G[a]*G[b],axis=0)*dx)
    # P2 weak residual for each harmonic basis is zero on free rows to solver precision.
    residuals={k:float(np.linalg.norm((A@v)[free])/max(np.linalg.norm((A@v)[fixed]),1e-300)) for k,v in fields.items()}
    out={'mesh':{'hmax_r0':hmax,'hcenter_r0':hcenter,'vertices':len(p),'tetrahedra':len(t),'P2_dofs':int(basis.N),'fixed_dofs':len(fixed),'missing_boundary_dofs':len(missing),
                 'terminal_boundary_dof_counts':{k:int(len(v)) for k,v in bc.items()},'terminal_tag_overlaps':overlaps,**mesh_meta},
         'geometry':'M2b ideal-conductor hyperbolic branches + annular endcaps; vacuum; polygonal surfaces; analyst-selected supports omitted',
         'fields':rows,'maxwell_capacitance_matrix_F':C.tolist(),'weak_residual_free_dof_by_terminal':residuals,
         'limitations':['single coarse P2 geometry level; not a convergence claim','no dielectric stack or measured conductivity','no curved high-order CAD boundary','no RF port/full-wave coupling']}
    if dc_surface_h is not None:
        tag=f'{dc_surface_h:g}'.replace('.','p')
        effective='_effective' if rf_surface_h is not None or axial_end_h is not None else ''
        np.savez_compressed(ROOT/'results'/f'm2b2_dc_surface_refined_p2_fields{effective}_h{tag}.npz',
                            points=p,tetrahedra=t,p2_dof_coordinates=basis.doflocs.T,
                            rf_potential=rf,dc_potential=dc,weighting_potential=weighting)
    return out
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--hmax',type=float,default=.75); parser.add_argument('--hcenter',type=float,default=.18)
    parser.add_argument('--dc-surface-h',type=float);parser.add_argument('--rf-surface-h',type=float);parser.add_argument('--axial-end-h',type=float);parser.add_argument('--mesh-npz')
    args=parser.parse_args(); out=solve_level(args.hmax,args.hcenter,mesh_npz=args.mesh_npz,dc_surface_h=args.dc_surface_h,
                                               rf_surface_h=args.rf_surface_h,axial_end_h=args.axial_end_h)
    suffix=f'_dc{args.dc_surface_h:g}' if args.dc_surface_h else ''
    effective='_effective' if args.rf_surface_h is not None or args.axial_end_h is not None else ''
    dest=ROOT/'results'/f'm2b2_p2{effective}_h{args.hmax:g}_hc{args.hcenter:g}{suffix}.json'
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8'); print(json.dumps(out,indent=2))
