"""Cross-check DC Hessians from saved targeted-refinement P2 nodal fields."""
import json
from pathlib import Path
import numpy as np
from skfem import Basis,MeshTet,ElementTetP2
ROOT=Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.hessian import fit_local_hessian

def direct_element_fit(basis,u,mesh):
    centers=mesh.p[:,mesh.t].mean(axis=1).T;ie=int(np.argmin(np.linalg.norm(centers,axis=1)))
    d=basis.element_dofs[:,ie];p=basis.doflocs[:,d].T;v=u[d];c=p.mean(axis=0)
    s=max(np.max(np.linalg.norm(p-c,axis=1)),1e-15);q=(p-c)/s
    M=np.column_stack((np.ones(len(q)),q[:,0],q[:,1],q[:,2],.5*q[:,0]**2,.5*q[:,1]**2,.5*q[:,2]**2,q[:,0]*q[:,1],q[:,0]*q[:,2],q[:,1]*q[:,2]))
    coef,_,rank,sing=np.linalg.lstsq(M,v,rcond=None)
    h=np.array([[coef[4],coef[7],coef[8]],[coef[7],coef[5],coef[9]],[coef[8],coef[9],coef[6]]])/s**2
    return {'element':ie,'rank':int(rank),'condition':float(sing[0]/sing[-1]),'Hessian':h.tolist(),'trace':float(np.trace(h))}

def analyze(filename):
    D=np.load(filename);p=D['points'];t=D['tetrahedra'];mesh=MeshTet(p.T,t.T);basis=Basis(mesh,ElementTetP2())
    if not np.allclose(basis.doflocs.T,D['p2_dof_coordinates'],rtol=0,atol=1e-13):raise RuntimeError('P2 DOF mapping changed')
    u=D['dc_potential'];fits={}
    for r in (.25,.35,.5,.75):
        f=fit_local_hessian(basis.doflocs.T,u,radius=r,order=2)
        h=fit_local_hessian(basis.doflocs.T,u,radius=r,order=2,harmonic_quadratic=True)
        fits[str(r)]={'unconstrained':{'samples':f.samples,'rank':f.rank,'condition':f.condition_number,'residual':f.relative_residual,'Hessian':f.hessian.tolist(),'trace':float(np.trace(f.hessian))},
                      'harmonic_constrained':{'samples':h.samples,'rank':h.rank,'condition':h.condition_number,'residual':h.relative_residual,'Hessian':h.hessian.tolist(),'trace':float(np.trace(h.hessian))}}
    # Pointwise potential symmetry audit on deterministic interior samples.
    # The tetra element finder scales poorly with points × cells; use a small,
    # fixed interior set to keep the audit below a few hundred MB locally.
    rng=np.random.default_rng(20261027);q=rng.normal(size=(32,3));q=q/np.linalg.norm(q,axis=1)[:,None]*rng.uniform(0,.45,(32,1))
    interp=basis.interpolator(u);base=interp(q.T);sym={}
    for ax in range(3):
        qm=q.copy();qm[:,ax]*=-1;other=interp(qm.T)
        sym['xyz'[ax]]=float(np.sqrt(np.mean((base-other)**2))/max(np.sqrt(np.mean(base**2)),1e-300))
    return {'file':Path(filename).name,'P2_dofs':int(basis.N),'direct_central_element_P2_derivative':direct_element_fit(basis,u,mesh),
            'regional_potential_value_fits':fits,'potential_reflection_residual':sym}

if __name__=='__main__':
    files=[ROOT/'results'/'m2b2_dc_surface_refined_p2_fields_h0p1.npz']
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('files',nargs='*');args=ap.parse_args()
    target=args.files or [str(f) for f in files]
    out={'analyses':[analyze(f) for f in target],
         'note':'Local element derivative and regional fits are independent postprocessings of the same P2 solutions; they are not independent mesh-error estimators.'}
    tag='_'.join(Path(f).stem.replace('m2b2_dc_surface_refined_p2_fields','').strip('_') for f in target)
    dest=ROOT/'results'/f'm2b2_dc_saved_p2_crosscheck_{tag or "default"}.json'
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
