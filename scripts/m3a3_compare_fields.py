"""Compare full-domain and symmetry-reduced octant P2 scalar fields."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from skfem import Basis,ElementTetP2,MeshTet
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
E_CHARGE=1.602176634e-19; R0_M=1e-6; VDC=.08697882845042744; VRF_DIFF=10.

def reflect_eval(field,xyz):
    p=np.asarray(xyz,float); signs=np.where(p<0,-1.,1.); pp=np.abs(p)
    phi=field.potential(pp.T); grad=field.gradient(pp.T)*signs.T
    return phi,grad

def load_oct(path,name):
    d=np.load(path);mesh=MeshTet(d['points'].T,d['tetrahedra'].T);basis=Basis(mesh,ElementTetP2())
    return basis,d[name]

def eval_p2(basis,u,pts):
    finder=basis.mesh.element_finder(mapping=basis.mapping);q=np.asarray(pts,float).T
    cells=np.asarray(finder(*q),int)
    if np.any(cells<0):raise ValueError('point outside mesh')
    X=basis.mapping.invF(q[:,:,None],tind=cells);ed=basis.element_dofs[:,cells];co=u[ed]
    val=np.zeros(q.shape[1]);grad=np.zeros((3,q.shape[1]))
    for k in range(basis.Nbfun):
        val+=co[k]*basis.elem.lbasis(X[:,:,0],k)[0]
        grad+=co[k][None,:]*basis.elem.gbasis(basis.mapping,X,k,tind=cells)[0].grad[:,:,0]
    return val,grad

def sample_points():
    out=[]
    for axis in range(3):
      for a in (.05,.10,.20):
       for s in (-1.,1.):
        p=np.zeros(3);p[axis]=s*a;out.append(p)
    out.extend([np.array(v,float) for v in ((.07,.11,.13),(.14,.08,.17),(.18,.12,.09),(.05,.19,.16),(.22,.06,.12),(.11,.21,.04))])
    return np.asarray(out)

def compare_fields():
    pts=sample_points();records=[]
    full_paths={'coarse':'results/m3a2_center_fields_hc0p09.npz','fine':'results/m3a2_center_fields_hc0p06.npz'}
    all_oct={tag:{name:load_oct(ROOT/'results'/f'm3a3_octant_fields_{tag}.npz',nm) for name,nm in [('dc','phi_dc'),('rf','phi_rf')]} for tag in ('plateau030','plateau018')}
    # Full npz names are dc_potential/rf_potential; use the same P2 basis API.
    def full_field(path,key):
        d=np.load(ROOT/path); mesh=MeshTet(d['points'].T,d['tetrahedra'].T); b=Basis(mesh,ElementTetP2())
        return b,d[key]
    ff={lv:{kind:full_field(path,kind+'_potential') for kind in ('dc','rf')} for lv,path in full_paths.items()}
    factors={'dc':VDC,'rf':VRF_DIFF}
    for lv,by_kind in ff.items():
      for kind in ('dc','rf'):
       bf,uf=by_kind[kind]
       _,gf=eval_p2(bf,uf,pts)
       for tag in ('plateau030','plateau018'):
        bo,uo=all_oct[tag][kind];_,go=eval_p2(bo,uo,abs(pts))
        go*=np.where(pts.T<0,-1.,1.)
        scale=E_CHARGE*factors[kind]/R0_M
        Ff=-scale*gf.T;Fo=-scale*go.T;d=Fo-Ff
        for i,p in enumerate(pts):
          denom=np.linalg.norm(Ff[i]);records.append({'kind':kind,'full_level':lv,'octant_level':tag,'point_r0':p.tolist(),
            'F_full_N':Ff[i].tolist(),'F_octant_N':Fo[i].tolist(),'deltaF_N':d[i].tolist(),
            'relative_vector_pct':float(100*np.linalg.norm(d[i])/denom) if denom>1e-45 else None,
            'relative_magnitude_pct':float(100*abs(np.linalg.norm(Fo[i])-denom)/denom) if denom>1e-45 else None,
            'absolute_error_over_eV_per_r0':float(np.linalg.norm(d[i])/(E_CHARGE/R0_M))})
    summary={}
    for kind in ('dc','rf'):
      for lv in ('coarse','fine'):
       for tag in ('plateau030','plateau018'):
        a=[r for r in records if r['kind']==kind and r['full_level']==lv and r['octant_level']==tag]
        summary[f'{kind}_{lv}_{tag}']={'max_vector_pct':float(np.nanmax([r['relative_vector_pct'] or 0 for r in a])),
          'max_magnitude_pct':float(np.nanmax([r['relative_magnitude_pct'] or 0 for r in a])),
          'max_absolute_force_error_N':float(max(np.linalg.norm(r['deltaF_N']) for r in a))}
    # full-domain mirror residuals for the fixed voltage DC field.
    mirror={}
    for lv in ff:
      mirror[lv]={}
      for kind in ('dc','rf'):
       b,u=ff[lv][kind];data=[];potdiff=[]
       for a in (.05,.10,.20):
        for axis in range(3):
         p=np.zeros(3);p[axis]=a;pm=-p
         vp,gp=eval_p2(b,u,p[None,:]);vm,gm=eval_p2(b,u,pm[None,:])
         expect=gp[:,0].copy();expect[axis]*=-1
         data.append(float(np.linalg.norm(gm[:,0]-expect)/max(np.linalg.norm(gp[:,0]),1e-30)))
         potdiff.append(float(abs(vm[0]-vp[0])))
       mirror[lv][kind]={'max_axis_mirror_gradient_residual_pct':100*max(data),'all_residuals':data,
          'max_axis_mirror_potential_difference_per_V_basis':max(potdiff),'potential_differences':potdiff}
    out={'voltage_factors':factors,'force_units':'N, with dimensionless FEM coordinates xi=r/r0','points_count':len(pts),
         'records':records,'summary':summary,'full_domain_mirror_force_residuals':mirror,
         'interpretation':'octant mirror parity is imposed by construction; only full-domain residual and full-vs-octant comparisons are independent diagnostics'}
    (ROOT/'results/m3a3_full_octant_force_comparison.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2));print(json.dumps(mirror,indent=2))

if __name__=='__main__': compare_fields()
