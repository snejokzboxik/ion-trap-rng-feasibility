"""Mesh convergence and harmonic-potential recovery audit for M3a.3."""
from __future__ import annotations
import json,sys,itertools
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from skfem import Basis,ElementTetP2,MeshTet
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from m3a3_compare_fields import eval_p2,sample_points,load_oct,E_CHARGE,R0_M,VDC,VRF_DIFF

def harmonic_even_basis(maxorder):
    exps=[(i,j,k) for i in range(maxorder+1) for j in range(maxorder+1-i) for k in range(maxorder+1-i-j)
          if (i+j+k)<=maxorder and all(a%2==0 for a in (i,j,k))]
    rows=sorted({(i,j,k) for a,b,c in exps for i,j,k in [(a-2,b,c),(a,b-2,c),(a,b,c-2)] if min(i,j,k)>=0})
    C=np.zeros((len(rows),len(exps)))
    for col,(a,b,c) in enumerate(exps):
      for exp,factor in [((a-2,b,c),a*(a-1)),((a,b-2,c),b*(b-1)),((a,b,c-2),c*(c-1))]:
        if factor and exp in rows:C[rows.index(exp),col]+=factor
    N=null_space(C) if len(rows) else np.eye(len(exps))
    return exps,N

def eval_poly(q,exps,coef,R,deriv=False):
    q=np.asarray(q,float); s=q/R; val=np.zeros(q.shape[0]);g=np.zeros((q.shape[0],3))
    for c,(a,b,d) in zip(coef,exps):
      val+=c*s[:,0]**a*s[:,1]**b*s[:,2]**d
      for ax,powr in enumerate((a,b,d)):
       if powr:
        powers=[a,b,d];powers[ax]-=1
        g[:,ax]+=c*powr/R*np.prod(s**np.asarray(powers),axis=1)
    return val,g

def fit_harmonic(basis,u,order,R):
    exps,N=harmonic_even_basis(order);coords=basis.doflocs.T
    mask=np.linalg.norm(coords,axis=1)<R
    q=coords[mask]; vals=u[mask]; s=q/R
    M=np.column_stack([np.prod(s**np.array(e),axis=1) for e in exps])
    A=M@N
    coef_reduced,resid,rank,sv=np.linalg.lstsq(A,vals,rcond=None);coef=N@coef_reduced
    cond=float(sv[0]/sv[-1]) if len(sv) and sv[-1]>0 else float('inf')
    return exps,coef,{'fit_dofs':int(mask.sum()),'basis_terms':len(exps),'harmonic_basis_rank':N.shape[1],
      'matrix_rank':int(rank),'condition_number':cond,'potential_rms_residual_V_per_basis_V':float(np.sqrt(np.mean((M@coef-vals)**2)))}

def fit_free(basis,u,order,R):
    coords=basis.doflocs.T;mask=np.linalg.norm(coords,axis=1)<R;q=coords[mask];vals=u[mask];s=q/R
    exps=[(i,j,k) for i in range(order+1) for j in range(order+1-i) for k in range(order+1-i-j)
          if i+j+k<=order and all(a%2==0 for a in (i,j,k))]
    M=np.column_stack([np.prod(s**np.array(e),axis=1) for e in exps]);co,resid,rank,sv=np.linalg.lstsq(M,vals,rcond=None)
    hdiag=[]
    for axis in range(3):
      ee=[0,0,0];ee[axis]=2
      hdiag.append(float(2*co[exps.index(tuple(ee))]/R**2) if tuple(ee) in exps else 0.)
    cond=float(sv[0]/sv[-1]) if len(sv) and sv[-1]>0 else float('inf')
    return {'fit_dofs':int(mask.sum()),'order':order,'fit_radius_r0':R,'condition_number':cond,'matrix_rank':int(rank),
      'hessian_diagonal_per_basis_V_per_r0_squared':hdiag,'trace_per_basis_V_per_r0_squared':float(sum(hdiag)),
      'potential_rms_residual_V_per_basis_V':float(np.sqrt(np.mean((M@co-vals)**2)))}

def main():
  pts=sample_points();octs={tag:{kind:load_oct(ROOT/'results'/f'm3a3_octant_fields_{tag}.npz',key) for kind,key in [('dc','phi_dc'),('rf','phi_rf')]} for tag in ('plateau030','plateau018')}
  mesh_change={}
  for kind in ('dc','rf'):
    b0,u0=octs['plateau030'][kind];b1,u1=octs['plateau018'][kind]
    _,g0=eval_p2(b0,u0,abs(pts));_,g1=eval_p2(b1,u1,abs(pts));g0*=np.where(pts.T<0,-1.,1.);g1*=np.where(pts.T<0,-1.,1.)
    volt=VDC if kind=='dc' else VRF_DIFF;scale=E_CHARGE*volt/R0_M
    F0=-scale*g0.T;F1=-scale*g1.T;errs=[]
    for i,p in enumerate(pts):
      den=np.linalg.norm(F1[i]);errs.append({'point_r0':p.tolist(),'F_coarse_N':F0[i].tolist(),'F_fine_N':F1[i].tolist(),
       'delta_N':(F0[i]-F1[i]).tolist(),'vector_change_pct':float(100*np.linalg.norm(F0[i]-F1[i])/den) if den>1e-45 else None,
       'magnitude_change_pct':float(100*abs(np.linalg.norm(F0[i])-den)/den) if den>1e-45 else None,
       'absolute_error_eV_per_r0':float(np.linalg.norm(F0[i]-F1[i])/(E_CHARGE/R0_M))})
    mesh_change[kind]={'max_vector_pct':max(v['vector_change_pct'] or 0 for v in errs),'max_magnitude_pct':max(v['magnitude_change_pct'] or 0 for v in errs),
      'max_absolute_error_N':max(np.linalg.norm(v['delta_N']) for v in errs),'pointwise':errs}
  recovery={}; hessian_mesh_levels={}
  # Held-out points confined to positive octant and radius <= 0.24 r0.
  held=np.array([[.037,.081,.113],[.103,.047,.157],[.176,.091,.053],[.064,.174,.121],[.211,.039,.086],[.121,.132,.169]])
  for kind in ('dc','rf'):
    b,u=octs['plateau018'][kind];_,gref=eval_p2(b,u,held);factor=VDC if kind=='dc' else VRF_DIFF
    rows=[]
    for order in (2,4,6):
      for radius in (.35,.50):
        ex,co,meta=fit_harmonic(b,u,order,radius);_,gp=eval_poly(held,ex,co,radius)
        # Compare force from potential gradient, not independently filtered components.
        Fre=-E_CHARGE*factor/R0_M*gref.T;Ffit=-E_CHARGE*factor/R0_M*gp
        denom=np.linalg.norm(Fre,axis=1)
        meta.update({'order':order,'fit_radius_r0':radius,'max_relative_vector_error_pct':float(np.max(100*np.linalg.norm(Ffit-Fre,axis=1)/denom)),
          'rms_relative_vector_error_pct':float(np.sqrt(np.mean((100*np.linalg.norm(Ffit-Fre,axis=1)/denom)**2))),
          'max_relative_magnitude_error_pct':float(np.max(100*np.abs(np.linalg.norm(Ffit,axis=1)-denom)/denom)),
          'heldout_points_r0':held.tolist()})
        rows.append(meta)
    unconst=[fit_free(b,u,o,R) for o in (2,4,6) for R in (.35,.50)]
    recovery[kind]={'harmonic_constrained_heldout':rows,'unconstrained_even_hessian':unconst}
    hessian_mesh_levels[kind]={}
    for tag in ('plateau030','plateau018'):
      bm,um=octs[tag][kind]
      hessian_mesh_levels[kind][tag]=[fit_free(bm,um,2,R) for R in (.35,.50)]
  out={'fixed_voltage_V':{'dc_each_ring':VDC,'rf_differential_peak':VRF_DIFF},'mesh_to_mesh':mesh_change,
    'harmonic_even_potential_recovery':recovery,'unconstrained_quadratic_hessian_mesh_levels':hessian_mesh_levels,
    'note':'Potential fits satisfy Laplace by construction; held-out agreement tests representation fidelity but is not an independent proof of FEM PDE accuracy.'}
  (ROOT/'results/m3a3_octant_recovery.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
  print(json.dumps({'mesh_to_mesh_summary':{k:{a:b for a,b in v.items() if a!='pointwise'} for k,v in mesh_change.items()},'recovery':recovery},indent=2))
if __name__=='__main__':main()
