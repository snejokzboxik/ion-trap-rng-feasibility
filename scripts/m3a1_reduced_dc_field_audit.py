"""Test local polynomial DC-field reductions against saved P2 FEM fields."""
import json,sys
from itertools import product
from pathlib import Path
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
from iontrap_rng.hessian import fit_local_hessian
CASES={'h015':'m2b2_dc_surface_refined_p2_fields_h0p15.npz',
       'h010':'m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',
       'h0075':'m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz'}
def powers(order):return [e for e in product(range(order+1),repeat=3) if sum(e)<=order]
def design(q,exps):return np.column_stack([np.prod(q**np.asarray(e)[None,:],axis=1) for e in exps])
def grad_model(c,exps,x,R):
    q=np.asarray(x,float)/R;g=np.zeros(3)
    for k,e in enumerate(exps):
      for j in range(3):
        if e[j]:
          p=list(e);p[j]-=1;g[j]+=c[k]*e[j]*np.prod(q**np.asarray(p))/R
    return g
u=qmc.Sobol(3,scramble=False).random_base2(8)
z=2*u[:,0]-1;ph=2*np.pi*u[:,1];rad=.3*np.cbrt(u[:,2]);ss=np.sqrt(1-z*z)
evalpts=(rad[:,None]*np.column_stack((ss*np.cos(ph),ss*np.sin(ph),z))).T
coords={}
for axis in range(3):
  for s in (-.2,-.1,-.05,.05,.1,.2):
    x=np.zeros(3);x[axis]=s;coords[f'{"xyz"[axis]}{s:+.2f}']=x
coords.update({'off1':np.array([.1,.1,0.]),'off2':np.array([.1,0,.1]),'off3':np.array([0,.1,.1]),'off4':np.array([.08,-.06,.1])})
out={'fit_radius_r0':.5,'validation_ball_r0':.3,'fixed_dc_V_per_ring':.08697882845042744,'models':{},
     'warning':'fit is to saved P2 nodal potential; same discrete PDE solution, not independent solver truth'}
for mesh,npz in CASES.items():
  d=np.load(ROOT/'results'/npz);field=P2ScalarField(ROOT/'results'/npz,'dc_potential')
  p=d['p2_dof_coordinates'];v=d['dc_potential'];use=np.linalg.norm(p,axis=1)<=.5;R=.5;q=p[use]/R
  rows={}
  direct_eval=field.gradient(evalpts)
  for order in (2,3,4):
    exps=powers(order);X=design(q,exps);coef,_,rank,s=np.linalg.lstsq(X,v[use],rcond=None)
    res=float(np.linalg.norm(X@coef-v[use])/np.linalg.norm(v[use]));cond=float(s[0]/s[-1])
    gmodel=np.column_stack([grad_model(coef,exps,x,R) for x in evalpts.T])
    dd=gmodel-direct_eval
    rms_delta=float(np.sqrt(np.mean(np.sum(dd**2,axis=0))))
    rms_direct=float(np.sqrt(np.mean(np.sum(direct_eval**2,axis=0))))
    pointrows=[]
    for name,x in coords.items():
      gd=field.gradient(x[:,None])[:,0];gp=grad_model(coef,exps,x,R)
      pointrows.append({'point':name,'xi':x.tolist(),'direct_gradient':gd.tolist(),'polynomial_gradient':gp.tolist(),
        'gradient_difference_norm_per_V':float(np.linalg.norm(gd-gp)),
        'relative_to_direct_norm':float(np.linalg.norm(gd-gp)/max(np.linalg.norm(gd),1e-30))})
    rows[f'order{order}']={'terms':len(exps),'samples':int(use.sum()),'rank':int(rank),'condition':cond,
      'relative_potential_fit_residual_nodes':res,'validation_ball_gradient_difference_rms_per_V':rms_delta,
      'validation_ball_direct_gradient_rms_per_V':rms_direct,'relative_gradient_difference_rms':rms_delta/max(rms_direct,1e-30),
      'point_comparisons':pointrows}
  f2=fit_local_hessian(p,v,radius=.5,order=2,harmonic_quadratic=False)
  fh=fit_local_hessian(p,v,radius=.5,order=2,harmonic_quadratic=True)
  rows['quadratic_fit_hessians']={'unconstrained':{'H':f2.hessian.tolist(),'g':f2.gradient.tolist(),'trace':float(np.trace(f2.hessian))},
    'harmonic_constrained':{'H':fh.hessian.tolist(),'g':fh.gradient.tolist(),'trace':float(np.trace(fh.hessian))}}
  out['models'][mesh]=rows
(ROOT/'results'/'m3a1_reduced_dc_field_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({m:{k:{z:v for z,v in d.items() if z not in ('point_comparisons','quadratic_fit_hessians')} for k,d in r.items() if k.startswith('order')} for m,r in out['models'].items()},indent=2))
