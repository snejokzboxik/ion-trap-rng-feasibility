"""Fixed-voltage DC-force and local-polynomial audit of saved M2b.2 meshes."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
from iontrap_rng.hessian import fit_local_hessian

r0=1e-6;Q=elementary_charge;vdc=.08697882845042744;vrf=10.
cases={
 'h015':('m2b2_dc_surface_refined_p2_fields_h0p15.npz','m2b2_p2_h0.75_hc0.18_dc0.15.json'),
 'h010':('m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz','m2b2_p2_effective_h0.75_hc0.18_dc0.1.json'),
 'h0075':('m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz','m2b2_p2_effective_h0.75_hc0.18_dc0.075.json')}
points={}
for axis in range(3):
  for s in (-.3,-.2,-.1,-.05,.05,.1,.2,.3):
    x=np.zeros(3);x[axis]=s;points[f'axis_{"xyz"[axis]}_{s:+.2f}']=x
points.update({'off_xy':[.1,.1,0.],'off_xz':[.1,0.,.1],'off_yz':[0.,.1,.1],
               'off_xyz':[.08,-.06,.10]})
fields={};fits={};meta={};hzz_per_v={}
for key,(npz,jsonname) in cases.items():
  d=json.loads((ROOT/'results'/jsonname).read_text(encoding='utf-8'));meta[key]=d['mesh']
  hzz_per_v[key]=float(d['fields']['dc'][1]['Hessian'][2][2])
  ff=ROOT/'results'/npz
  fields[key]={'dc':P2ScalarField(ff,'dc_potential'),'rf':P2ScalarField(ff,'rf_potential')}
  raw=np.load(ff);coords=raw['p2_dof_coordinates'];u=raw['dc_potential']
  pts=coords[np.linalg.norm(coords,axis=1)<=.5]
  ids=np.linalg.norm(coords,axis=1)<=.5
  fits[key]={}
  for hcon in (False,True):
    f=fit_local_hessian(coords,u,radius=.5,order=2,harmonic_quadratic=hcon)
    fits[key]['harmonic' if hcon else 'unconstrained']={'fit':f,'samples':f.samples,'rank':f.rank,
      'condition':f.condition_number,'residual':f.relative_residual,
      'trace':float(np.trace(f.hessian))}

rows=[]
for name,x in points.items():
  x=np.asarray(x,float); row={'point':name,'xi':x.tolist(),'radius_r0':float(np.linalg.norm(x)),'meshes':{}}
  for key in cases:
    dcg=fields[key]['dc'].gradient(x[:,None])[:,0]
    rfg=fields[key]['rf'].gradient(x[:,None])[:,0]
    md={}
    for mode in ('unconstrained','harmonic'):
      f=fits[key][mode]['fit'];md[mode+'_gradient']= (f.gradient+f.hessian@x).tolist()
    md.update({'dc_gradient_direct':dcg.tolist(),'rf_gradient_direct':rfg.tolist(),
      'dc_force_N':(-Q*vdc/r0*dcg).tolist(),'dc_force_magnitude_N':float(Q*vdc/r0*np.linalg.norm(dcg)),
      'rf_force_peak_N':(-Q*vrf/r0*rfg).tolist()})
    row['meshes'][key]=md
  fine=row['meshes']['h0075'];coarse=row['meshes']['h010']
  dd=np.asarray(coarse['dc_gradient_direct'])-np.asarray(fine['dc_gradient_direct'])
  Fd=Q*vdc/r0*np.linalg.norm(dd);Ffine= fine['dc_force_magnitude_N']
  pred=np.asarray(fine['unconstrained_gradient'])
  dfine=np.asarray(fine['dc_gradient_direct'])-pred
  row['comparison_h010_h0075']={'dc_force_difference_N':float(Fd),
    'dc_force_difference_over_fine_dc_force':float(Fd/max(Ffine,1e-30)),
    'dc_gradient_difference_per_V':dd.tolist(),
    'finest_direct_minus_unconstrained_fit_gradient_per_V':dfine.tolist(),
    'finest_direct_minus_fit_force_N':float(Q*vdc/r0*np.linalg.norm(dfine)),
    'force_error_over_fine_force':float(Q*vdc/r0*np.linalg.norm(dfine)/max(Ffine,1e-30)),
    'force_error_over_total_force_rf_phase0':float(Q*vdc/r0*np.linalg.norm(dd)/max(np.linalg.norm(np.asarray(fine['rf_force_peak_N'])+np.asarray(fine['dc_force_N'])),1e-30))}
  rows.append(row)

fd_rows=[]
fine_field=fields['h0075']['dc']
for name in ('axis_x_+0.10','axis_y_+0.10','axis_z_+0.10','axis_z_-0.10','off_xy','off_xz','off_yz','off_xyz'):
  x=np.asarray(points[name],float);direct=fine_field.gradient(x[:,None])[:,0]
  entry={'point':name,'direct_gradient':direct.tolist(),'finite_difference':{}}
  for eps in (.005,.01,.02):
    g=[]
    for axis in range(3):
      xp=x.copy();xm=x.copy();xp[axis]+=eps;xm[axis]-=eps
      g.append(float((fine_field.potential(xp[:,None])[0]-fine_field.potential(xm[:,None])[0])/(2*eps)))
    g=np.asarray(g)
    entry['finite_difference'][str(eps)]={'gradient':g.tolist(),
      'difference_norm_per_V':float(np.linalg.norm(g-direct)),
      'relative_to_direct_norm':float(np.linalg.norm(g-direct)/max(np.linalg.norm(direct),1e-30))}
  fd_rows.append(entry)

fit_summary={k:{mode:{q:v for q,v in vrow.items() if q!='fit'} for mode,vrow in modes.items()} for k,modes in fits.items()}
mass=39.96259098*1.66053906660e-27
freqs={k:float(np.sqrt(Q*vdc*h/(mass*r0*r0))/(2*np.pi)) for k,h in hzz_per_v.items()}
out={'purpose':'fixed-Vdc audit before any extended integration; no trajectories',
 'field_normalization':{'rf':'1 V differential peak pattern','dc':'1 V on each DC ring','coordinates':'r/r0'},
 'fixed_physical_dc_voltage_per_ring_V':vdc,'rf_differential_peak_V':vrf,'r0_m':r0,
 'meshes':meta,'dc_Hzz_per_V_r0_squared':hzz_per_v,'axial_frequency_fixed_same_VDC_Hz':freqs,
 'local_value_fit_R_r0':.5,'local_fits':fit_summary,
 'positions':rows,'finite_difference_gradient_crosscheck_h0075':fd_rows,
 'interpretation_rules':{'relative_force_errors':'reported relative to local fine-mesh DC force only; near-zero force values must be read with absolute force difference',
    'harmonic_fit':'trace-constrained comparator; not independent proof of Laplace satisfaction',
    'unconstrained_fit':'local potential-value polynomial fit from P2 DOFs; compared against elementwise direct P2 gradient'},
 'limitations':['meshes refine DC surfaces only; central, RF and axial-end mesh controls remain fixed',
    'P2 element gradients are discontinuous across tetrahedron faces',
    'local fit and FEM solution use the same discrete potential; comparison tests derivative recovery, not independent PDE error']}
(ROOT/'results'/'m3a1_dc_force_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'fixed_vdc':vdc,'fits':fit_summary,'points':[
 {'point':r['point'],'r':r['radius_r0'],'force_h075_N':r['meshes']['h0075']['dc_force_N'],
  'mesh_delta_over_force':r['comparison_h010_h0075']['dc_force_difference_over_fine_dc_force'],
  'fit_delta_over_force':r['comparison_h010_h0075']['force_error_over_fine_force']}
 for r in rows]},indent=2))
