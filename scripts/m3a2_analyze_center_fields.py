"""Compare saved central-refinement P2 fields at identical physical points."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as Q
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
from iontrap_rng.hessian import fit_local_hessian
sys.path.insert(0,str(ROOT/'scripts'))
from run_m2b_3d_preflight import boundary_values
r0=1e-6;vdc=.08697882845042744;vrf=10.;mass=39.96259098*1.66053906660e-27
files={
 'hc018':'results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz',
 'hc009':'results/m3a2_center_fields_hc0p09.npz',
 'hc006':'results/m3a2_center_fields_hc0p06.npz'}
fields={k:{'dc':P2ScalarField(ROOT/f,'dc_potential'),'rf':P2ScalarField(ROOT/f,'rf_potential')} for k,f in files.items()}
fits={};mesh_meta={}
for k,fpath in files.items():
    d=np.load(ROOT/fpath); xyz=d['p2_dof_coordinates']; fits[k]={}
    for radius in (.35,.50):
        fits[k][str(radius)]={}
        for key in ('dc','rf'):
            result=fit_local_hessian(xyz,d[f'{key}_potential'],radius=radius,order=2,harmonic_quadratic=False)
            fits[k][str(radius)][key]={'hessian':result.hessian,'gradient':result.gradient,
              'condition':result.condition_number,'samples':result.samples,'rank':result.rank,'residual':result.relative_residual}
    if k=='hc018': mesh_meta[k]=json.loads((ROOT/'results/m2b2_p2_effective_h0.75_hc0.18_dc0.075.json').read_text())['mesh']
    else:
        htag=str(int(k[2:])/100).replace('.','p')
        meta=json.loads((ROOT/'results'/f'm3a2_center_mesh_hc{htag}.json').read_text())
        meta['P2_dofs']=int(len(d['dc_potential']))
        mesh_meta[k]=meta

points={}
for ax in range(3):
  for s in (-.20,-.10,-.05,.05,.10,.20):
    x=np.zeros(3);x[ax]=s;points[f'axis_{"xyz"[ax]}_{s:+.2f}']=x
points.update({'off_xy':[.1,.1,0.],'off_xz':[.1,0.,.1],'off_yz':[0.,.1,.1],
               'off_xyz':[.08,-.06,.10]})
def field_force(g,voltage): return -Q*voltage/r0*np.asarray(g)
point_rows=[]
for name,x in points.items():
 row={'point':name,'xi':np.asarray(x).tolist(),'by_mesh':{}}
 for m in files:
  dcg=fields[m]['dc'].gradient(np.asarray(x)[:,None])[:,0]
  rfg=fields[m]['rf'].gradient(np.asarray(x)[:,None])[:,0]
  row['by_mesh'][m]={'dc_gradient_V_per_xi_per_V':dcg.tolist(),'dc_force_N':field_force(dcg,vdc).tolist(),
    'dc_force_magnitude_N':float(np.linalg.norm(field_force(dcg,vdc))),
    'rf_gradient_V_per_xi_per_Vdiff':rfg.tolist(),'rf_peak_force_N':field_force(rfg,vrf).tolist(),
    'polynomial_gradient_predictions':{}}
  for rad in (.35,.50):
   f=fits[m][str(rad)]['dc'];g=np.asarray(f['gradient'])+np.asarray(f['hessian'])@x
   row['by_mesh'][m]['polynomial_gradient_predictions'][str(rad)]={'gradient':g.tolist(),
     'force_N':field_force(g,vdc).tolist(),'absolute_force_difference_N':float(np.linalg.norm(field_force(dcg-g,vdc)))}
 for c,f in (('hc018','hc009'),('hc009','hc006'),('hc018','hc006')):
  fc=np.asarray(row['by_mesh'][c]['dc_force_N']);ff=np.asarray(row['by_mesh'][f]['dc_force_N']);den=np.linalg.norm(ff)
  row.setdefault('mesh_comparisons',{})[c+'_vs_'+f]={'vector_difference_N':float(np.linalg.norm(fc-ff)),
   'relative_vector_force_difference':float(np.linalg.norm(fc-ff)/den) if den>1e-25 else None,
   'force_magnitude_difference_N':float(abs(np.linalg.norm(fc)-den)),
   'relative_force_magnitude_difference':float(abs(np.linalg.norm(fc)-den)/den) if den>1e-25 else None,
   'component_difference_N':(fc-ff).tolist()}
 point_rows.append(row)

# Mirror parity: DC potential should be even, and its gradient transforms as a vector under reflection.
symmetry={}
for mesh in files:
 rows=[]
 for ax in range(3):
  S=np.eye(3);S[ax,ax]=-1
  for name,x0 in points.items():
   x=np.asarray(x0,float)
   if abs(x[ax])<1e-15: continue
   xp=x.copy();xm=x.copy();xp[ax]=abs(x[ax]);xm[ax]=-abs(x[ax])
   f=fields[mesh]['dc'];vp=float(f.potential(xp[:,None])[0]);vm=float(f.potential(xm[:,None])[0])
   gp=f.gradient(xp[:,None])[:,0];gm=f.gradient(xm[:,None])[:,0]
   delta=gm-S@gp;den=max(np.linalg.norm(gm),np.linalg.norm(gp),1e-30)
   rows.append({'axis':'xyz'[ax],'point':name,'positive_xi':xp.tolist(),'negative_xi':xm.tolist(),
      'potential_even_abs_residual_V_per_Vring':abs(vp-vm),
      'potential_even_relative_residual':abs(vp-vm)/max(abs(vp),abs(vm),1e-12),
      'gradient_parity_abs_residual_V_per_xi':float(np.linalg.norm(delta)),
      'gradient_parity_relative_residual':float(np.linalg.norm(delta)/den),
      'gradient_plus':gp.tolist(),'gradient_minus':gm.tolist()})
 symmetry[mesh]={'pairs':rows,
   'max_potential_even_abs_residual':max(x['potential_even_abs_residual_V_per_Vring'] for x in rows),
   'max_gradient_parity_abs_residual':max(x['gradient_parity_abs_residual_V_per_xi'] for x in rows),
   'rms_gradient_parity_relative_residual':float(np.sqrt(np.mean([x['gradient_parity_relative_residual']**2 for x in rows])))}

# Independent local Hessian estimate from sampled direct gradients (not a
# harmonic-constrained fit). This tests derivative recovery of the discrete field.
sob=qmc.Sobol(d=3,scramble=False).random_base2(9)
rad=sob[:,0]**(1./3.)*.30
az=2*np.pi*sob[:,2]; zdir=1.-2.*sob[:,1]; rxy=np.sqrt(np.maximum(0.,1.-zdir*zdir))
samples=np.column_stack((rad*rxy*np.cos(az),rad*rxy*np.sin(az),rad*zdir))
gradient_recovery={}
for mesh in files:
    g=fields[mesh]['dc'].gradient(samples.T)
    M=np.column_stack((np.ones(len(samples)),samples/.30))
    coef,_,rank,svals=np.linalg.lstsq(M,g.T,rcond=None)
    H=coef[1:,:].T/.30
    Hsym=.5*(H+H.T)
    potential_H=np.asarray(fits[mesh]['0.35']['dc']['hessian'])
    gradient_recovery[mesh]={'samples':int(len(samples)),'radius_r0':.30,'rank':int(rank),
      'condition_number':float(svals[0]/svals[-1]),'gradient_fit_relative_residual':float(np.linalg.norm(M@coef-g.T)/np.linalg.norm(g.T)),
      'gradient_at_origin_V_per_xi_per_V':coef[0].tolist(),'raw_hessian_V_per_r0_squared':H.tolist(),
      'symmetrized_hessian_V_per_r0_squared':Hsym.tolist(),'trace_raw':float(np.trace(H)),
      'antisymmetric_fraction':float(np.linalg.norm(H-H.T)/max(np.linalg.norm(H),1e-30)),
      'difference_from_nodal_potential_hessian_frobenius_relative':float(np.linalg.norm(Hsym-potential_H)/np.linalg.norm(potential_H))}

boundary_audit={}
for mesh in files:
    basis=fields[mesh]['dc'].basis
    boundary=np.unique(basis.get_dofs(facets=basis.mesh.boundary_facets()).all())
    mask=np.zeros(basis.N,dtype=bool);mask[boundary]=True
    raw=boundary_values(basis.doflocs.T,p2_surface_tolerance=.02)
    tagged={k:np.asarray(v,int)[mask[np.asarray(v,int)]] for k,v in raw.items()}
    fixed=np.unique(np.concatenate(list(tagged.values())))
    overlaps=[];keys=list(tagged)
    for i,a in enumerate(keys):
      for b in keys[i+1:]:
        n=len(np.intersect1d(tagged[a],tagged[b]))
        if n: overlaps.append({'terminals':[a,b],'count':int(n)})
    boundary_audit[mesh]={'boundary_dofs':int(len(boundary)),'missing_boundary_dofs':int(len(np.setdiff1d(boundary,fixed))),
      'tagged_dof_counts':{k:int(len(v)) for k,v in tagged.items()},'terminal_overlap':overlaps}

# Element-neighborhood sensitivity at the highlighted +/-z points; no force is symmetrized.
boundary_sensitivity={}
point_cell_clearance={}
for mesh in files:
 f=fields[mesh]['dc'];items={};clearance={}
 for x0 in ([0,0,.1],[0,0,-.1],[.1,0,0],[-.1,0,0]):
  x=np.array(x0,float);g0=f.gradient(x[:,None])[:,0];eps_rows=[]
  cell=int(np.asarray(f.find_element(*x[:,None])).ravel()[0])
  verts=f.mesh.p[:,f.mesh.t[:,cell]].T
  lam=np.linalg.solve((verts[1:]-verts[0]).T,x-verts[0]); bary=np.r_[1.-lam.sum(),lam]
  face_heights=[]
  for omit in range(4):
   face=np.delete(verts,omit,axis=0)
   area=.5*np.linalg.norm(np.cross(face[1]-face[0],face[2]-face[0]))
   volume=abs(np.linalg.det((verts[1:]-verts[0]).T))/6.
   face_heights.append(3.*volume/area)
  clearance[str(x0)]={'selected_tetrahedron':cell,'barycentric_coordinates':bary.tolist(),
    'minimum_distance_to_selected_tetrahedron_face_r0':float(np.min(bary*np.asarray(face_heights)))}
  for eps in (1e-8,1e-6,1e-4):
   variations=[]
   for ax in range(3):
    for sign in (-1,1):
     xp=x.copy();xp[ax]+=sign*eps;gp=f.gradient(xp[:,None])[:,0]
     variations.append(np.linalg.norm(gp-g0))
   eps_rows.append({'offset_r0':eps,'max_gradient_change_V_per_xi_per_V':float(max(variations)),
      'relative_to_center_gradient':float(max(variations)/max(np.linalg.norm(g0),1e-30))})
  items[str(x0)]=eps_rows
 boundary_sensitivity[mesh]=items
 point_cell_clearance[mesh]=clearance

summary={}
for m in files:
 summary[m]={}
 for rad,by in fits[m].items():
  summary[m][rad]={}
  for key,v in by.items():
   H=np.asarray(v['hessian']); freq=float(np.sqrt(Q*vdc*H[2,2]/(mass*r0*r0))/(2*np.pi)) if key=='dc' and H[2,2]>0 else None
   summary[m][rad][key]={'hessian_V_per_r0_squared':H.tolist(),'trace':float(np.trace(H)),
    'fit_residual':float(v['residual']),'condition':float(v['condition']),'samples':int(v['samples']),'rank':int(v['rank']),
    'axial_frequency_from_Hzz_fixed_VDC_Hz':freq}

out={'purpose':'central mesh refinement force and symmetry validation; fixed voltages; no trajectories',
 'force_convergence_gate_passed':False,
 'conditional_trajectory_runs':[],
 'decision':'NO-GO for 0.1 or 1 axial-period integration: successive central refinements retain multi-percent vector-force error and sampled gradient-parity residual',
 'voltage_and_units':{'VDC_each_ring_V':vdc,'VRF_differential_peak_V':vrf,'r0_m':r0,'Ca40_mass_kg':mass,
   'field_gradient_coordinates':'dimensionless xi=r/r0; force=-Q*V/r0*grad_xi(phi_basis)'},
 'mesh_files':files,'mesh_metadata':mesh_meta,'hessian_fits':summary,'force_point_comparisons':point_rows,
 'dc_mirror_symmetry':symmetry,'gradient_based_hessian_recovery':gradient_recovery,
 'element_neighborhood_gradient_sensitivity':boundary_sensitivity,'point_cell_clearance':point_cell_clearance,
 'boundary_tag_audit':boundary_audit,
 'limits':['all compared solutions share discrete P2 solver and polygonized geometry','elementwise P2 gradient is discontinuous','potential and field mirror symmetry are sampled, not a proof of continuum PDE convergence','no symmetry is imposed in postprocessing']}
(ROOT/'results/m3a2_center_field_comparison.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'hessian_fits':summary,'symmetry_summary':{k:{q:v for q,v in d.items() if q!='pairs'} for k,d in symmetry.items()},
 'pair_comparison_count':len(point_rows)},indent=2))
