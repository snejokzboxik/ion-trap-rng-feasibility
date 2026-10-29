"""Audit geometry and excitation parity before octant electrostatic FEM."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.readout_fem import truncated_electrode_polygon

polys={b:truncated_electrode_polygon(b,1.,.10,1.5,80).T for b in ('right','left','top','bottom')}
def hausdorff(a,b):
    ta,tb=cKDTree(a),cKDTree(b)
    return float(max(ta.query(b)[0].max(),tb.query(a)[0].max()))
tests={
 'reflect_x_right_to_left':hausdorff(polys['right']*[-1,1],polys['left']),
 'reflect_y_top_to_bottom':hausdorff(polys['top']*[1,-1],polys['bottom']),
 'reflect_y_right_to_right':hausdorff(polys['right']*[1,-1],polys['right']),
 'reflect_x_top_to_top':hausdorff(polys['top']*[-1,1],polys['top'])}
hyperbola={}
from iontrap_rng.readout_fem import hyperbola_branch
for b in polys:
    p=hyperbola_branch(np.linspace(-np.arcsinh(1.5),np.arcsinh(1.5),161),b,1.)
    sign=1. if b in ('right','left') else -1.
    hyperbola[b]=float(np.max(np.abs(p[0]**2-p[1]**2-sign)))

# RF physical basis is +1/2 on right and left; -1/2 on top and bottom.
rf={'right':.5,'left':.5,'top':-.5,'bottom':-.5}
rf_reflections={
 'x':{'right':'left','left':'right','top':'top','bottom':'bottom'},
 'y':{'right':'right','left':'left','top':'bottom','bottom':'top'},
 'z':{b:b for b in rf}}
rf_invariant={axis:all(rf[key]==rf[image] for key,image in perm.items()) for axis,perm in rf_reflections.items()}
dc={'dc_plus':1.,'dc_minus':1.}
dc_reflections={'x':{'dc_plus':'dc_plus','dc_minus':'dc_minus'},
 'y':{'dc_plus':'dc_plus','dc_minus':'dc_minus'},
 'z':{'dc_plus':'dc_minus','dc_minus':'dc_plus'}}
dc_invariant={axis:all(dc[key]==dc[image] for key,image in perm.items()) for axis,perm in dc_reflections.items()}
out={'purpose':'preflight geometry and excitation symmetry; no FEM',
 'geometry':{'rf_profile_hausdorff_residual_r0':tests,'central_hyperbola_max_residual_r0_squared':hyperbola,
   'rf_rail_z_extent_r0':[-5.,5.],'dc_ring_centers_z_r0':[-5.5,5.5],
   'dc_ring_thickness_r0':.2,'dc_ring_inner_outer_r0':[.5,1.5],
   'grounded_enclosure_half_extents_r0':[6.,6.,7.]},
 'excitation_patterns':{'rf_terminal_values_V_per_Vdiff':rf,'rf_invariant_under_reflection':rf_invariant,
   'dc_terminal_values_V_per_Vring':dc,'dc_invariant_under_reflection':dc_invariant,
   'individual_terminal_basis_fields_are_not_all_invariant':True},
 'octant_boundary_conditions':{
   'open_vacuum_parts_of_x0_y0_z0':'homogeneous Neumann by even scalar-potential symmetry',
   'electrode_cross_sections_on_symmetry_planes':'Dirichlet at that conductor terminal voltage, not Neumann',
   'grounded_outer_faces':'Dirichlet zero',
   'physical_RF_and_DC_surfaces':'Dirichlet at paired terminal voltage',
   'individual_terminal_weighting_fields':'do not reuse octant BC unless that specific excitation is proven even'},
 'interpretation':'Expected exact geometric/electrical symmetry permits octant reduction for the paired RF and symmetric DC patterns. Numerical symmetry enforced on octant planes is not an accuracy check.'}
(ROOT/'results/m3a3_symmetry_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
