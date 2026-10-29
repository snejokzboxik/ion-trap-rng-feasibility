"""Generate one targeted DC-electrode-surface-refined mesh, no field solve."""
import json
from pathlib import Path
import sys
import argparse
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_m2b_3d_preflight import make_mesh
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--dc-surface-h',type=float,default=.15);args=parser.parse_args()
surface_h=args.dc_surface_h
p,t,meta=make_mesh(.75,.18,dc_surface_h=surface_h,rf_surface_h=.30,axial_end_h=.30,return_metadata=True)
edges=np.concatenate((t[:,[0,1]],t[:,[0,2]],t[:,[0,3]],t[:,[1,2]],t[:,[1,3]],t[:,[2,3]]))
edges=np.unique(np.sort(edges,axis=1),axis=0); edgecount=len(edges)
mid=p[edges].mean(axis=1); lengths=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
rho=np.linalg.norm(mid[:,:2],axis=1); az=np.abs(mid[:,2])
near_wall=(np.abs(az-5.5)<.15)&(np.minimum(np.abs(rho-.5),np.abs(rho-1.5))<.12)
near_face=(np.minimum(np.abs(az-5.4),np.abs(az-5.6))<.025)&(rho>.45)&(rho<1.55)
local=near_wall|near_face
# Boundary facets are faces appearing once in the tetrahedral complex.
faces=np.concatenate((t[:,[1,2,3]],t[:,[0,2,3]],t[:,[0,1,3]],t[:,[0,1,2]]))
sf,counts=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True); bf=sf[counts==1]
bc=p[bf].mean(axis=1);br=np.linalg.norm(bc[:,:2],axis=1);bz=np.abs(bc[:,2])
face_tag=(np.minimum(np.abs(bz-5.4),np.abs(bz-5.6))<.02)&(br>=.49)&(br<=1.51)
wall_tag=(np.abs(bz-5.5)<.101)&(np.minimum(np.abs(br-.5),np.abs(br-1.5))<.03)
ring_faces=bf[face_tag|wall_tag]
ring_edges=np.unique(np.sort(np.concatenate((ring_faces[:,[0,1]],ring_faces[:,[0,2]],ring_faces[:,[1,2]])),axis=1),axis=0)
ring_l=np.linalg.norm(p[ring_edges[:,0]]-p[ring_edges[:,1]],axis=1)
tri_area=.5*np.linalg.norm(np.cross(p[ring_faces[:,1]]-p[ring_faces[:,0]],p[ring_faces[:,2]]-p[ring_faces[:,0]]),axis=1)
def stats(v):return {'count':int(len(v)),'min':float(np.min(v)),'median':float(np.median(v)),'p90':float(np.percentile(v,90)),'p99':float(np.percentile(v,99)),'max':float(np.max(v))}
out={'mesh_parameters':{'hmax_r0':.75,'hcenter_r0':.18,'dc_surface_h_r0':surface_h,'rf_surface_h_r0':.30,'axial_end_h_r0':.30},
     'geometry_metadata':meta,'vertices':len(p),'tetrahedra':len(t),'unique_edges':int(edgecount),
     'estimated_P2_dofs':int(len(p)+edgecount),
     'local_edge_lengths_r0':{'DC_surface_adjacent':stats(lengths[local]),'global':stats(lengths)},
     'DC_surface_triangle_edge_lengths_r0':stats(ring_l),
     'DC_surface_triangle_areas_r0_squared':stats(tri_area),'DC_boundary_triangle_count':int(len(ring_faces)),
     'purpose':'bounded mesh-size probe; no matrix assembly or electrostatic solve',
     'limitations':['still polygonal/linear tetrahedral geometry','ideal conductors in vacuum','local DC/RF surface and axial-end fields; far field remains coarse']}
tag=f'{surface_h:g}'.replace('.','p')
(ROOT/'results'/f'm2b2_effective_mesh_probe_h{tag}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
np.savez_compressed(ROOT/'results'/f'm2b2_effective_mesh_probe_h{tag}.npz',points=p,tetrahedra=t)
print(json.dumps(out,indent=2))
