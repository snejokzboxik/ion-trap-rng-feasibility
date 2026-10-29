"""Mesh-only central refinement probes for M3a.2; no FEM assembly or solve."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_m2b_3d_preflight import make_mesh

def stats(a):
    a=np.asarray(a,dtype=float)
    return {"count":int(a.size),"min":float(a.min()),"p10":float(np.percentile(a,10)),
            "median":float(np.median(a)),"p90":float(np.percentile(a,90)),
            "p99":float(np.percentile(a,99)),"max":float(a.max())}

def tetra_quality(p,t):
    q=p[t]
    e2=[]
    for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
        e2.append(np.sum((q[:,i]-q[:,j])**2,axis=1))
    sum_e2=np.sum(e2,axis=0)
    vol=np.abs(np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0])))/6.
    mr=12.*np.power(3.*vol,2./3.)/sum_e2
    return vol,mr

rows=[]
for hc in (.18,.09,.06):
    p,t,meta=make_mesh(.75,hc,dc_surface_h=.075,rf_surface_h=.30,axial_end_h=.30,return_metadata=True)
    edges=np.unique(np.sort(np.concatenate([t[:,[0,1]],t[:,[0,2]],t[:,[0,3]],t[:,[1,2]],t[:,[1,3]],t[:,[2,3]]]),axis=1),axis=0)
    ep=(p[edges[:,0]]+p[edges[:,1]])/2
    el=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
    r=np.linalg.norm(ep,axis=1)
    vol,mr=tetra_quality(p,t)
    cent=t[np.linalg.norm(p[t].mean(axis=1),axis=1)<.5]
    _,cmr=tetra_quality(p,cent)
    # A strict interior central region avoids counting surface-refined edges.
    region={str(R):stats(el[r<R]) for R in (.15,.30,.50)}
    out={"hcenter_target_r0":hc,"hmax_r0":.75,"dc_surface_h_target_r0":.075,
         "rf_surface_h_target_r0":.30,"axial_end_h_target_r0":.30,
         "mesh_size_min_effective_r0":meta['mesh_size_min_effective_r0'],
         "mesh_size_max_r0":meta['mesh_size_max_r0'],"vertices":len(p),"tetrahedra":len(t),
         "unique_edges":len(edges),"estimated_P2_dofs":int(len(p)+len(edges)),
         "all_unique_edge_lengths_r0":stats(el),"central_edge_lengths_by_radius_r0":region,
         "tetrahedron_volume_r0_cubed":stats(vol),"tetrahedron_mean_ratio_quality_0_1":stats(mr),
         "central_tet_quality_centroid_radius_lt_0p5r0":stats(cmr),
         "geometry_metadata":meta,
         "purpose":"mesh-only probe; no P2 matrix assembly or field solve"}
    tag=str(hc).replace('.','p')
    np.savez_compressed(ROOT/'results'/f'm3a2_center_mesh_hc{tag}.npz',points=p,tetrahedra=t)
    (ROOT/'results'/f'm3a2_center_mesh_hc{tag}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    rows.append(out)
print(json.dumps(rows,indent=2))
