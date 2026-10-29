"""Generate one-octant vacuum mesh for the symmetric RF/DC patterns.

Only mesh generation/inspection occurs here; no FEM system is assembled.
All dimensions are normalized by r0.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
import gmsh
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.readout_fem import truncated_electrode_polygon

RAIL_L=5.; BOX=(6.,6.,7.); RING_Z=5.5; RING_T=.2; RIN=.5; ROUT=1.5
BRANCHES=('right','left','top','bottom')

def _stats(a):
    a=np.asarray(a,float)
    return {'count':int(a.size),'min':float(a.min()),'median':float(np.median(a)),
            'p90':float(np.percentile(a,90)),'p99':float(np.percentile(a,99)),'max':float(a.max())}

def _quality(p,t):
    q=p[t]; l2=[]
    for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
        l2.append(np.sum((q[:,i]-q[:,j])**2,axis=1))
    v=np.abs(np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0])))/6
    mr=12*(3*v)**(2/3)/np.sum(l2,axis=0)
    return v,mr

def make_octant_mesh(hcenter=.025, hmax=.75, dc_h=.075, rf_h=.30, end_h=.30,
                     transition=.55, save_tag=None):
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
    gmsh.option.setNumber('Mesh.Algorithm3D',1)
    gmsh.option.setNumber('Mesh.MeshSizeMin',min(hcenter,dc_h,rf_h,end_h))
    gmsh.option.setNumber('Mesh.MeshSizeMax',hmax)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
    gmsh.model.add('m3a3_octant_hyperbolic_vacuum')
    occ=gmsh.model.occ
    box=occ.addBox(0,0,0,*BOX)
    solids=[]
    for branch in BRANCHES:
        poly=truncated_electrode_polygon(branch,1.,.10,1.5,80)
        pts=[occ.addPoint(float(x),float(y),-RAIL_L) for x,y in poly.T]
        lines=[occ.addLine(pts[i],pts[(i+1)%len(pts)]) for i in range(len(pts))]
        face=occ.addPlaneSurface([occ.addCurveLoop(lines)])
        ext=occ.extrude([(2,face)],0,0,2*RAIL_L)
        solids.append((3,ext[1][1]))
    for z0 in (-RING_Z-RING_T/2,RING_Z-RING_T/2):
        outer=occ.addCylinder(0,0,z0,0,0,RING_T,ROUT)
        inner=occ.addCylinder(0,0,z0-.01,0,0,RING_T+.02,RIN)
        ann,_=occ.cut([(3,outer)],[(3,inner)],removeObject=True,removeTool=True)
        solids.extend(ann)
    vac,_=occ.cut([(3,box)],solids,removeObject=True,removeTool=True)
    occ.synchronize()
    if not vac: raise RuntimeError('octant Boolean subtraction returned no vacuum volume')

    # Classify CAD boundary surfaces for separate RF-side, RF-end and DC sizing.
    polys={b:truncated_electrode_polygon(b,1.,.10,1.5,80).T for b in BRANCHES}
    from matplotlib.path import Path as MPath
    def poly_distance(xy,poly):
        a=poly;b=np.roll(a,-1,axis=0);v=b-a;w=xy[None,:]-a
        t=np.clip(np.sum(w*v,axis=1)/np.sum(v*v,axis=1),0.,1.)
        d=xy-(a+t[:,None]*v)
        return float(np.sqrt(np.min(np.sum(d*d,axis=1))))
    dc_surfs=[];rf_surfs=[];end_surfs=[]
    occ.synchronize()
    for dim,tag in gmsh.model.getEntities(2):
        bb=gmsh.model.getBoundingBox(dim,tag)
        xmin,ymin,zmin,xmax,ymax,zmax=bb
        cx,cy,cz=occ.getCenterOfMass(dim,tag)
        bounded=max(abs(xmin),abs(ymin),abs(xmax),abs(ymax))<=1.8
        if (5.34<=zmin and zmax<=5.66 and bounded): dc_surfs.append(tag)
        if (zmin>=-1e-6 and zmax<=RAIL_L+1e-6 and zmax-zmin>.25 and bounded):
            if min(poly_distance(np.array([cx,cy]),poly) for poly in polys.values())<.02:
                rf_surfs.append(tag)
        if (abs(zmin-RAIL_L)<.02 and abs(zmax-RAIL_L)<.02 and bounded):
            if any(MPath(poly).contains_point((cx,cy)) for poly in polys.values()): end_surfs.append(tag)
    if not dc_surfs or not rf_surfs:
        gmsh.finalize();raise RuntimeError(f'surface selection failed: DC={len(dc_surfs)}, RF={len(rf_surfs)}')

    # Smooth local field: characteristic size rises linearly from hcenter at
    # the origin to hmax by `transition`. The true achieved edge lengths are
    # measured after generation; they are not equated with these targets.
    p0=occ.addPoint(0,0,0);occ.synchronize()
    dist=gmsh.model.mesh.field.add('Distance');gmsh.model.mesh.field.setNumbers(dist,'PointsList',[p0])
    cent=gmsh.model.mesh.field.add('Threshold');gmsh.model.mesh.field.setNumber(cent,'InField',dist)
    gmsh.model.mesh.field.setNumber(cent,'SizeMin',hcenter);gmsh.model.mesh.field.setNumber(cent,'SizeMax',hmax)
    gmsh.model.mesh.field.setNumber(cent,'DistMin',.30);gmsh.model.mesh.field.setNumber(cent,'DistMax',transition)
    fields=[cent]
    for surfids,h,distmax in ((dc_surfs,dc_h,.45),(rf_surfs,rf_h,.45)):
        df=gmsh.model.mesh.field.add('Distance');gmsh.model.mesh.field.setNumbers(df,'SurfacesList',surfids)
        th=gmsh.model.mesh.field.add('Threshold');gmsh.model.mesh.field.setNumber(th,'InField',df)
        gmsh.model.mesh.field.setNumber(th,'SizeMin',h);gmsh.model.mesh.field.setNumber(th,'SizeMax',hmax)
        gmsh.model.mesh.field.setNumber(th,'DistMin',0.);gmsh.model.mesh.field.setNumber(th,'DistMax',distmax);fields.append(th)
    endpt=occ.addPoint(0,0,RAIL_L);occ.synchronize()
    ed=gmsh.model.mesh.field.add('Distance');gmsh.model.mesh.field.setNumbers(ed,'PointsList',[endpt])
    et=gmsh.model.mesh.field.add('Threshold');gmsh.model.mesh.field.setNumber(et,'InField',ed)
    gmsh.model.mesh.field.setNumber(et,'SizeMin',end_h);gmsh.model.mesh.field.setNumber(et,'SizeMax',hmax)
    gmsh.model.mesh.field.setNumber(et,'DistMin',0.);gmsh.model.mesh.field.setNumber(et,'DistMax',.45);fields.append(et)
    comb=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(comb,'FieldsList',fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(comb)
    gmsh.model.mesh.generate(3)
    node_tags,coords,_=gmsh.model.mesh.getNodes();allp=np.asarray(coords).reshape(-1,3)
    tagmap={int(tag):i for i,tag in enumerate(node_tags)};tet=[]
    ty,_,bl=gmsh.model.mesh.getElements(3)
    for et,b in zip(ty,bl):
        if et==4:tet.extend([[tagmap[int(i)] for i in row] for row in np.asarray(b).reshape(-1,4)])
    gmsh.finalize();tet=np.asarray(tet,int)
    used=np.unique(tet);remap=np.full(len(allp),-1,int);remap[used]=np.arange(len(used));p=allp[used];t=remap[tet]
    edge=np.unique(np.sort(np.concatenate([t[:,[0,1]],t[:,[0,2]],t[:,[0,3]],t[:,[1,2]],t[:,[1,3]],t[:,[2,3]]]),axis=1),axis=0)
    ep0=p[edge[:,0]];ep1=p[edge[:,1]];el=np.linalg.norm(ep0-ep1,axis=1)
    edge_regions={}
    for R in (.15,.30):
        fully_inside=(np.linalg.norm(ep0,axis=1)<R)&(np.linalg.norm(ep1,axis=1)<R)
        edge_regions[f'radius_lt_{R:.2f}']={'both_endpoints_inside':_stats(el[fully_inside]),
          'midpoint_inside_diagnostic':_stats(el[np.linalg.norm((ep0+ep1)/2,axis=1)<R])}
    vol,mr=_quality(p,t)
    q=p[t];local_edge_lengths=[]
    for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
        local_edge_lengths.append(np.linalg.norm(q[:,i]-q[:,j],axis=1))
    local_edge_lengths=np.column_stack(local_edge_lengths)
    tmid=q.mean(axis=1);tcent=np.linalg.norm(tmid,axis=1)
    tets_by_radius={f'radius_lt_{R:.2f}':_stats(local_edge_lengths[tcent<R].ravel()) for R in (.15,.30)}
    central=tcent<.5
    out={'model':'positive octant of unchanged M2b.2 geometry; vacuum only',
      'mesh_targets_r0':{'hcenter':hcenter,'hmax_outer':hmax,'dc_surface':dc_h,'rf_surface':rf_h,'rail_end':end_h,
       'central_size_plateau_radius':.30,'central_transition_end_radius':transition,
       'Mesh.MeshSizeMin':min(hcenter,dc_h,rf_h,end_h)},
      'cad_surface_entity_counts':{'dc':len(dc_surfs),'rf_side':len(rf_surfs),'rail_end_point_refinement':1},
      'vertices':len(p),'tetrahedra':len(t),'unique_edges':len(edge),'estimated_P2_dofs':int(len(p)+len(edge)),
      'edge_lengths_r0':{'global':_stats(el),'region_specific_edges':edge_regions,'tetrahedra_by_centroid_radius':tets_by_radius},
      'tetra_quality_mean_ratio_0_1':{'global':_stats(mr),'central_centroid_radius_lt_0p5':_stats(mr[central])},
      'tetra_volume_r0_cubed':_stats(vol),'purpose':'mesh-only; no P2 assembly or solve'}
    if save_tag:
        np.savez_compressed(ROOT/'results'/f'm3a3_octant_mesh_{save_tag}.npz',points=p,tetrahedra=t)
        (ROOT/'results'/f'm3a3_octant_mesh_{save_tag}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return p,t,out

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--hcenter',type=float,required=True);ap.add_argument('--tag',required=True)
    a=ap.parse_args();_,_,d=make_octant_mesh(hcenter=a.hcenter,save_tag=a.tag);print(json.dumps(d,indent=2))
