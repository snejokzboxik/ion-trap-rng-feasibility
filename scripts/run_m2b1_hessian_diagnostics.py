"""Diagnose M2b P1 Hessian postprocessing without changing the archived run."""
from __future__ import annotations
import json
from pathlib import Path
import sys

import numpy as np
from skfem import MeshTet, Basis, ElementTetP1

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from iontrap_rng.hessian import fit_local_hessian,tetra_element_gradients


def gradient_hessian_fit(centers,gradients,volumes,radius,order=1):
    use=np.linalg.norm(centers,axis=1)<radius
    q=centers[use]/radius; g=gradients[use]; w=np.sqrt(volumes[use])
    if order==1:
        X=np.column_stack((np.ones(len(q)),q))
    else:
        X=np.column_stack((np.ones(len(q)),q,q[:,0]**2,q[:,1]**2,q[:,2]**2,q[:,0]*q[:,1],q[:,0]*q[:,2],q[:,1]*q[:,2]))
    Xw=X*w[:,None]; Yw=g*w[:,None]
    c,res,rank,sing=np.linalg.lstsq(Xw,Yw,rcond=None)
    H=c[1:4,:].T/radius
    rel=float(np.linalg.norm(Xw@c-Yw)/max(np.linalg.norm(Yw),1e-300))
    cond=float(sing[0]/sing[-1]) if sing[-1]>0 else float("inf")
    return {"Hessian":H,"samples":int(use.sum()),"rank":int(rank),"condition_number":cond,"relative_residual":rel}


def polyhedral_gauss_flux(points,tets,gradient,centers,radius):
    inside=np.linalg.norm(centers,axis=1)<radius
    face_map={}
    for e,tet in enumerate(tets):
        for omit in range(4):
            face=tuple(sorted(np.delete(tet,omit).tolist()))
            face_map.setdefault(face,[]).append(e)
    flux=absflux=0.
    for face,adj in face_map.items():
        sel=[e for e in adj if inside[e]]
        if len(sel)!=1 or (len(adj)==2 and inside[adj[0]]==inside[adj[1]]): continue
        e=sel[0]; xyz=points[list(face)]; normal=np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0]); area=np.linalg.norm(normal)/2
        normal=normal/(2*area)
        fc=xyz.mean(axis=0)
        if np.dot(normal,fc-centers[e])<0: normal=-normal
        F=-gradient[e]
        val=float(np.dot(F,normal)*area);flux+=val;absflux+=abs(val)
    return {"radius_r0":radius,"net_flux_normalized":flux,"sum_absolute_face_flux":absflux,
            "relative_net_flux":abs(flux)/max(absflux,1e-300)}


def interpolate_p1_in_mesh(points,tets,potential,q,mesh,finder):
    elem=finder(q[:,0],q[:,1],q[:,2])
    mapping=mesh._mapping()
    xi=mapping.invF(q.T[:,None,:],elem)[:,0,:]
    bary=np.vstack((1-xi.sum(axis=0),xi)).T
    # MeshTet may reorient local vertex ordering; use its own topology with
    # the reference coordinates returned by the same mapping.
    nodes=mesh.t[:,elem].T
    return np.sum(bary*potential[nodes],axis=1)


def symmetry_audit(potential,points,tets,kind,mesh,finder):
    interpolator=Basis(mesh,ElementTetP1()).interpolator(potential)
    ax=np.linspace(-.45,.45,9);X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij")
    q=np.column_stack((X.ravel(),Y.ravel(),Z.ravel()))
    # Remove points too near electrodes; this audit is an interior symmetry check.
    q=q[np.linalg.norm(q,axis=1)<.70]
    base=interpolator(q.T)
    rows={}
    for axis in range(3):
        mirrored=q.copy();mirrored[:,axis]*=-1
        other=interpolator(mirrored.T)
        parity=-1 if (kind=="weighting" and axis==0) else 1
        good=np.isfinite(base)&np.isfinite(other)
        rows[f"axis_{'xyz'[axis]}"]=float(np.sqrt(np.mean((base[good]-parity*other[good])**2))/max(np.sqrt(np.mean(base[good]**2)),1e-30))
    return rows


def gaussian_sphere_flux(points,tets,element_grad,radius,mesh,finder,n_mu=12,n_phi=24):
    """Integrate P1 FEM E·n on a true spherical Gaussian surface by sampling elements."""
    from numpy.polynomial.legendre import leggauss
    mu,wmu=leggauss(n_mu); phi=2*np.pi*np.arange(n_phi)/n_phi
    M,P=np.meshgrid(mu,phi,indexing="ij")
    s=np.sqrt(1-M*M)
    q=np.column_stack((radius*s.ravel()*np.cos(P).ravel(),radius*s.ravel()*np.sin(P).ravel(),radius*M.ravel()))
    n=q/radius
    e=finder(q[:,0],q[:,1],q[:,2]); E=-element_grad[e]
    weights=np.repeat(wmu,n_phi)*(2*np.pi/n_phi)*radius**2
    f=np.sum(E*n,axis=1)
    signed=float(np.sum(f*weights));absolute=float(np.sum(np.abs(f)*weights))
    return {"radius_r0":radius,"quadrature_points":int(len(q)),"net_flux":signed,
            "integral_abs_normal_flux":absolute,"relative_net_flux":abs(signed)/max(absolute,1e-300)}


def main():
    d=np.load(ROOT/"results"/"m2b_3d_fields.npz")
    p,t=d["points"],d["tetrahedra"]
    mesh=MeshTet(p.T,t.T); finder=mesh.element_finder()
    out={"source":"archived results/m2b_3d_fields.npz from M2b coarse P1 tetra mesh; no old results overwritten",
         "mesh":{"nodes":len(p),"tetrahedra":len(t)},"unconstrained_local_fits":[],"gradient_recovery":[],"gauss_flux":[],"symmetry":{},
         "direct_elementwise_P1_Hessian":"identically zero within each affine tetra; it cannot be used as physical Hessian estimator"}
    # Lumped nodal volume weights reduce bias from the center-refined vertex cloud.
    dummy=np.zeros(len(p));_,vol,cent=tetra_element_gradients(p,t,dummy)
    nodalw=np.zeros(len(p))
    np.add.at(nodalw,t.ravel(),np.repeat(vol/4,4))
    for name in ("rf","dc","weighting"):
        u=d[name]
        fits=[]
        for r in (.35,.50,.75,1.0):
            for order in (2,3):
                for weighted in (False,True):
                    f=fit_local_hessian(p,u,radius=r,order=order,weights=nodalw if weighted else None)
                    # Raw-coordinate conditioning diagnostic on exactly the same selected nodes.
                    use=np.linalg.norm(p,axis=1)<=r; q=p[use]
                    from iontrap_rng.hessian import _powers
                    powers=_powers(order)
                    Xraw=np.column_stack([q[:,0]**i*q[:,1]**j*q[:,2]**k for i,j,k in powers])
                    rawcond=float(np.linalg.cond(Xraw))
                    fitrow={"radius_r0":r,"order":order,"weighted_by_lumped_volume":weighted,
                            "samples":f.samples,"rank":f.rank,"condition_normalized":f.condition_number,
                            "condition_raw_coordinates":rawcond,"relative_fit_residual":f.relative_residual,
                            "Hessian":f.hessian.tolist(),"trace":float(np.trace(f.hessian))}
                    fits.append(fitrow)
            if r in (.5,1.0):
                hfit=fit_local_hessian(p,u,radius=r,order=2,weights=nodalw,harmonic_quadratic=True)
                fits.append({"radius_r0":r,"order":"quadratic harmonic-constrained","weighted_by_lumped_volume":True,
                             "samples":hfit.samples,"rank":hfit.rank,"condition_normalized":hfit.condition_number,
                             "relative_fit_residual":hfit.relative_residual,"Hessian":hfit.hessian.tolist(),"trace":float(np.trace(hfit.hessian))})
        eg,ev,ec=tetra_element_gradients(p,t,u)
        gradrows=[]
        for r in (.5,.75,1.0):
            row=gradient_hessian_fit(ec,eg,ev,r,1)
            gradrows.append({"radius_r0":r,"samples":row["samples"],"rank":row["rank"],"condition_number":row["condition_number"],
                             "relative_fit_residual":row["relative_residual"],"Hessian":row["Hessian"].tolist(),"trace":float(np.trace(row["Hessian"]))})
            if r in (.5,.75):
                out["gauss_flux"].append({"field":name,**gaussian_sphere_flux(p,t,eg,r,mesh,finder)})
        out["unconstrained_local_fits"].append({"field":name,"fits":fits})
        out["gradient_recovery"].append({"field":name,"method":"volume-weighted least-squares recovery of piecewise-constant P1 element gradients", "fits":gradrows})
        out["symmetry"][name]=symmetry_audit(u,p,t,name,mesh,finder)
    dest=ROOT/"results"/"m2b1_old_grid_diagnostics.json"
    dest.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(dest)
    for item in out["unconstrained_local_fits"]:
        print(item["field"],[(f["radius_r0"],f["order"],f["samples"],round(f["condition_normalized"],2),round(f["trace"],5)) for f in item["fits"]])
    print("gradient recovery",[(x["field"],[(f["radius_r0"],f["samples"],np.round(f["Hessian"],4).tolist(),round(f["trace"],4)) for f in x["fits"]]) for x in out["gradient_recovery"]])

if __name__=="__main__":main()
