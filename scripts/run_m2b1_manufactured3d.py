"""Independent 3-D harmonic manufactured-solution test for P1/P2 tetrahedra."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from skfem import Basis,MeshTet,ElementTetP1,ElementTetP2,asm,solve,condense
from skfem.models.poisson import laplace

ROOT=Path(__file__).resolve().parents[1]
H=np.array([[1.0,.20,.10],[.20,-.40,.15],[.10,.15,-.60]])
assert abs(np.trace(H)) < 1e-14

def phi(x):
    return .5*np.einsum('i...,ij,j...->...',x,H,x)

def grad_exact(x):return np.einsum('ij,j...->i...',H,x)

def element_hessian_rms(mesh,basis,u):
    mats=[];volumes=[]
    for edofs in basis.element_dofs.T:
        xyz=basis.doflocs[:,edofs].T
        cen=xyz.mean(axis=0)
        scale=max(np.max(np.linalg.norm(xyz-cen,axis=1)),1e-15)
        q=(xyz-cen)/scale
        X=np.column_stack((np.ones(len(q)),q[:,0],q[:,1],q[:,2],
                           .5*q[:,0]**2,.5*q[:,1]**2,.5*q[:,2]**2,
                           q[:,0]*q[:,1],q[:,0]*q[:,2],q[:,1]*q[:,2]))
        c=np.linalg.lstsq(X,u[edofs],rcond=None)[0]
        h=np.array([[c[4],c[7],c[8]],[c[7],c[5],c[9]],[c[8],c[9],c[6]]])/scale**2
        mats.append(h)
        verts=mesh.p[:,mesh.t[:,len(mats)-1]].T
        a,b,cross=verts[1]-verts[0],verts[2]-verts[0],verts[3]-verts[0]
        volumes.append(abs(np.dot(a,np.cross(b,cross)))/6)
    mats=np.asarray(mats);volumes=np.asarray(volumes)
    err=np.sqrt(np.sum(volumes[:,None,None]*(mats-H)**2)/np.sum(volumes))/np.linalg.norm(H)
    return float(err),float(np.max(np.linalg.norm(mats-H,axis=(1,2))/np.linalg.norm(H,axis=(0,1))))

def run(n,order):
    axis=np.linspace(-1,1,n+1);mesh=MeshTet.init_tensor(axis,axis,axis)
    elem=ElementTetP1() if order==1 else ElementTetP2()
    basis=Basis(mesh,elem,intorder=8);A=asm(laplace,basis);D=basis.get_dofs().all()
    bc=phi(basis.doflocs);u=solve(*condense(A,np.zeros(basis.N),x=bc,D=D))
    X=basis.global_coordinates();dx=basis.dx
    e=basis.interpolate(u).grad-grad_exact(X)
    eg=float(np.sqrt(np.sum(np.sum(e*e,axis=0)*dx)/np.sum(np.sum(grad_exact(X)**2,axis=0)*dx)))
    val=basis.interpolate(u); exact=phi(X)
    ep=float(np.sqrt(np.sum((val-exact)**2*dx)/np.sum(exact**2*dx)))
    eh,emax=element_hessian_rms(mesh,basis,u)
    return {"n_cells_per_axis":int(n),"order":int(order),"dofs":int(basis.N),"tetrahedra":int(mesh.nelements),
            "global_relative_potential_L2":ep,"global_relative_gradient_L2":eg,
            "element_recovered_relative_hessian_volume_RMS":eh,
            "max_element_hessian_relative_error":emax,"analytic_hessian":H.tolist()}

def main():
    rows=[run(n,p) for n in (2,4,8) for p in (1,2)]
    payload={"scope":"independent 3-D harmonic manufactured solution, cube [-1,1]^3, exact Dirichlet data",
             "potential":"0.5*x^T*H*x","Hessian_trace":float(np.trace(H)),"mesh":"structured tetrahedral cube",
             "limits":["no physical electrodes", "P2 exactness applies only to this quadratic harmonic reference"],"results":rows}
    (ROOT/"results"/"m2b1_manufactured3d.json").write_text(json.dumps(payload,indent=2),encoding="utf-8")
    print(json.dumps(payload,indent=2))

if __name__=="__main__":main()
