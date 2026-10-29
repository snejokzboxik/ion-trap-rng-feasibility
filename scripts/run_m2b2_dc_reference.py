"""Small exact axial-quadrupole manufactured test for derivative estimators."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from skfem import Basis,MeshTet,ElementTetP2,asm,solve,condense
from skfem.models.poisson import laplace
ROOT=Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.hessian import fit_local_hessian

A_AXIAL=1.0
H=A_AXIAL*np.diag([-1.,-1.,2.])

def phi(x): return .5*np.einsum('i...,ij,j...->...',x,H,x)
def grad(x): return np.einsum('ij,j...->i...',H,x)

def local_quadratic_hessian(points,values):
    center=points.mean(axis=0);scale=max(np.max(np.linalg.norm(points-center,axis=1)),1e-15)
    q=(points-center)/scale
    X=np.column_stack((np.ones(len(q)),q[:,0],q[:,1],q[:,2],.5*q[:,0]**2,.5*q[:,1]**2,.5*q[:,2]**2,
                       q[:,0]*q[:,1],q[:,0]*q[:,2],q[:,1]*q[:,2]))
    c,_,rank,s=np.linalg.lstsq(X,values,rcond=None)
    h=np.array([[c[4],c[7],c[8]],[c[7],c[5],c[9]],[c[8],c[9],c[6]]])/scale**2
    return h,{'rank':int(rank),'condition':float(s[0]/s[-1]),'residual':float(np.linalg.norm(X@c-values)/np.linalg.norm(values))}

def run(n=8):
    axis=np.linspace(-1,1,n+1);mesh=MeshTet.init_tensor(axis,axis,axis)
    basis=Basis(mesh,ElementTetP2(),intorder=8);K=asm(laplace,basis)
    fixed=basis.get_dofs().all();u=solve(*condense(K,np.zeros(basis.N),x=phi(basis.doflocs),D=fixed))
    # Direct second derivative of the local P2 polynomial in the central tetrahedron.
    centers=mesh.p[:,mesh.t].mean(axis=1).T;ie=int(np.argmin(np.linalg.norm(centers,axis=1)))
    edofs=basis.element_dofs[:,ie];Hloc,localdiag=local_quadratic_hessian(basis.doflocs[:,edofs].T,u[edofs])
    # Full regional SVD/QR least-squares fits to nodal P2 values, both free and trace-constrained.
    point=basis.doflocs.T;values=u
    valuefits={}
    for radius in (.25,.5,.75):
        row={}
        for constrained in (False,True):
            f=fit_local_hessian(point,values,radius=radius,order=2,harmonic_quadratic=constrained)
            row['harmonic_constrained' if constrained else 'unconstrained']={
                'samples':f.samples,'rank':f.rank,'condition':f.condition_number,'relative_fit_residual':f.relative_residual,
                'Hessian':f.hessian.tolist(),'absolute_frobenius_error':float(np.linalg.norm(f.hessian-H)),
                'relative_frobenius_error':float(np.linalg.norm(f.hessian-H)/np.linalg.norm(H)),'trace':float(np.trace(f.hessian))}
        valuefits[str(radius)]=row
    # Independent recovered-gradient regression at P2 quadrature points.
    X=basis.global_coordinates();dx=basis.dx;g=basis.interpolate(u).grad
    xyz=np.moveaxis(X,0,-1).reshape(-1,3);gv=np.moveaxis(g,0,-1).reshape(-1,3);weights=np.abs(dx.ravel())
    use=np.linalg.norm(xyz,axis=1)<.5; q=xyz[use]/.5;sw=np.sqrt(weights[use]);M=np.column_stack((np.ones(len(q)),q));Mw=M*sw[:,None]
    c,_,rank,s=np.linalg.lstsq(Mw,gv[use]*sw[:,None],rcond=None);Hgrad=c[1:4,:].T/.5
    gradfit={'samples':int(use.sum()),'rank':int(rank),'condition':float(s[0]/s[-1]),
             'relative_gradient_fit_residual':float(np.linalg.norm(Mw@c-gv[use]*sw[:,None])/np.linalg.norm(gv[use]*sw[:,None])),
             'Hessian':Hgrad.tolist(),'absolute_frobenius_error':float(np.linalg.norm(Hgrad-H)),
             'relative_frobenius_error':float(np.linalg.norm(Hgrad-H)/np.linalg.norm(H)),'trace':float(np.trace(Hgrad))}
    # Central finite differences through the actual P2 FEM interpolator.
    interp=basis.interpolator(u)
    def evalp(q):return interp(np.asarray(q,float).reshape(3,-1)).reshape(-1)
    f0=float(evalp([[0,0,0]])[0]);fds=[]
    for d in (.2,.1,.05,.025):
        hf=np.empty((3,3))
        for i in range(3):
            e=np.zeros(3);e[i]=d
            hf[i,i]=(evalp(e)-2*f0+evalp(-e))[0]/d**2
        for i,j in ((0,1),(0,2),(1,2)):
            ei=np.zeros(3);ej=np.zeros(3);ei[i]=d;ej[j]=d
            hf[i,j]=hf[j,i]=(evalp(ei+ej)-evalp(ei-ej)-evalp(-ei+ej)+evalp(-ei-ej))[0]/(4*d**2)
        fds.append({'step':d,'Hessian':hf.tolist(),'absolute_frobenius_error':float(np.linalg.norm(hf-H)),
                    'relative_frobenius_error':float(np.linalg.norm(hf-H)/np.linalg.norm(H)),'trace':float(np.trace(hf))})
    val=basis.interpolate(u);exact=phi(X);fielderr=np.sqrt(np.sum(np.sum((g-grad(X))**2,axis=0)*dx)/np.sum(np.sum(grad(X)**2,axis=0)*dx))
    return {'problem':'Phi=A[z^2-(x^2+y^2)/2]','A':A_AXIAL,'analytic_Hessian':H.tolist(),'analytic_trace':float(np.trace(H)),
            'mesh':{'n_per_axis':n,'tetrahedra':int(mesh.nelements),'P2_dofs':int(basis.N)},
            'global_relative_potential_L2':float(np.sqrt(np.sum((val- exact)**2*dx)/np.sum(exact**2*dx))),
            'global_relative_field_L2':float(fielderr),
            'direct_central_element_P2_derivative':{'element':ie,'Hessian':Hloc.tolist(),'diagnostics':localdiag,
                'absolute_frobenius_error':float(np.linalg.norm(Hloc-H)),'relative_frobenius_error':float(np.linalg.norm(Hloc-H)/np.linalg.norm(H)),'trace':float(np.trace(Hloc))},
            'local_polynomial_value_fits':valuefits,'recovered_gradient_fit':gradfit,'finite_difference_on_P2_interpolant':fds,
            'interpretation':'Exact harmonic manufactured problem; validates estimators for a quadratic field only, not the finite DC ring geometry.'}

if __name__=='__main__':
    out=run();dest=ROOT/'results'/'m2b2_dc_reference.json';dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
