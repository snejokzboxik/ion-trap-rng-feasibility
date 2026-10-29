"""Bounded Gaussian-sphere flux audit by finite differences of saved P2 fields."""
import json,sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from skfem import Basis,MeshTet,ElementTetP2
ROOT=Path(__file__).resolve().parents[1]
def flux(filename,r=.5,nmu=8,nphi=16,eps=.02,batch=16):
    d=np.load(filename);mesh=MeshTet(d['points'].T,d['tetrahedra'].T);basis=Basis(mesh,ElementTetP2());u=d['dc_potential']
    interp=basis.interpolator(u);mu,wm=leggauss(nmu);phi=2*np.pi*np.arange(nphi)/nphi;M,P=np.meshgrid(mu,phi,indexing='ij');s=np.sqrt(1-M*M)
    n=np.column_stack((s.ravel()*np.cos(P).ravel(),s.ravel()*np.sin(P).ravel(),M.ravel()))
    w=np.repeat(wm,nphi)*(2*np.pi/nphi)*r*r;dr=[]
    for i in range(0,len(n),batch):
        nn=n[i:i+batch];qp=(r+eps)*nn;qm=(r-eps)*nn
        vals=interp(np.concatenate((qp.T,qm.T),axis=1));m=len(nn)
        dr.extend(((vals[:m]-vals[m:])/(2*eps)).tolist())
    radial=np.asarray(dr);signed=-float(np.dot(radial,w));absolute=float(np.dot(np.abs(radial),w))
    return {'file':Path(filename).name,'radius_r0':r,'n_mu':nmu,'n_phi':nphi,'points':len(n),'fd_step_r0':eps,
            'net_flux_normalized':signed,'integral_abs_normal_field':absolute,
            'relative_gauss_residual':abs(signed)/max(absolute,1e-300)}
if __name__=='__main__':
    inputs=[('results/m2b2_dc_surface_refined_p2_fields_h0p15.npz',8,16,.02),
            ('results/m2b2_dc_surface_refined_p2_fields_h0p15.npz',12,24,.02),
            ('results/m2b2_dc_surface_refined_p2_fields_h0p15.npz',12,24,.01),
            ('results/m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',8,16,.02),
            ('results/m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',12,24,.02),
            ('results/m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',12,24,.01),
            ('results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz',12,24,.02)]
    out={'method':'central radial finite difference of P2 interpolated potential; bounded batches of 16; normalized model units',
         'scope':'DC field at one sphere radius; finite-difference/quadrature diagnostic, not a high-order flux estimator',
         'results':[flux(ROOT/f,r=.5,nmu=nm,nphi=np_,eps=e) for f,nm,np_,e in inputs]}
    (ROOT/'results'/'m2b2_gaussian_flux.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
