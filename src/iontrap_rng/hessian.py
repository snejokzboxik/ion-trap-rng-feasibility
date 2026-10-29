"""Condition-aware local polynomial and recovered-gradient Hessian tools."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass(frozen=True)
class HessianFit:
    hessian: np.ndarray
    gradient: np.ndarray
    value_at_center: float
    rank: int
    condition_number: float
    relative_residual: float
    samples: int
    fit_order: int
    weighted: bool
    harmonic_constraint: bool


def _powers(max_order: int):
    return [(i, j, k) for degree in range(max_order + 1)
            for i in range(degree, -1, -1)
            for j in range(degree - i, -1, -1)
            for k in [degree - i - j]]


def fit_local_hessian(points: np.ndarray, values: np.ndarray, *, radius: float,
                      center: np.ndarray | None = None, order: int = 2,
                      weights: np.ndarray | None = None,
                      harmonic_quadratic: bool = False) -> HessianFit:
    """Fit a locally scaled Taylor polynomial with SVD-backed least squares.

    Coordinates are centered and divided by `radius`; returned derivatives
    are in the input coordinate units. Harmonic mode is restricted to a
    quadratic, trace-free Hessian and must be compared with the unconstrained fit.
    """
    p = np.asarray(points, float)
    v = np.asarray(values, float)
    c0 = np.zeros(3) if center is None else np.asarray(center, float)
    if p.ndim != 2 or p.shape[1] != 3 or v.shape != (len(p),):
        raise ValueError("points must be (n,3) and values (n,)")
    if radius <= 0 or order < 1:
        raise ValueError("radius and order must be positive")
    use = np.linalg.norm(p-c0, axis=1) <= radius
    q = (p[use]-c0)/radius
    y = v[use]
    if harmonic_quadratic and order != 2:
        raise ValueError("harmonic_quadratic currently supports order=2 only")
    if harmonic_quadratic:
        # c + g·q + .5*a(qx²-qz²)+.5*b(qy²-qz²)+cxy+cxz+cyz
        X = np.column_stack((np.ones(len(q)), q, .5*(q[:,0]**2-q[:,2]**2),
                             .5*(q[:,1]**2-q[:,2]**2), q[:,0]*q[:,1],
                             q[:,0]*q[:,2], q[:,1]*q[:,2]))
        terms = None
    else:
        terms = _powers(order)
        X = np.column_stack([q[:,0]**i*q[:,1]**j*q[:,2]**k for i,j,k in terms])
    if weights is None:
        sw = np.ones(len(p))[use]
        is_weighted = False
    else:
        w = np.asarray(weights, float)
        if w.shape != (len(p),) or np.any(w < 0):
            raise ValueError("weights must be nonnegative and have shape (n,)")
        sw = np.sqrt(w[use])
        is_weighted = True
    good = sw > 0
    Xw = X[good]*sw[good,None]
    yw = y[good]*sw[good]
    coef, _, rank, sing = np.linalg.lstsq(Xw, yw, rcond=None)
    cond = float(sing[0]/sing[-1]) if len(sing) and sing[-1] > 0 else float("inf")
    fit = X@coef
    denom = np.linalg.norm(yw)
    residual = float(np.linalg.norm((fit[good]-y[good])*sw[good]) / denom) if denom else float("nan")
    g = np.zeros(3)
    H = np.zeros((3,3))
    if harmonic_quadratic:
        g[:] = coef[1:4]/radius
        a,b=coef[4:6]
        H=np.array([[a,coef[6],coef[7]],[coef[6],b,coef[8]],[coef[7],coef[8],-a-b]])/radius**2
        value=float(coef[0])
    else:
        lookup={term:i for i,term in enumerate(terms)}
        value=float(coef[lookup[(0,0,0)]])
        for i in range(3):
            exponent=[0,0,0]; exponent[i]=1
            g[i]=coef[lookup[tuple(exponent)]]/radius
            exponent2=[0,0,0]; exponent2[i]=2
            if tuple(exponent2) in lookup:
                H[i,i]=2*coef[lookup[tuple(exponent2)]]/radius**2
        for i,j in ((0,1),(0,2),(1,2)):
            e=[0,0,0];e[i]=e[j]=1
            if tuple(e) in lookup:
                H[i,j]=H[j,i]=coef[lookup[tuple(e)]]/radius**2
    return HessianFit(H,g,value,int(rank),cond,residual,int(np.count_nonzero(good)),order,is_weighted,harmonic_quadratic)


def tetra_element_gradients(points: np.ndarray, tetrahedra: np.ndarray,
                            potential: np.ndarray):
    """Recover exact piecewise-constant P1 gradients and tetra volumes."""
    tetp=np.asarray(points)[np.asarray(tetrahedra)]
    vals=np.asarray(potential)[np.asarray(tetrahedra)]
    M=np.concatenate((tetp,np.ones((*tetp.shape[:2],1))),axis=2)
    coeff=np.linalg.solve(M,vals[...,None])[...,0]
    v0=tetp[:,1]-tetp[:,0];v1=tetp[:,2]-tetp[:,0];v2=tetp[:,3]-tetp[:,0]
    det=np.einsum("ei,ei->e",v0,np.cross(v1,v2))
    volume=np.abs(det)/6
    return coeff[:,:3],volume,tetp.mean(axis=1)
