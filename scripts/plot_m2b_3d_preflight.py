"""Make deterministic XY field slices and a schematic of the 3-D preflight."""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from matplotlib.path import Path as MPath
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from iontrap_rng.readout_fem import truncated_electrode_polygon

ROOT=Path(__file__).resolve().parents[1]
d=np.load(ROOT/"results"/"m2b_3d_fields.npz")
p=d["points"]; zidx=np.abs(p[:,2])<0.08
xg=np.linspace(-2.0,2.0,240); yg=np.linspace(-2.0,2.0,240)
X,Y=np.meshgrid(xg,yg); Q=np.column_stack((X.ravel(),Y.ravel(),np.zeros(X.size)))
mask=np.zeros(len(Q),dtype=bool)
for b in ("right","left","top","bottom"):
    poly=truncated_electrode_polygon(b,1.,.10,1.5,80)
    mask |= MPath(poly.T).contains_points(Q[:,:2])

fig,axs=plt.subplots(1,2,figsize=(10.5,4.5),constrained_layout=True)
for ax,field,title in ((axs[0],d["rf"],"RF potential slice, z=0"),(axs[1],d["weighting"],"Differential weighting potential, z=0")):
    interp=LinearNDInterpolator(p[zidx],field[zidx],fill_value=np.nan)
    V=interp(Q).reshape(X.shape); V[:, :][mask.reshape(X.shape)]=np.nan
    cont=ax.contourf(X,Y,V,levels=25,cmap="RdBu_r")
    dx=xg[1]-xg[0]; dy=yg[1]-yg[0]
    dVdy,dVdx=np.gradient(V,dy,dx)
    stride=18
    ax.quiver(X[::stride,::stride],Y[::stride,::stride],-dVdx[::stride,::stride],-dVdy[::stride,::stride],
              color="k",alpha=.55,pivot="mid",scale=10)
    for b in ("right","left","top","bottom"):
        poly=truncated_electrode_polygon(b,1.,.10,1.5,80)
        ax.fill(poly[0],poly[1],color="0.25",zorder=3)
    ax.set(xlim=(-2,2),ylim=(-2,2),xlabel=r"$x/r_0$",ylabel=r"$y/r_0$",title=title,aspect="equal")
    fig.colorbar(cont,ax=ax,label="potential (normalized V)")
fig.suptitle("M2b ideal-conductor vacuum preflight (P1 tetrahedral mesh)")
fig.savefig(ROOT/"results"/"figures"/"m2b_3d_potential_slices.png",dpi=180)
plt.close(fig)

# Geometry illustration is explicitly schematic; plotted rails use the exact
# analytic hyperbola profile, while the solve used its polygonal approximation.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
fig=plt.figure(figsize=(7.5,6)); ax=fig.add_subplot(111,projection="3d")
zz=np.linspace(-5,5,80)
for branch in ("right","left","top","bottom"):
    poly=truncated_electrode_polygon(branch,1.,.10,1.5,80)
    for j in (0,len(poly.T)//2):
        ax.plot(np.full_like(zz,poly[0,j]),np.full_like(zz,poly[1,j]),zz,lw=1.2)
    ax.plot(poly[0],poly[1],np.full(poly.shape[1],-5),lw=.5,color="0.5")
    ax.plot(poly[0],poly[1],np.full(poly.shape[1],5),lw=.5,color="0.5")
theta=np.linspace(0,2*np.pi,120)
for zc in (-5.5,5.5):
    for r in (.5,1.5): ax.plot(r*np.cos(theta),r*np.sin(theta),np.full_like(theta,zc),lw=1.1)
ax.set(xlim=(-2,2),ylim=(-2,2),zlim=(-7,7),xlabel=r"$x/r_0$",ylabel=r"$y/r_0$",zlabel=r"$z/r_0$",title="Preflight terminal geometry (idealized)")
fig.tight_layout(); fig.savefig(ROOT/"results"/"figures"/"m2b_3d_geometry.png",dpi=180); plt.close(fig)
