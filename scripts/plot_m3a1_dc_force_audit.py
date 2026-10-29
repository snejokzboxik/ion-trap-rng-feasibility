"""Plot fixed-voltage direct P2 and local fitted restoring-force cuts."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
from iontrap_rng.hessian import fit_local_hessian
files={'hDC=.15':'m2b2_dc_surface_refined_p2_fields_h0p15.npz',
       'hDC=.10':'m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',
       'hDC=.075':'m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz'}
vdc=.08697882845042744;r0=1e-6;Q=elementary_charge
fig,ax=plt.subplots(1,2,figsize=(10,4))
for axis,name in ((0,'x'),(2,'z')):
    coord=np.linspace(-.25,.25,121);xi=np.zeros((3,len(coord)));xi[axis]=coord
    for label,path in files.items():
        fp=ROOT/'results'/path;field=P2ScalarField(fp,'dc_potential');g=field.gradient(xi)
        raw=np.load(fp);fit=fit_local_hessian(raw['p2_dof_coordinates'],raw['dc_potential'],radius=.5,order=2)
        # Force is -Q*Vdc/r0 times the dimensionless potential gradient.
        fdir=-Q*vdc/r0*g[axis]*1e18
        fpred=-Q*vdc/r0*(fit.gradient[axis]+fit.hessian[axis]@xi)*1e18
        ax[axis//2].plot(coord,fdir,label=f'{label}, direct P2',lw=1.0)
        if label=='hDC=.075':ax[axis//2].plot(coord,fpred,'k--',label='hDC=.075, nodal quadratic fit',lw=1.4)
    ax[axis//2].set_xlabel(f'{name} / r0');ax[axis//2].set_ylabel('DC force (aN)')
    ax[axis//2].set_title(f'Fixed VDC={vdc:.5f} V/ring, {name}-axis');ax[axis//2].grid(alpha=.3);ax[axis//2].legend(fontsize=7)
fig.suptitle('M3a.1: DC restoring-force mesh and interpolation audit')
fig.tight_layout();out=ROOT/'results'/'figures'/'m3a1_dc_force_audit.png';out.parent.mkdir(exist_ok=True,parents=True);fig.savefig(out,dpi=180);print(out)
