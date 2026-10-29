"""Create lightweight M3a deterministic trajectory/pseudopotential figure."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge,atomic_mass
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
res=json.loads((ROOT/'results'/'m3a_deterministic_results.json').read_text(encoding='utf-8'))
rf=P2ScalarField(ROOT/res['field_identity']['field_npz'],'rf_potential')
dc=P2ScalarField(ROOT/res['field_identity']['field_npz'],'dc_potential')
par=res['physical_parameters'];r0=par['r0_m'];omega=par['omega_rad_s'];Vrf=par['Vrf_differential_peak_V'];Vdc=par['Vdc_each_ring_V']
Q=elementary_charge;M=39.96259098*atomic_mass
fig,ax=plt.subplots(1,3,figsize=(14,4.2))
for tr in res['trajectory_trace']:
    t=np.asarray(tr['tauf'])/np.pi;x=np.asarray(tr['xi'])*r0*1e9
    ax[0].plot(t,np.linalg.norm(x[:,:2],axis=1),label=f"{tr['case']}: radial")
    if tr['case']=='combined':ax[0].plot(t,x[:,2],label='combined: z',ls='--')
ax[0].set_xlabel('RF periods');ax[0].set_ylabel('Position (nm)');ax[0].set_title('DOP853 trajectories, phase 0');ax[0].legend(fontsize=8);ax[0].grid(alpha=.3)
coord=np.linspace(-.2,.2,101)
for axis,label,color,panel in [(0,'radial x','#2878b5',1),(2,'axial z','#d95f02',2)]:
    xi=np.zeros((3,len(coord)));xi[axis]=coord
    g=rf.gradient(xi)/r0;phidc=dc.potential(xi)
    Us=Q**2*Vrf**2*np.sum(g*g,axis=0)/(4*M*omega**2)+Q*Vdc*phidc
    g0=rf.gradient(np.zeros((3,1)))[:,0]/r0
    center=Q**2*Vrf**2*np.dot(g0,g0)/(4*M*omega**2)+Q*Vdc*dc.potential(np.zeros((3,1)))[0]
    ax[panel].plot(coord,((Us-center)/Q)*1e6,label=label,color=color)
    ax[panel].set_xlabel('Displacement / r0');ax[panel].set_ylabel('ΔU (µeV)');ax[panel].set_title(f'Local pseudopotential: {label}');ax[panel].grid(alpha=.3)
fig.suptitle('M3a: deterministic 40Ca+ in saved M2b.2 FEM fields')
fig.tight_layout();out=ROOT/'results'/'figures'/'m3a_deterministic.png';out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=180);print(out)
