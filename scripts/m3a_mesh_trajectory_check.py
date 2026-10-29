"""Two-case, two-mesh trajectory sensitivity check for M3a."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge,atomic_mass
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
Q=elementary_charge;M=39.96259098*atomic_mass;r0=1e-6;Vrf=10.;q=.2;k_rf=.9998641069136932
omega=np.sqrt(2*Q*Vrf*k_rf/(M*r0*r0*q));gamma=4*Q/(M*r0*r0*omega*omega)
Vdc=.08697882845042744;T=20*np.pi;N=1600;dt=T/N
inputs={'h010':'m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz',
        'h0075':'m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz'}
cases={'near_center':np.array([.02,0.,0.]),'combined':np.array([.05,.04,.08])}
results={}
for mesh,filename in inputs.items():
    rf=P2ScalarField(ROOT/'results'/filename,'rf_potential');dc=P2ScalarField(ROOT/'results'/filename,'dc_potential')
    for case,x0 in cases.items():
        y=np.r_[x0,np.zeros(3)];maxr=np.linalg.norm(x0)
        def acc(tau,x):
            xx=x[:,None]
            g=Vrf*np.cos(2*tau)*rf.gradient(xx)[:,0]+Vdc*dc.gradient(xx)[:,0]
            return -gamma*g
        tau=0.
        for _ in range(N):
            x=y[:3];v=y[3:]
            k1x=v;k1v=acc(tau,x)
            k2x=v+dt*k1v/2;k2v=acc(tau+dt/2,x+dt*k1x/2)
            k3x=v+dt*k2v/2;k3v=acc(tau+dt/2,x+dt*k2x/2)
            k4x=v+dt*k3v;k4v=acc(tau+dt,x+dt*k3x)
            y[:3]+=dt*(k1x+2*k2x+2*k3x+k4x)/6
            y[3:]+=dt*(k1v+2*k2v+2*k3v+k4v)/6
            tau+=dt;maxr=max(maxr,float(np.linalg.norm(y[:3])))
            if maxr>.30:break
        results[f'{mesh}_{case}']={'mesh':mesh,'field':filename,'case':case,'rf_phase_rad':0.,
            'completed_rf_periods':tau/np.pi,'final_xi':y[:3].tolist(),'max_radius_r0':float(maxr)}
out={'method':'fixed RK4, 80 steps per RF period, 20 periods or stop at radius 0.30r0',
     'drive':{'Vrf_differential_peak_V':Vrf,'q_target':q,'Vdc_per_ring_V':Vdc,'f_rf_Hz':omega/(2*np.pi)},
     'results':results,'limitations':['DC-ring surface h is the only mesh refinement; center h is fixed',
       'short observation (20 RF periods); not axial-cycle or lifetime evidence',
       'P2 gradient is elementwise discontinuous; comparison estimates discretization sensitivity only']}
(ROOT/'results'/'m3a_mesh_trajectory_check.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
