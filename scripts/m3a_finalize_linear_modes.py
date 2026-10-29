"""Recompute M3a linear-mode diagnostics from saved FEM Hessians."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge,atomic_mass
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results'/'m3a_deterministic_results.json';out=json.loads(p.read_text(encoding='utf-8'))
base=json.loads((ROOT/'results'/'m2b2_p2_effective_h0.75_hc0.18_dc0.075.json').read_text(encoding='utf-8'))
rfrow=base['fields']['rf'][1];dcrow=base['fields']['dc'][1]
Hrf0=np.asarray(rfrow['Hessian'],float);Hdc0=np.asarray(dcrow['Hessian'],float)
Hrf=.5*(Hrf0+Hrf0.T);Hdc=.5*(Hdc0+Hdc0.T)
q=float(out['physical_parameters']['q_effective']);Vrf=10.;Vdc=float(out['physical_parameters']['Vdc_each_ring_V'])
R0=1e-6;Q=elementary_charge;M=39.96259098*atomic_mass;omega=float(out['physical_parameters']['omega_rad_s'])
gamma=4*Q/(M*R0**2*omega**2)
def rhs(t,z):
    X=z.reshape(6,6);A=np.zeros((6,6));A[:3,3:]=np.eye(3)
    A[3:,:3]=-gamma*(Vrf*Hrf*np.cos(2*t)+Vdc*Hdc)
    return (A@X).ravel()
mon=solve_ivp(rhs,(0,np.pi),np.eye(6).ravel(),method='DOP853',rtol=2e-12,atol=2e-14).y[:,-1].reshape(6,6)
ev=eigvals(mon);betas=sorted(set(round(float(abs(np.angle(x))/np.pi),10) for x in ev if abs(np.angle(x))/np.pi<.5))
freq=[b*omega/(4*np.pi) for b in betas]
g_rf=np.asarray(rfrow['gradient_at_center'],float);g_dc=np.asarray(dcrow['gradient_at_center'],float)
rf_null_xy=-np.linalg.solve(Hrf[:2,:2],g_rf[:2])
h_eff=(Q**2*Vrf**2/(2*M*omega**2*R0**2))*(Hrf.T@Hrf)+Q*Vdc*Hdc
g_eff=(Q**2*Vrf**2/(2*M*omega**2*R0**2))*(Hrf.T@g_rf)+Q*Vdc*g_dc
eq=-np.linalg.solve(h_eff,g_eff)
local={'Hrf_symmetric':Hrf.tolist(),'Hdc_symmetric':Hdc.tolist(),
 'Hrf_antisymmetric_relative_norm':float(np.linalg.norm(Hrf0-Hrf0.T)/np.linalg.norm(Hrf0)),
 'Hdc_antisymmetric_relative_norm':float(np.linalg.norm(Hdc0-Hdc0.T)/np.linalg.norm(Hdc0)),
 'rf_transverse_null_xi_from_local_quadratic':rf_null_xy.tolist(),
 'rf_axial_gradient_residual_per_V':float(g_rf[2]),
 'pseudopotential_equilibrium_xi_from_local_quadratic':eq.tolist(),
 'excess_micromotion_peak_estimate_m':float(.5*q*np.linalg.norm(eq[:2])*R0),
 'monodromy_matrix':mon.tolist(),'multipliers_real_imag':[[float(x.real),float(x.imag)] for x in ev],
 'floquet_beta_principal':betas,'floquet_frequencies_Hz':freq,
 'stability_check_max_abs_multiplier_modulus_error':float(np.max(np.abs(np.abs(ev)-1.)))}
out['local_modes_final']=local
p.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(local,indent=2))
