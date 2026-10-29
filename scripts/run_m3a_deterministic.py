"""M3a bounded deterministic Ca+ trajectories using saved M2b.2 P2 fields."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge, atomic_mass
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField

R0=1e-6; Q=elementary_charge; M=39.96259098*atomic_mass
VRF=10.0; Q_TARGET=.2; FZ_TARGET=1e6
FIELD='results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz'
FIELD_JSON='results/m2b2_p2_effective_h0.75_hc0.18_dc0.075.json'
OUT=ROOT/'results'

def main():
    meta=json.loads((ROOT/FIELD_JSON).read_text(encoding='utf-8'))
    if meta['geometry'] != 'M2b ideal-conductor hyperbolic branches + annular endcaps; vacuum; polygonal surfaces; analyst-selected supports omitted':
        raise RuntimeError('saved field identity does not match the approved M2b.2 ideal-vacuum model')
    meshmeta=meta['mesh']
    if (meshmeta['P2_dofs']!=215157 or meshmeta['dc_surface_h_r0']!=.075 or
            meshmeta['hcenter_r0']!=.18 or meshmeta['rf_surface_h_r0']!=.3 or meshmeta['axial_end_h_r0']!=.3):
        raise RuntimeError('M2b.2 geometry/mesh identity mismatch')
    rf=P2ScalarField(ROOT/FIELD,'rf_potential')
    dc=P2ScalarField(ROOT/FIELD,'dc_potential')
    rfrow=meta['fields']['rf'][1]; dcrow=meta['fields']['dc'][1]
    Hrf_raw=np.asarray(rfrow['Hessian'],float); Hdc_raw=np.asarray(dcrow['Hessian'],float)
    Hrf=.5*(Hrf_raw+Hrf_raw.T); Hdc=.5*(Hdc_raw+Hdc_raw.T)
    # RF basis field has 1 V differential between the opposing electrode pairs;
    # DC basis is 1 V on each of the two ring terminals.
    k_rf=float(Hrf[0,0]); k_z=float(Hdc[2,2])
    omega=np.sqrt(2*Q*VRF*k_rf/(M*R0**2*Q_TARGET)); f_rf=omega/(2*np.pi)
    vdc=(M*(2*np.pi*FZ_TARGET)**2*R0**2)/(Q*k_z)
    q_eff=2*Q*VRF*k_rf/(M*R0**2*omega**2)
    gamma=4*Q/(M*R0**2*omega**2)
    g_rf=np.asarray(rfrow['gradient_at_center'],float);g_dc=np.asarray(dcrow['gradient_at_center'],float)
    rf_null=-np.linalg.solve(Hrf,g_rf)
    h_eff=(Q**2*VRF**2/(2*M*omega**2*R0**2))*(Hrf.T@Hrf)+Q*vdc*Hdc
    g_eff=(Q**2*VRF**2/(2*M*omega**2*R0**2))*(Hrf.T@g_rf)+Q*vdc*g_dc
    pseudo_equilibrium=-np.linalg.solve(h_eff,g_eff)
    excess_micro=.5*q_eff*np.linalg.norm(pseudo_equilibrium[:2])*R0

    def grad_total(x,tau,phase=0.):
        # P2 potential and gradient are represented in coordinates xi=r/r0.
        g=VRF*np.cos(2*tau+phase)*rf.gradient(x)[...,0]+vdc*dc.gradient(x)[...,0]
        return -gamma*g
    def pot_total(x,tau,phase=0.):
        return VRF*np.cos(2*tau+phase)*rf.potential(x)[0]+vdc*dc.potential(x)[0]
    def rhs(tau,y,phase=0.):
        x=y[:3,None]
        acc=grad_total(x,tau,phase)
        workdot=-2*gamma*VRF*rf.potential(x)[0]*np.sin(2*tau+phase)
        return np.r_[y[3:6],acc,workdot]
    def energy(y,tau,phase=0.):
        return .5*np.dot(y[3:6],y[3:6])+gamma*pot_total(y[:3,None],tau,phase)

    # Poincare/Floquet linearization of the local recovered Hessians.
    def mat_rhs(tau,z,Hrf_use=Hrf,Hdc_use=Hdc):
        X=z.reshape(6,6)
        A=np.zeros((6,6)); A[:3,3:]=np.eye(3)
        A[3:,:3]=-gamma*(VRF*Hrf_use*np.cos(2*tau)+vdc*Hdc_use)
        return (A@X).ravel()
    mon=solve_ivp(mat_rhs,(0,np.pi),np.eye(6).ravel(),method='DOP853',rtol=2e-12,atol=2e-14).y[:,-1].reshape(6,6)
    mult=eigvals(mon); beta=np.sort(np.abs(np.angle(mult))/(np.pi))
    beta_pos=sorted(set(round(float(v),9) for v in beta if v>1e-7 and v<.999999))
    floquet_freq=[float(b*omega/(4*np.pi)) for b in beta_pos]

    # Ideal Mathieu benchmark with a=0 and |q| selected above.
    gamma_i=gamma
    def ideal_rhs(tau,z):
        X=z.reshape(2,2)
        A=np.array([[0.,1.],[-2*Q_TARGET*np.cos(2*tau),0.]])
        return (A@X).ravel()
    imon=solve_ivp(ideal_rhs,(0,np.pi),np.eye(2).ravel(),method='DOP853',rtol=2e-13,atol=2e-15).y[:,-1].reshape(2,2)
    imult=eigvals(imon); ibeta=float(np.min(np.abs(np.angle(imult)))/np.pi)
    f_ideal=ibeta*omega/(4*np.pi)

    # Short selected 3D trajectories: 20 RF cycles; run both RK4 and DOP853.
    periods=20; tau_end=periods*np.pi; dt_ref=np.pi/80
    ics={'near_center':[.02,0,0],'moderate_radial':[.08,0,0],
         'axial':[0,0,.10],'combined':[.05,.04,.08]}
    phases=[0.,np.pi/2]
    records=[]; traces=[]
    for name,x0 in ics.items():
      for phase in phases:
        y0=np.r_[np.asarray(x0,float),np.zeros(3),0.]
        e0=energy(y0,0.,phase)
        for nsub in (40,80,160):
          if nsub!=80 and not (name=='near_center' and phase==0.):
            continue
          dt=np.pi/nsub; n=int(round(tau_end/dt)); y=y0.copy(); tau=0.; maxr=float(np.linalg.norm(y[:3]))
          def f(t,u): return rhs(t,u,phase)
          for _ in range(n):
            k1=f(tau,y); k2=f(tau+dt/2,y+dt*k1/2)
            k3=f(tau+dt/2,y+dt*k2/2); k4=f(tau+dt,y+dt*k3)
            y=y+dt*(k1+2*k2+2*k3+k4)/6; tau+=dt
            maxr=max(maxr,float(np.linalg.norm(y[:3])))
            if maxr>.30:
                break
          residual=energy(y,tau,phase)-e0-y[6]
          records.append({'case':name,'phase_rad':phase,'method':'RK4','steps_per_rf_period':nsub,
              'completed_rf_periods':tau/np.pi,'final_xi':y[:3].tolist(),'max_norm_xi':maxr,
              'energy_balance_residual_scaled':float(residual),'initial_energy_scaled':float(e0),
              'work_rf_scaled':float(y[6]),'escaped_validated_sphere':bool(maxr>.30)})
          if nsub==80 and phase==0.:
            sol=solve_ivp(f,(0,tau_end),y0,method='DOP853',rtol=2e-10,atol=2e-12,max_step=dt_ref,
                          t_eval=np.linspace(0,tau_end,401))
            yf=sol.y[:,-1]; e_res=energy(yf,tau_end,phase)-e0-yf[6]
            records.append({'case':name,'phase_rad':phase,'method':'DOP853','steps_per_rf_period_max':80,
                'function_evaluations':sol.nfev,'completed_rf_periods':tau_end/np.pi,'final_xi':yf[:3].tolist(),
                'max_norm_xi':float(np.max(np.linalg.norm(sol.y[:3],axis=0))),
                'energy_balance_residual_scaled':float(e_res),'initial_energy_scaled':float(e0),
                'work_rf_scaled':float(yf[6]),'escaped_validated_sphere':bool(np.max(np.linalg.norm(sol.y[:3],axis=0))>.30)})
            traces.append({'case':name,'phase_rad':phase,'tauf':sol.t.tolist(),'xi':sol.y[:3].T.tolist()})

    # Difference between integrators/time steps for all sampled conditions.
    comparisons=[]
    for name in ics:
      for phase in (0.,):
        rows=[r for r in records if r['case']==name and r['phase_rad']==phase]
        ref=np.asarray(next(r['final_xi'] for r in rows if r['method']=='DOP853'))
        for r in rows:
          comparisons.append({'case':name,'phase_rad':phase,'method':r['method'],
             'steps_per_rf_period':r.get('steps_per_rf_period',r.get('steps_per_rf_period_max')),
             'final_position_error_r0_vs_DOP853':float(np.linalg.norm(np.asarray(r['final_xi'])-ref))})

    # Estimate RF-period oscillation around the mean position in each RF cycle.
    micromotion=[]
    for tr in traces:
        x=np.asarray(tr['xi']); groups=np.array_split(np.arange(len(x)),periods)
        amps=[]
        for ids in groups:
            cyc=x[ids]; amps.append(np.sqrt(np.mean(np.sum((cyc-cyc.mean(axis=0))**2,axis=1))))
        micromotion.append({'case':tr['case'],'phase_rad':tr['phase_rad'],
            'rms_within_rf_period_r0_median':float(np.median(amps)),
            'rms_within_rf_period_r0_max':float(np.max(amps)),
            'interpretation':'includes secular drift over the 20 samples/RF period; indicative micromotion scale, not a spectral separation'})

    # Local pseudopotential acceptance on a validated-radius sphere.
    nmu,nphi=15,32; mu=np.linspace(-1,1,nmu); phi=np.arange(nphi)*2*np.pi/nphi
    dirs=[]
    for mz in mu:
      for ph in phi:
        sr=np.sqrt(1-mz*mz);dirs.append([sr*np.cos(ph),sr*np.sin(ph),mz])
    dirs=np.asarray(dirs).T; pp=.20*dirs
    U=[]
    for i in range(pp.shape[1]):
      x=pp[:,i:i+1]; Eamp=-VRF*rf.gradient(x)[:,0]/R0
      us=Q**2*np.dot(Eamp,Eamp)/(4*M*omega**2)+Q*vdc*dc.potential(x)[0]
      U.append(us)
    U=np.asarray(U)
    g0=rf.gradient(np.zeros((3,1)))[:,0]/R0
    U0=Q**2*VRF**2*np.dot(g0,g0)/(4*M*omega**2)+Q*vdc*dc.potential(np.zeros((3,1)))[0]
    safe=float(np.min(U-U0)); idx=int(np.argmin(U-U0))
    out={'stage':'M3a','model':'deterministic single 40Ca+; saved M2b.2 P2 fields; ideal PEC vacuum; no stochastic forces',
      'field_identity':{'field_npz':FIELD,'field_json':FIELD_JSON,'geometry':meta['geometry'],'mesh':meshmeta,
        'normalization':{'rf_basis':'1 V differential between opposite RF electrode pairs','dc_basis':'1 V on each DC ring','potential_coordinate':'r/r0','physical_gradient':'(1/r0)*dPhi/d(r/r0)'}},
      'physical_parameters':{'r0_m':R0,'mass_kg':M,'charge_C':Q,'Vrf_differential_peak_V':VRF,'target_q_abs':Q_TARGET,
        'rf_coefficient_Hxx_per_V':k_rf,'omega_rad_s':omega,'f_rf_Hz':f_rf,'q_effective':q_eff,
        'dc_curvature_Hzz_per_V_V_m2':k_z/R0**2,'Vdc_each_ring_V':vdc,'target_fz_Hz':FZ_TARGET,
        'gamma_dimensionless':gamma},
      'local_modes':{'Hrf_normalized_symmetric':Hrf.tolist(),'Hdc_per_V_normalized_symmetric':Hdc.tolist(),
        'Hrf_antisymmetric_norm':float(np.linalg.norm(Hrf_raw-Hrf_raw.T)/np.linalg.norm(Hrf_raw)),
        'Hdc_antisymmetric_norm':float(np.linalg.norm(Hdc_raw-Hdc_raw.T)/np.linalg.norm(Hdc_raw)),
        'rf_null_estimate_xi':rf_null.tolist(),'pseudopotential_equilibrium_xi':pseudo_equilibrium.tolist(),
        'excess_micromotion_peak_estimate_m':float(excess_micro),
        'monodromy_multipliers_real_imag':[[float(v.real),float(v.imag)] for v in mult],
        'finite_field_floquet_frequencies_Hz':floquet_freq,'ideal_a0_mathieu_q02_frequency_Hz':f_ideal,
        'small_q_approx_radial_frequency_Hz':float(.2/ (2*np.sqrt(2))*f_rf)},
      'pseudopotential_local':{'sphere_radius_r0':.20,'min_delta_U_J_vs_center':safe,'min_delta_U_eV':safe/Q,'escape_direction':dirs[:,idx].tolist(),
        'interpretation':'safe-radius local acceptance energy only; not the actual finite-geometry escape barrier or trap depth'},
      'trajectory_runs':records,'integrator_comparison':comparisons,'micromotion_estimates':micromotion,'trajectory_trace':traces,
      'limitations':['selected trajectories stop at radius 0.30r0 and duration 20 RF periods; axial 1 MHz motion is not directly period-resolved',
        'P2 gradients are elementwise discontinuous; mesh-to-mesh pointwise DC-gradient differences are significant relative to the weak DC field near center',
        'only DC surface size varied among the field meshes; center h=0.18r0 remains fixed',
        'electrode escape barriers, end fringing, curved boundaries, material interfaces and RF circuit are not validated',
        'linear Floquet stability does not prove nonlinear long-term confinement']}
    (OUT/'m3a_deterministic_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    np.savez_compressed(OUT/'m3a_trajectory_traces.npz',traces=np.array(traces,dtype=object))
    print(json.dumps({k:v for k,v in out.items() if k not in ('trajectory_runs','trajectory_trace')},indent=2))
    print('Trajectory records:',len(records),'integrator comparisons:',len(comparisons))

if __name__=='__main__':main()
