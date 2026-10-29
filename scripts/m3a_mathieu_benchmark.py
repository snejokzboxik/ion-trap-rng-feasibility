"""Independent integrator check on the ideal q=0.2 Mathieu equation."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
ROOT=Path(__file__).resolve().parents[1]
q=.2; period=np.pi
def A(t): return np.array([[0.,1.],[-2*q*np.cos(2*t),0.]])
def rhs(t,z): return (A(t)@z.reshape(2,2)).ravel()
ref=solve_ivp(rhs,(0,period),np.eye(2).ravel(),method='DOP853',rtol=2e-13,atol=2e-15).y[:,-1].reshape(2,2)
rows=[]
for n in (20,40,80,160):
    dt=period/n; M=np.eye(2);t=0.
    for _ in range(n):
        k1=A(t)@M;k2=A(t+dt/2)@(M+dt*k1/2);k3=A(t+dt/2)@(M+dt*k2/2);k4=A(t+dt)@(M+dt*k3)
        M+=dt*(k1+2*k2+2*k3+k4)/6;t+=dt
    ev=eigvals(M); beta=float(np.min(np.abs(np.angle(ev)))/np.pi)
    evref=eigvals(ref); bref=float(np.min(np.abs(np.angle(evref)))/np.pi)
    rows.append({'steps_per_rf_period':n,'monodromy_relative_frobenius_error':float(np.linalg.norm(M-ref)/np.linalg.norm(ref)),
        'floquet_beta':beta,'beta_error':abs(beta-bref),'determinant':float(np.linalg.det(M)),
        'max_multiplier_modulus_error':float(np.max(np.abs(np.abs(ev)-1.)) )})
out={'equation':'x_ddot + 2 q cos(2 tau) x = 0; q=0.2; RF period pi in tau',
    'reference_method':'DOP853 rtol=2e-13, atol=2e-15 over one period',
    'reference_monodromy':ref.tolist(),'reference_beta':bref,'results':rows}
(ROOT/'results'/'m3a_mathieu_benchmark.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
